"""
    prepare(c::Campaign)

Phase 1 (single process): enumerate centers/pairs/batches and write
`centers.tsv` + `manifest.tsv` under `c.root`. Idempotent: existing files are kept
(delete them to re-prepare).
"""
_params_string(c::Campaign) =
    "nx=$(c.nx) ny=$(c.ny) cutoff=$(c.neighbor_cutoff) n_per_batch=$(c.n_centers_per_batch)"

function prepare(c::Campaign)
    params_file = joinpath(c.root, "manifest.params")
    if isfile(manifest_path(c)) && isfile(centers_path(c))
        @info "prepare: manifest exists, skipping" manifest_path(c)
        if isfile(params_file)
            stored = strip(read(params_file, String))
            current = _params_string(c)
            stored == current || error(
                "Campaign parameters changed since last prepare.\n" *
                "  stored : $stored\n" *
                "  current: $current\n" *
                "Delete the manifest files in $(c.root) to re-prepare.")
        else
            @warn "prepare: manifest.params missing (old campaign); skipping parameter check"
        end
        return read_manifest(manifest_path(c))
    end
    temps = load_templates!(c)
    st1 = temps[1][1]
    centroids = [BoundaryIntegral.density_centroid(t[2]) for t in temps]
    s1 = lattice_grid_steps(temps[1][2], st1.primvec, (1, 0, 0))
    s2 = lattice_grid_steps(temps[1][2], st1.primvec, (0, 1, 0))
    centers = enumerate_centers(c.nx, c.ny, st1.primvec,
        [ntuple(d -> Float64(ct[d]), 3) for ct in centroids], (s1, s2))
    pairs = enumerate_pairs(centers, c.neighbor_cutoff)
    batches = build_batches(pairs, c.n_centers_per_batch)
    write_centers(centers_path(c), centers)
    write_manifest(manifest_path(c), batches)
    mkpath(dirname(params_file))
    write(params_file, _params_string(c))
    @info "prepare: wrote manifest" n_centers=length(centers) n_pairs=length(pairs) n_batches=length(batches)
    return batches
end

"""
    pending_batches(c, phase::Symbol) -> Vector{Int}

Batch ids still to do, derived from files on disk (spec §3: no mutable status).
`:solve` → no complete batch file; `:eval` → no complete V file.
"""
# presence == completeness: batch/V files are written atomically (tmp+rename); deep validation happens at actual load sites.
function pending_batches(c::Campaign, phase::Symbol)
    isfile(manifest_path(c)) || error("manifest not found at $(manifest_path(c)); run prepare(c) first")
    batches = read_manifest(manifest_path(c))
    ids = getfield.(batches, :batch_id)
    if phase === :solve
        return [id for id in ids if !isfile(batch_path(c, id))]
    elseif phase === :eval
        return [id for id in ids if !_is_complete_v(v_path(c, id))]
    end
    error("unknown phase $phase")
end

# Temporary stub — Task 12 replaces this with a validating loader for V files.
_is_complete_v(path::AbstractString) = isfile(path)

# OrbitalInstances for the centers referenced by a batch (template grids shared via cache)
function _batch_instances(c::Campaign, spec::BatchSpec)
    centers = read_centers(centers_path(c))
    byid = Dict(ct.id => ct for ct in centers)
    need = unique(reduce(vcat, [[p[1], p[2]] for p in spec.pairs]; init = Int[]))
    return Dict(id => OrbitalInstance(id, byid[id].template_id, byid[id].steps) for id in need)
end

"""
    solve_batch(c::Campaign, batch_id) -> path | nothing

Solve phase for one batch (spec §5): assemble pair densities on the global grid,
envelope-refine ONE shared interface, block-GMRES, write the BatchResult atomically.
Skips (returns nothing) if the output already exists.
"""
function solve_batch(c::Campaign, batch_id::Int)
    out = batch_path(c, batch_id)
    if isfile(out)
        @info "solve_batch: already complete, skipping" batch_id
        return nothing
    end
    t0 = time()
    found = filter(b -> b.batch_id == batch_id, read_manifest(manifest_path(c)))
    isempty(found) && error("batch_id=$batch_id not found in manifest $(manifest_path(c))")
    spec = only(found)
    temps = load_templates!(c)
    grids = [t[2] for t in temps]
    insts = _batch_instances(c, spec)
    t_setup = time() - t0

    t1 = time()
    b = assemble_lattice_batch(grids, insts, spec.pairs;
        support_rtol = c.solve["support_rtol"])
    t_asm = time() - t1

    t2 = time()
    res = solve_dielectric_lattice_batch(c.boxes, c.epses, c.eps_out, b;
        n_quad = Int(c.solve["n_quad"]), rhs_atol = c.solve["rhs_tol"],
        l_ec = campaign_l_ec(c), fmm_tol = c.solve["lhs_tol"],
        up_tol = c.solve["lhs_tol"], max_order = Int(c.solve["max_order"]),
        gmres_rtol = c.solve["gmres_rtol"], max_depth = Int(c.solve["max_depth"]))
    t_solve = time() - t2
    t_total = time() - t0

    stats = Dict{String,Any}(
        "t_setup" => t_setup, "t_assemble" => t_asm, "t_solve" => t_solve, "t_total" => t_total,
        "niter" => res.stats.niter, "dof" => size(res.sigma, 1),
        "n_support" => length(b.gidx), "K" => length(spec.pairs),
        "hostname" => gethostname())
    br = BatchResult(BoundaryIntegral.BATCH_FORMAT_VERSION, batch_id, b.pair_ids,
        b.gidx, b.weights, b.densities, res.interface, res.sigma, stats)
    save_batch_result(out, br)
    @info "solve_batch: done" batch_id dof=stats["dof"] K=stats["K"] t_total
    return out
end

"""
    consolidate(c::Campaign)

Between solve and eval (single process; plan deviation from spec §3 documented there):
build
- `targets.jls`: `(; gidx, positions)` — the sorted union of ALL stored batch supports
  (the shared eval target set T), positions from the template's affine map;
- `rho_store.jls`: `(; pair_ids, t_idx, tw)` — per pair, its support's row indices in T
  and the contraction vector `w .* ρ` (raw density), in manifest batch order.
Built from the SOLVED batch files, so T is exactly consistent with the stored ρ.
Errors if any batch is missing. Idempotent (overwrites deterministically).

Memory note: this loads every BatchResult (σ + interface included) to read only
gidx/weights/densities. For the production campaign (~200 batches × ~100 MB) that is
~20 GB on a 1.5 TB node — acceptable; a leaner loader that drops σ/interface could be
added if memory becomes tight.
"""
function consolidate(c::Campaign)
    batches = read_manifest(manifest_path(c))
    missing_ids = [b.batch_id for b in batches if !isfile(batch_path(c, b.batch_id))]
    isempty(missing_ids) || error("consolidate: unsolved batches: $missing_ids")

    brs = [load_batch_result(batch_path(c, b.batch_id)) for b in batches]
    gall = NTuple{3,Int}[]
    for br in brs
        append!(gall, br.gidx)
    end
    gidx = sort(unique(gall))
    rowofg = Dict(g => r for (r, g) in enumerate(gidx))

    temps = load_templates!(c)
    dg = temps[1][2]
    positions = Matrix{Float64}(undef, 3, length(gidx))
    for (r, g) in enumerate(gidx)
        p = BoundaryIntegral.grid_point(dg, g[1], g[2], g[3])
        positions[1, r] = p[1]; positions[2, r] = p[2]; positions[3, r] = p[3]
    end
    _atomic_serialize(targets_path(c), (; gidx, positions))

    pair_ids = Tuple{Int,Int}[]
    t_idx = Vector{Vector{Int}}()
    tw = Vector{Vector{Float64}}()
    for br in brs
        rows = [rowofg[g] for g in br.gidx]
        for k in 1:length(br.pair_ids)
            push!(pair_ids, br.pair_ids[k])
            push!(t_idx, rows)
            push!(tw, br.weights .* br.densities[:, k])
        end
    end
    _atomic_serialize(rho_store_path(c), (; pair_ids, t_idx, tw))
    @info "consolidate: done" n_targets=length(gidx) n_pairs=length(pair_ids)
    return nothing
end

function _atomic_serialize(path::AbstractString, obj)
    d = dirname(path)
    isempty(d) || mkpath(d)
    tmp = string(path, ".tmp.", getpid(), "_", rand(UInt32))
    open(io -> serialize(io, obj), tmp, "w")
    mv(tmp, path; force = true)
end

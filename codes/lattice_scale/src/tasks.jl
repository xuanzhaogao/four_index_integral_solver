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

const V_FORMAT_VERSION = 1

function save_v_rows(path::AbstractString, batch_id::Int,
        source_pairs::Vector{Tuple{Int,Int}}, target_pairs::Vector{Tuple{Int,Int}},
        V::Matrix{Float64}, stats::Dict{String,Any})
    _atomic_serialize(path, (; version = V_FORMAT_VERSION, batch_id,
        source_pairs, target_pairs, V, stats))
end

function load_v_rows(path::AbstractString)
    vr = open(deserialize, path)
    vr.version == V_FORMAT_VERSION || error("$path: V format version mismatch")
    size(vr.V) == (length(vr.target_pairs), length(vr.source_pairs)) ||
        error("$path: V shape mismatch")
    return vr
end

function _is_complete_v(path::AbstractString)
    isfile(path) || return false
    try
        load_v_rows(path)
        return true
    catch
        return false
    end
end

"""
    eval_batch(c::Campaign, batch_id) -> path | nothing

Post-eval phase for one batch (spec §6): rebuild the K sources from the BatchResult,
evaluate Φ_a = u_inc[ρ_a] + u[σ_a] ONCE at the shared target set T (near/far split in
`evaluate_batch_potential`), contract against every stored pair density, write the
K columns of V atomically. Skips if the V file is already complete.
"""
function eval_batch(c::Campaign, batch_id::Int)
    out = v_path(c, batch_id)
    if _is_complete_v(out)
        @info "eval_batch: already complete, skipping" batch_id
        return nothing
    end
    t0 = time()
    br = load_batch_result(batch_path(c, batch_id))

    isfile(targets_path(c)) && isfile(rho_store_path(c)) ||
        error("eval_batch: targets.jls / rho_store.jls not found under $(c.root); run consolidate(c) first")

    T = open(deserialize, targets_path(c))
    store = open(deserialize, rho_store_path(c))
    temps = load_templates!(c)
    dg = temps[1][2]

    K = length(br.pair_ids)
    pos = Matrix{Float64}(undef, 3, length(br.gidx))
    for (r, g) in enumerate(br.gidx)
        p = BoundaryIntegral.grid_point(dg, g[1], g[2], g[3])
        pos[1, r] = p[1]; pos[2, r] = p[2]; pos[3, r] = p[3]
    end
    sources = [VolumeSource(copy(pos), copy(br.weights), br.densities[:, k]) for k in 1:K]
    At, Bt, Ct = BoundaryIntegral.true_cell_vectors(dg)
    max_step = maximum((norm(At) / dg.nx, norm(Bt) / dg.ny, norm(Ct) / dg.nz))
    far_pad = c.far_pad_steps * max_step

    t_phi_start = time()
    Φ = evaluate_batch_potential(br.interface, br.sigma, sources, T.positions;
        lhs_tol = c.solve["lhs_tol"], volume_tol = c.solve["volume_tol"],
        far_pad = far_pad)
    t_phi = time() - t_phi_start
    t_setup = t_phi_start - t0

    nP = length(store.pair_ids)
    V = Matrix{Float64}(undef, nP, K)
    for kl in 1:nP, a in 1:K
        V[kl, a] = dot(store.tw[kl], view(Φ, store.t_idx[kl], a))
    end
    stats = Dict{String,Any}("t_setup" => t_setup, "t_phi" => t_phi, "t_total" => time() - t0,
        "n_targets" => size(T.positions, 2), "hostname" => gethostname())
    save_v_rows(out, batch_id, br.pair_ids, store.pair_ids, V, stats)
    @info "eval_batch: done" batch_id n_targets=size(T.positions, 2) t_total=stats["t_total"]
    return out
end

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
    gset = Set{NTuple{3,Int}}()
    for br in brs
        union!(gset, br.gidx)
    end
    gidx = sort!(collect(gset))
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
            push!(t_idx, copy(rows))
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

"""
    assemble_v(c::Campaign) -> (; max_rel_asym, n)

Final phase (single process, spec §4/§6): gather all V files into the dense
`n_pairs × n_pairs` matrix (columns keyed by source pair, rows by target pair — the
SAME pair ordering, from the rho store), write `V_full.jls` and `report.txt` with the
campaign's built-in accuracy diagnostic max |V - Vᵀ| / max |V| (each entry is computed
twice, from the two independently adapted interfaces).
"""
function assemble_v(c::Campaign)
    batches = read_manifest(manifest_path(c))
    store = open(deserialize, rho_store_path(c))
    pair_ids = store.pair_ids
    col = Dict(p => i for (i, p) in enumerate(pair_ids))
    n = length(pair_ids)
    V = fill(NaN, n, n)
    for b in batches
        vr = load_v_rows(v_path(c, b.batch_id))
        vr.target_pairs == pair_ids || error("V_$(b.batch_id): target ordering mismatch")
        for (k, sp) in enumerate(vr.source_pairs)
            V[:, col[sp]] = vr.V[:, k]
        end
    end
    any(isnan, V) && error("assemble_v: missing columns (run eval for all batches first)")

    scale = maximum(abs.(V))
    max_rel_asym = maximum(abs.(V .- transpose(V))) / scale
    _atomic_serialize(joinpath(c.root, "V_full.jls"), (; pair_ids, V))
    open(joinpath(c.root, "report.txt"), "w") do io
        println(io, "campaign: $(c.name)")
        println(io, "pairs: $n   batches: $(length(batches))")
        println(io, "max|V|: $scale")
        println(io, "max rel asymmetry |V - V'|/max|V|: $max_rel_asym")
    end
    @info "assemble_v: done" n max_rel_asym
    return (; max_rel_asym, n)
end

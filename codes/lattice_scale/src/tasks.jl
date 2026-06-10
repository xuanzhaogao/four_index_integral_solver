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

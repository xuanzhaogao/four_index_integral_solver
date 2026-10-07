# Onsite-U eval: for each solved batch, evaluate the screened potential Φ_i of the onsite
# density ρ_ii at the orbital's OWN support and contract → U_i = ⟨ρ_ii | Φ_i⟩ (in eV).
# This is the cheap, onsite-ONLY path: 1792 local evals (each at the orbital's own ~few-k
# support), NOT the campaign's full-matrix eval (Φ at the shared union of all supports).
# It REPLACES consolidate/eval/assemble — run it after `prepare` + `solve`.
#
# Resumable: appends one row per orbital to <root>/onsite_U.tsv and skips finished ones.
#
# Run (after prepare+solve):
#   julia --project=codes/article -t8 codes/article/sec53_lattice/scripts/onsite_eval.jl campaigns/onsite_line.toml

using BoundaryIntegral, Serialization, LinearAlgebra
const BI = BoundaryIntegral
const E2 = 14.3996      # e^2/(4πε0) in eV·Å (matches assemble_v's eV normalization)

c       = BI.load_campaign(abspath(ARGS[1]))
temps   = BI.load_templates!(c)
dg      = temps[1][2]                       # shared grid geometry (all templates compatible)
batches = BI.read_manifest(BI.manifest_path(c))
outpath = joinpath(c.root, "onsite_U.tsv")

done = Set{Int}()
if isfile(outpath)
    for ln in eachline(outpath)
        (isempty(ln) || startswith(ln, "orbital")) && continue
        push!(done, parse(Int, first(split(ln, '\t'))))
    end
end

open(outpath, "a") do io
    isempty(done) && println(io, "orbital\tx\ty\tz\tU_eV\tnorm2")
    for b in batches
        isfile(BI.batch_path(c, b.batch_id)) || continue   # skip not-yet-solved batches
        br  = BI.load_batch_result(BI.batch_path(c, b.batch_id))
        ids = unique(p[1] for p in br.pair_ids if p[1] == p[2])
        all(in(done), ids) && continue                       # batch already evaluated

        pos = BI.grid_positions(dg, br.gidx)                 # the batch's OWN support
        At, Bt, Ct = BI.true_cell_vectors(dg)
        max_step = maximum((norm(collect(At)) / dg.nx,
                            norm(collect(Bt)) / dg.ny,
                            norm(collect(Ct)) / dg.nz))
        far_pad = c.far_pad_steps * max_step

        K = length(br.pair_ids)
        sources = [BI.VolumeSource(copy(pos), copy(br.weights), br.densities[:, k]) for k in 1:K]
        Φ = BI.evaluate_batch_potential(br.interface, br.sigma, sources, pos;
                lhs_tol = c.solve["lhs_tol"], volume_tol = c.solve["volume_tol"],
                far_pad = far_pad,
                screen_boxes = c.boxes, screen_epses = c.epses, screen_eps_out = c.eps_out)

        for k in 1:K
            (i, j) = br.pair_ids[k]
            (i == j && !(i in done)) || continue              # onsite only
            wρ  = br.weights .* br.densities[:, k]
            n2  = sum(wρ)                                      # ∫ρ_ii = ‖φ_i‖²
            U   = dot(wρ, view(Φ, :, k)) * 4π * E2 / n2^2      # onsite U_i in eV
            p   = c.orbitals[i].pos
            println(io, i, '\t', p[1], '\t', p[2], '\t', p[3], '\t', U, '\t', n2)
            flush(io); push!(done, i)
        end
        @info "onsite_eval: batch done" batch_id=b.batch_id n_done=length(done)
    end
end
println("onsite_eval: $(length(done)) onsite U → $(outpath)")

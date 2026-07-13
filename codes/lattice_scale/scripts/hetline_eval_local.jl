# Per-orbital LOCAL onsite-U eval for a BATCHED het-line solve. The stock onsite_eval.jl
# evaluates Φ over the batch's full union support in one call — for a wide batched line that
# union spans ~244 Å and TKM3D's FINUFFT fails to allocate. Here we reuse the cached batched
# σ but evaluate each orbital at ITS OWN support only (small extent → no FINUFFT blowup), exactly
# like the n_centers_per_batch=1 path. Onsite U_i in eV. Resumable-overwrite (rewrites the TSV).
#
# Run (after prepare+solve):  julia --project -t <n> scripts/hetline_eval_local.jl <campaign.toml>
using BoundaryIntegral, LinearAlgebra
const BI = BoundaryIntegral
const E2 = 14.3996

c   = BI.load_campaign(abspath(ARGS[1]))
dg  = BI.load_templates!(c)[1][2]
br  = BI.load_batch_result(BI.batch_path(c, 1))
pos = BI.grid_positions(dg, br.gidx)
At, Bt, Ct = BI.true_cell_vectors(dg)
far_pad = c.far_pad_steps * maximum((norm(collect(At)) / dg.nx,
                                     norm(collect(Bt)) / dg.ny,
                                     norm(collect(Ct)) / dg.nz))
outpath = joinpath(c.root, "onsite_U.tsv")

open(outpath, "w") do io
    println(io, "orbital\tx\ty\tz\tU_eV\tnorm2")
    K = length(br.pair_ids)
    for k in 1:K
        (i, j) = br.pair_ids[k]
        i == j || continue
        supp = findall(!iszero, view(br.densities, :, k))     # this orbital's own support
        posk = pos[:, supp]; wk = br.weights[supp]; dk = br.densities[supp, k]
        src = BI.VolumeSource(copy(posk), copy(wk), copy(dk))
        Φ = BI.evaluate_batch_potential(br.interface, br.sigma[:, k:k], [src], posk;
                lhs_tol = c.solve["lhs_tol"], volume_tol = c.solve["volume_tol"], far_pad = far_pad,
                screen_boxes = c.boxes, screen_epses = c.epses, screen_eps_out = c.eps_out)
        wρ = wk .* dk; n2 = sum(wρ)
        U = dot(wρ, view(Φ, :, 1)) * 4π * E2 / n2^2
        p = c.orbitals[i].pos
        println(io, i, '\t', p[1], '\t', p[2], '\t', p[3], '\t', U, '\t', n2); flush(io)
        @info "local eval" orbital=i x=p[1] U_eV=U
    end
end
println("hetline_eval_local: → $(outpath)")

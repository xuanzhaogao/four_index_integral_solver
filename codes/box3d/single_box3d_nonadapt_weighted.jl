include(joinpath(@__DIR__, "single_box3d_utils.jl"))

function single_box3d_nonadpt_weighted()

    ref = load(joinpath(@__DIR__, "data/single_box3d_ref_non_reduce.jld2"))
    trgs = BI.sphere_points(2.0, 100, 100)
    S_ref = BI.laplace3d_pottarg_fmm3d(ref["dbox"], trgs, 1e-6)
    pot_ref = S_ref * ref["sigma"]

    df = joinpath(@__DIR__, "data/single_box3d_nonadapt_weighted.csv")
    CSV.write(df, DataFrame(n_quad = [], n_val = [], n_iter = [], flux = [], pot_abserr = [], pot_relerr = [], L1_err = [], L2_err = [], Linf_err = []))

    n_boxes = 1
    eps = 4.0
    reduce_quad = false
    n_edge = 0

    for n_quad in [1] ∪ [4:4:72...]

        dbox, sigma, flux, n_val, n_iter = solve_single_box3d_weighted(eps, n_boxes, n_quad, reduce_quad, n_edge, n_edge, (0.2, 0.3, 0.4))

        S_trg = BI.laplace3d_pottarg_fmm3d(dbox, trgs, 1e-6)
        pot_trg = S_trg * sigma
        pot_abserr = maximum(abs.(pot_trg - pot_ref))
        pot_relerr = maximum(abs.((pot_trg - pot_ref) ./ pot_ref))

        D_fmm3d = BI.laplace3d_D_fmm3d(dbox, 1e-6)
        o = ones(size(D_fmm3d, 1))
        w = BI.all_weights(dbox)
        gi = abs.(D_fmm3d * o .+ 0.5)
        L1_err = sum(w .* (gi .+ 0.5))
        L2_err = sqrt(sum(w .* (gi .+ 0.5).^2))
        Linf_err = norm(gi, Inf)

        println("n_quad = $n_quad, n_val = $n_val, n_iter = $n_iter, flux = $flux, pot_abserr = $pot_abserr, pot_relerr = $pot_relerr, L1_err = $L1_err, L2_err = $L2_err, Linf_err = $Linf_err")

        CSV.write(df, DataFrame(n_quad = [n_quad], n_val = [n_val], n_iter = [n_iter], flux = [flux], pot_abserr = [pot_abserr], pot_relerr = [pot_relerr], L1_err = [L1_err], L2_err = [L2_err], Linf_err = [Linf_err]), append = true)
    end
end

single_box3d_nonadpt_weighted()
# single_box2d_nonadpt()
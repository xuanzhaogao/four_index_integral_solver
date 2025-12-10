include(joinpath(@__DIR__, "single_box3d_utils.jl"))

function single_box3d_nonadpt()

    ref = load(joinpath(@__DIR__, "data/single_box3d_ref_non_reduce.jld2"))
    trgs = BI.sphere_points(2.0, 100, 100)
    S_ref = BI.laplace3d_pottarg_fmm3d(ref["dbox"], trgs, 1e-6)
    pot_ref = S_ref * ref["sigma"]

    df = joinpath(@__DIR__, "data/single_box3d_nonadapt.csv")
    CSV.write(df, DataFrame(n_quad = [], n_val = [], n_iter = [], flux = [], pot_abserr = [], pot_relerr = [], L1_err = [], L2_err = [], Linf_err = []))

    n_boxes = 1
    eps = 4.0
    reduce_quad = false
    n_edge = 0

    for n_quad in [1] ∪ [4:4:72...]

        dbox, sigma, flux, n_val, n_iter = solve_single_box3d(eps, n_boxes, n_quad, reduce_quad, n_edge, n_edge, (0.2, 0.3, 0.4))

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

function single_box2d_nonadpt()
    df2d = joinpath(@__DIR__, "data/single_box2d_nonadapt.csv")
    CSV.write(df2d, DataFrame(n_quad = [], L1_err = [], L2_err = [], Linf_err = [], pt_flux = []))
    for n_quad in [1, 4, 8, 12, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
        dbox = BI.dielectric_box2d(1, n_quad)
        D_fmm2d = BI.laplace2d_D(dbox)
        o = ones(size(D_fmm2d, 1))
        gi = abs.(D_fmm2d * o .- 0.5)

        L1_err = norm(gi, 1)
        L2_err = norm(gi, 2)
        Linf_err = norm(gi, Inf)

        pt_flux = BI.l2d_point_gi(dbox, (0.2, 0.3))

        @show n_quad, L1_err, L2_err, Linf_err, pt_flux
        CSV.write(df2d, DataFrame(n_quad = [n_quad], L1_err = [L1_err], L2_err = [L2_err], Linf_err = [Linf_err], pt_flux = [pt_flux]), append = true)
    end
end

single_box3d_nonadpt()
# single_box2d_nonadpt()
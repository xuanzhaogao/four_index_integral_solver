include(joinpath(@__DIR__, "single_box3d_utils.jl"))

function main()
    res_ref = load(joinpath(@__DIR__, "data/single_box3d_ref_non_reduce.jld2"))
    pot_ref = res_ref["pot_trg"]
    trgs = load(joinpath(@__DIR__, "data/sphere_trgs.jld"))["trgs"]

    df = joinpath(@__DIR__, "data/single_box3d_reduce_convergence.csv")
    CSV.write(df, DataFrame(n_quad_max = [], n_quad_min = [], n_edge = [], gi = [], n_val = [], n_iter = [], pot_abserr = [], pot_relerr = []))

    n_quad_maxs = [4, 6, 8, 10, 12, 14]
    n_quad_mins = [2, 4]
    n_edges = 0:2:10
    src = (0.2, 0.3, 0.4)

    for n_quad_max in n_quad_maxs
        for n_quad_min in n_quad_mins
            for n_edge in n_edges
                dbox, sigma, gi, n_val, n_iter = solve_single_box3d(4.0, 1, n_quad_max, n_quad_min, n_edge, n_edge, (0.2, 0.3, 0.4))
                S_trg = BI.laplace3d_pottarg_fmm3d(dbox, trgs, 1e-6)
                pot_trg = S_trg * sigma
                pot_abserr = maximum(abs.(pot_trg - pot_ref))
                pot_relerr = maximum(abs.((pot_trg - pot_ref) ./ pot_ref))

                println("n_quad_max = $n_quad_max, n_quad_min = $n_quad_min, n_edge = $n_edge, gi = $gi, n_val = $n_val, n_iter = $n_iter, pot_abserr = $pot_abserr, pot_relerr = $pot_relerr")

                CSV.write(df, DataFrame(n_quad_max = [n_quad_max], n_quad_min = [n_quad_min], n_edge = [n_edge], gi = [gi], n_val = [n_val], n_iter = [n_iter], pot_abserr = [pot_abserr], pot_relerr = [pot_relerr]), append = true)
            end
        end
    end
end

main()
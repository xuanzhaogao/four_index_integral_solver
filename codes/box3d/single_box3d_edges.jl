# simply do edge refinement for all surfaces
include(joinpath(@__DIR__, "single_box3d_utils.jl"))

function main()

    res_ref = load(joinpath(@__DIR__, "data/single_box3d_ref_non_reduce.jld2"))
    pot_ref = res_ref["pot_trg"]
    trgs = load(joinpath(@__DIR__, "data/sphere_trgs.jld"))["trgs"]

    df = joinpath(@__DIR__, "data/single_box3d_edges.csv")
    CSV.write(df, DataFrame(n_boxes = [], n_quad = [], reduce_quad = [], n_edge = [], n_corner = [], gi = [], n_val = [], n_iter = [], pot_abserr = [], pot_relerr = []))

    n_boxes = 1
    eps = 4.0

    for reduce_quad in [false]
        for n_edge in 0:2:8
            for n_quad in [1, 4, 8, 12, 16]
                dbox, sigma, gi, n_val, n_iter = solve_single_box3d(eps, n_boxes, n_quad, reduce_quad, n_edge, n_edge, (0.2, 0.3, 0.4))
                S_trg = BI.laplace3d_pottarg_fmm3d(dbox, trgs, 1e-6)
                pot_trg = S_trg * sigma
                pot_abserr = maximum(abs.(pot_trg - pot_ref))
                pot_relerr = maximum(abs.((pot_trg - pot_ref) ./ pot_ref))

                println("n_quad = $n_quad, reduce_quad = $reduce_quad, n_edge = $n_edge, gi = $gi, n_val = $n_val, n_iter = $n_iter, pot_abserr = $pot_abserr, pot_relerr = $pot_relerr")

                r_flag = reduce_quad ? 1 : 0
                CSV.write(df, DataFrame(n_boxes = [n_boxes], n_quad = [n_quad], reduce_quad = [r_flag], n_edge = [n_edge], n_corner = [n_edge], gi = [gi], n_val = [n_val], n_iter = [n_iter], pot_abserr = [pot_abserr], pot_relerr = [pot_relerr]), append = true)
            end
        end
    end
end

main()
include(joinpath(@__DIR__, "single_box3d_utils.jl"))

function main()

    df = joinpath(@__DIR__, "data/single_box3d_gi.csv")
    CSV.write(df, DataFrame(n_boxes = [], n_quad = [], reduce_quad = [], n_edge = [], n_corner = [], L1_err = [], L2_err = [], Linf_err = []))

    n_boxes = 1
    eps = 4.0

    for reduce_quad in [false]
        for n_edge in 0:2:12
            for n_quad in [1, 4, 8, 12, 16]
                dbox = BI.dielectric_box3d(eps, 1.0, 1, n_quad, n_quad, n_edge, n_edge)
                D = BI.laplace3d_D_fmm3d(dbox, 1e-6)
                o = ones(size(D, 1))
                gi = abs.(D * o .+ 0.5)

                L1_err = norm(gi, 1)
                L2_err = norm(gi, 2)
                Linf_err = norm(gi, Inf)

                @show n_quad, reduce_quad, n_edge, L1_err, L2_err, Linf_err

                CSV.write(df, DataFrame(n_boxes = [n_boxes], n_quad = [n_quad], reduce_quad = [reduce_quad], n_edge = [n_edge], n_corner = [n_edge], L1_err = [L1_err], L2_err = [L2_err], Linf_err = [Linf_err]), append = true)
            end
        end
    end

end
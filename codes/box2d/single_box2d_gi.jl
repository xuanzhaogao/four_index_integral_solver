using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
using CSV, DataFrames

function gi_quad_points()
    rect = BI.square(-0.5, -0.5)
    df = CSV.write("data/single_box2d_gi.csv", DataFrame(n_quad = [], n_adapt = [], L1_err = [], L2_err = [], Linf_err = []))

    for n_quad in [1, 2, 4, 8, 16, 32]
        for n_adapt in [1, 2, 4, 8, 16, 32]
            dbox = BI.dielectric_mbox2d([4.0], [rect], 2, n_quad, n_adapt)
            D_fmm2d = BI.laplace2d_D_fmm2d(dbox, 1e-12)

            gi = D_fmm2d * ones(BI.num_points(dbox))

            w = BI.all_weights(dbox)
            L1_err = sum(w .* abs.(gi .+ 0.5))
            L2_err = sqrt(sum(w .* (gi .+ 0.5).^2))
            Linf_err = norm(gi .+ 0.5, Inf)

            @show n_quad, n_adapt, L1_err, L2_err, Linf_err

            CSV.write("data/single_box2d_gi.csv", DataFrame(n_quad = n_quad, n_adapt = n_adapt, L1_err = L1_err, L2_err = L2_err, Linf_err = Linf_err), append = true)
        end
    end
end


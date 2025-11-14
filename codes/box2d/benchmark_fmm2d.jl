# benchmark the performance of the BIE solver with FMM2D.jl, including the time for applying the linear operator and the time for the GMRES solver

using BoundaryIntegral
import BoundaryIntegral as BI
using Krylov
using CSV, DataFrames
using LinearAlgebra, Random
using BenchmarkTools

function benchmark()
    rects = [BI.square(-1.25, -0.25), BI.square(-0.25, 0.25), BI.square(0.25, -0.75), BI.square(-0.75, -1.25)]
    eps_boxes = [2.0, 3.0, 4.0, 5.0]
    r_src = (0.1, 0.2)
    eps_src = 1.0

    df = "data/benchmark_fmm2d.csv"
    CSV.write(df, DataFrame(n_adapt = [], n_panels = [], n_quad = [], n_points = [], time_mapD = [], time_map = [], time_singleiter = [], time_alliter = []))

    n_adapt = 20

    for n_panel in 4:4:32
        n_quad = n_panel * 2
        mbox = BI.dielectric_mbox2d(eps_boxes, rects, n_panel, n_quad, n_adapt)

        n_points = BI.num_points(mbox)
        @show n_panel, n_quad, n_points

        D = BI.laplace2d_DT_fmm2d(mbox, 1e-4)
        time_mapD = @belapsed $(D) * $(ones(n_points))

        lhs = BI.Lhs_dielectric_mbox2d_fmm2d(mbox, 1e-4)
        rhs = BI.Rhs_dielectric_mbox2d(mbox, r_src, eps_src)

        time_map = @belapsed $(lhs) * $(ones(n_points))

        x, status = Krylov.gmres(lhs, rhs)
        n_iter = status.niter

        time_alliter = @belapsed Krylov.gmres($(lhs), $(rhs), atol = 1e-4)
        time_singleiter = time_alliter / n_iter

        @show time_map, time_singleiter, time_alliter

        CSV.write(df, DataFrame(n_adapt = n_adapt, n_panel = n_panel, n_quad = n_quad, n_points = n_points, time_mapD = time_mapD, time_map = time_map, time_singleiter = time_singleiter, time_alliter = time_alliter), append = true)
    end
end

benchmark()
# benchmark the performance of the BIE solver with KernelAbstractions.jl, including the time for applying the linear operator and the time for the GMRES solver

using CUDA

using BoundaryIntegral
import BoundaryIntegral as BI
using Krylov
using CSV, DataFrames
using LinearAlgebra, Random
using BenchmarkTools

rects = [BI.square(-1.25f0, -0.25f0), BI.square(-0.25f0, 0.25f0), BI.square(0.25f0, -0.75f0), BI.square(-0.75f0, -1.25f0)]
eps_boxes = [2.0f0, 3.0f0, 4.0f0, 5.0f0]
r_src = (0.1f0, 0.2f0)
eps_src = 1.0f0

df = CSV.write("data/benchmark_ka.csv", DataFrame(n_adapt = [], n_panels = [], n_quad = [], n_points = [], time_map = [], time_singleiter = [], time_alliter = []))

n_adapt = 20

for n_panel in 4:4:32
    n_quad = n_panel * 2
    mbox = BI.dielectric_mbox2d(eps_boxes, rects, n_panel, n_quad, n_adapt)

    n_points = BI.num_points(mbox)
    @show n_panel, n_quad, n_points

    lhs = BI.Lhs_dielectric_mbox2d_ka(mbox)
    rhs = BI.Rhs_dielectric_mbox2d_ka(mbox, r_src, eps_src)

    x = CUDA.ones(n_points)
    time_map = @belapsed $(lhs) * $(x)

    x, status = Krylov.gmres(lhs, rhs)
    n_iter = status.niter

    time_alliter = @belapsed Krylov.gmres($(lhs), $(rhs))
    time_singleiter = time_alliter / n_iter

    @show time_map, time_singleiter, time_alliter

    CSV.write("data/benchmark_ka.csv", DataFrame(n_adapt = n_adapt, n_panel = n_panel, n_quad = n_quad, n_points = n_points, time_map = time_map, time_singleiter = time_singleiter, time_alliter = time_alliter), append = true)
end
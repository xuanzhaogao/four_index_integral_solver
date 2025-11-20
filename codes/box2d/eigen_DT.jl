using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
using CairoMakie

rects = [BI.square(-1.0, -0.5), BI.square(0.0, -0.5)]
eps_boxes = [2.0, 3.0]

n_panel = 8
n_quad = 16

n_adapt = 10

mbox = BI.dielectric_mbox2d(eps_boxes, rects, n_panel, n_quad, n_adapt)
DT = BI.laplace2d_DT(mbox)

res = eigen(DT)

begin
    fig = Figure(size = (500, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = "Re", ylabel = "Im")
    scatter!(ax, real.(res.values), imag.(res.values))
    xlims!(ax, (-1.0, 1.0))
    fig
end

save("figs/eigen_DT.svg", fig)
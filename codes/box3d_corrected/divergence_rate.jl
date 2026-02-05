# use a cubic box excited by a far away point source
# see the divergence rate of the edge singularities

include("utils.jl")
using CairoMakie

L = 1.0
l = 0.0049
p = 6
atol = 1e-6
max_order = 128
source = PointSource((10.0, 10.0, 10.0), 10000.0)

box, sigma, _, _, _ = solve_single_box3d_adaptive_mesh_corrected(L, L, L, p, l, 4.0, 1.0, atol, atol, atol, max_order, source, include_edges_src = false, include_edges_trg = false)

sigma_approx = BI.interface_approx(box, sigma)

xs = range(-0.49, 0, length = 7)
ys = range(0.499, 0.4999, length = 100)
z = 0.5

vals = [[sigma_approx((x, y, z)) for y in ys] for x in xs]

begin
    fig = Figure(size = (600, 400), fontsize = 16)
    ax = Axis(fig[1, 1], xlabel = "y", yscale = log10, xscale = log10)

    for (i, x) in enumerate(xs)
        lines!(ax, abs.(0.5 .- ys), abs.(vals[i]), label = "x=$(round(x, digits=2))")
    end

    Legend(fig[1, 2], ax, orientation = :vertical)

    fig
end

save(joinpath(@__DIR__, "figs/divergence_rate.svg"), fig)
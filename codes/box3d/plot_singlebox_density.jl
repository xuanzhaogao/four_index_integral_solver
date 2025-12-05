using BoundaryIntegral
import BoundaryIntegral as BI
using JLD2
using CairoMakie

res = load(joinpath(@__DIR__, "data/single_box3d_ref_non_reduce.jld2"))

dbox = res["dbox"]
sigma = res["sigma"]

points = [p.point for p in BI.eachpoint(dbox)]

density = []
pts_x = []
pts_y = []
# select points with z == 1
for (i, pt) in enumerate(points)
    x, y, z = pt
    if (z ≈ 1.0) && (x >= 0.95) && (y > 0.95)
        push!(density, sigma[i])
        push!(pts_x, x)
        push!(pts_y, y)
    end
end

ys = unique(pts_y)

begin
    y = ys[10]
    mask = pts_y .== y
    xs = pts_x[mask]
    vals = density[mask]
    fig = Figure(size = (500, 400))
    ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y", title = "Density", xscale = log10, yscale = log10)
    scatter!(ax, 1 .- xs, abs.(vals))
    # alpha = 1 - 0.7028479561066333
    alpha = 0.1
    r = vals[end] / ((1 - xs[end]) ^(-alpha))
    lines!(ax, 1 .- xs, (1 .- xs) .^ (-alpha) .* r)
    fig
end

begin
    fig = Figure(size = (500, 400))
    ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y", title = "Density")
    hm = heatmap!(ax, pts_x, pts_y, density)
    Colorbar(fig[1, 2], hm)
    fig
end
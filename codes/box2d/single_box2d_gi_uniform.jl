# viz the error of the Green's identity on uniform points

using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
using CSV, DataFrames
using CairoMakie

rect = BI.square(-0.5, -0.5)

n_quads = [2, 4, 8]
n_adapts = [1, 2, 4]

pts = range(0.0, 0.5, length = 1000)
trgs = zeros(2, length(pts))
for i in eachindex(pts)
    trgs[1, i] = pts[i]
    trgs[2, i] = 0.5
end

gi_trgs = []
for n_quad in n_quads
    gt = []
    for n_adapt in n_adapts
        dbox = BI.dielectric_mbox2d([4.0], [rect], 2, n_quad, n_adapt)
        D_trgs = BI.laplace2d_D_trg_fmm2d(dbox, trgs, 1e-14)
        gi = D_trgs * ones(BI.num_points(dbox))
        push!(gt, gi)
    end
    push!(gi_trgs, gt)
end

begin
    fig = Figure(size = (1000, 800), fontsize = 20)
    ax_1 = Axis(fig[1, 1], xlabel = "x", ylabel = "value", title = "p = 8", xticks = 0.4:0.02:0.5, yticks = -0.55:0.05:-0.2) # fixed n_quad to 8
    ax_2 = Axis(fig[1, 2], xlabel = "x", ylabel = "error", title = "p = 8", yscale = log10, xticks = 0.4:0.02:0.5) 

    ax_3 = Axis(fig[2, 1], xlabel = "x", ylabel = "value", title = "r = 4", xticks = 0.4:0.02:0.5, yticks = -0.55:0.05:-0.2) # fixed n_adapt to 4
    ax_4 = Axis(fig[2, 2], xlabel = "x", ylabel = "error", title = "r = 4", yscale = log10, xticks = 0.4:0.02:0.5) # fixed n_adapt to 4

    for j in 1:3
        lines!(ax_1, pts, gi_trgs[3][j], label = "r = $(n_adapts[j])")
        lines!(ax_2, pts, abs.(0.5 .+ gi_trgs[3][j]))
    end

    for i in 1:3
        lines!(ax_3, pts, gi_trgs[i][3], label = "p = $(n_quads[i])")
        lines!(ax_4, pts, abs.(0.5 .+ gi_trgs[i][3]))
    end
    
    xlims!(ax_1, 0.4, 0.51)
    xlims!(ax_2, 0.4, 0.51)
    xlims!(ax_3, 0.4, 0.51)
    xlims!(ax_4, 0.4, 0.51)

    ylims!(ax_1, -0.54, -0.2)
    ylims!(ax_2, 1e-16, 1.0)
    ylims!(ax_3, -0.55, -0.2)
    ylims!(ax_4, 1e-16, 1.0)

    axislegend(ax_1, position = :lt)
    axislegend(ax_3, position = :lt)

    save("figs/single_box2d_gi_uniform.svg", fig)

    fig
end
using BoundaryIntegral
import BoundaryIntegral as BI
using JLD2
using CairoMakie

res = load(joinpath(@__DIR__, "data/single_box3d_ref.jld"))["res_ref"]
res_2 = load(joinpath(@__DIR__, "data/single_box3d_4_12_6_8.jld"))["res"]
trgs = load(joinpath(@__DIR__, "data/sphere_trgs.jld"))["trgs"]

X = trgs[1, :]
Y = trgs[2, :]
Z = trgs[3, :]

F_ref = res[4]
F_sub = res_2[4]

F_err = log10.(abs.((F_ref - F_sub) ./ F_ref))

begin
    θ = range(-π/2, π/2; length=100)
    φ = range(-π, π; length=100)

    t_g = reshape([t for t in θ, p in φ], 10000)
    p_g = reshape([p for t in θ, p in φ], 10000)
end

begin
    fig = Figure(size = (1200, 400), fontsize = 20)
    ax1 = Axis(fig[1, 1], xlabel = "φ", ylabel = "θ", title = "reference solution", xticks = ([-π, 0, π], ["-π", "0", "π"]), yticks = ([-π/2, 0, π/2], ["-π/2", "0", "π/2"]),aspect = DataAspect())
    ax2 = Axis(fig[1, 3], xlabel = "φ", ylabel = "θ", title = "log10 relative error", xticks = ([-π, 0, π], ["-π", "0", "π"]), yticks = ([-π/2, 0, π/2], ["-π/2", "0", "π/2"]),aspect = DataAspect())

    hm1 = heatmap!(ax1, p_g, t_g, F_ref; colormap=:viridis)
    tightlimits!(ax1)
    Colorbar(fig[1, 2], hm1)

    hm2 = heatmap!(ax2, p_g, t_g, F_err;
    colormap=:viridis)
    tightlimits!(ax2)
    Colorbar(fig[1, 4], hm2)
    fig
end

save(joinpath(@__DIR__, "figs/sphere_heatmap.svg"), fig)
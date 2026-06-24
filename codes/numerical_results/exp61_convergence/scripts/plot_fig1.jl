# 6.1 Fig.-1 system convergence figure: V-error vs DOF + N_iter companion.
# Output: figs/fig61_fig1_convergence.pdf

include(joinpath(@__DIR__, "..", "..", "common", "Lite.jl"))
using .Lite
using CairoMakie, LaTeXStrings, LinearAlgebra
include(joinpath(@__DIR__, "..", "..", "..", "fig_gen", "fig_style.jl"))

const DATA = joinpath(@__DIR__, "..", "data")
const FIGS = joinpath(@__DIR__, "..", "figs")
mkpath(FIGS)

ref = load_ref(joinpath(DATA, "raw", "fig1_ref.jls"))
load_series(p) = [load_ref(joinpath(DATA, "raw", "fig1_p$(p)_r$(r).jls")) for r in 1:5]

const _SC = sweep_colors(3)
const COL = Dict(2 => _SC[1], 4 => _SC[2], 6 => _SC[3])

fig = Figure(size = (FIG_W, FIG_H))
ax1 = Axis(fig[1, 1]; xscale = log10, yscale = log10,
    xlabel = L"\mathrm{DOF}\; N", ylabel = L"|V - V_\mathrm{ref}| / |V_\mathrm{ref}|",
    title = "(a) convergence, slab on two cubes (Fig. 1), ε = 10⁻⁴")
ax2 = Axis(fig[1, 2]; xscale = log10,
    xlabel = L"\mathrm{DOF}\; N", ylabel = L"N_\mathrm{iter}",
    title = "(b) GMRES iterations")

for p in (2, 4, 6)
    s = load_series(p)
    Ns = [d.N for d in s]
    ev = [abs(d.V - ref.V) / abs(ref.V) for d in s]
    scatterlines!(ax1, Ns, ev; color = COL[p], marker = :circle,
        linewidth = LW_DATA, markersize = MS, label = L"p = %$p")
    scatterlines!(ax2, Ns, [d.niter for d in s]; color = COL[p], marker = :circle,
        linewidth = LW_DATA, markersize = MS, label = L"p = %$p")
end
axislegend(ax1; position = :lb, framevisible = false)
ylims!(ax2, 0, 28)
axislegend(ax2; position = :lt, framevisible = false)

save(joinpath(FIGS, "fig61_fig1_convergence.pdf"), fig; px_per_unit = PX_PER_UNIT)
save(joinpath(FIGS, "fig61_fig1_convergence.png"), fig; px_per_unit = 2)
println("wrote figs/fig61_fig1_convergence.pdf")

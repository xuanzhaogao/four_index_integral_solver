# 6.1 Fig.-1 system convergence figure: V-error vs DOF + N_iter companion.
# Output: figs/fig61_fig1_convergence.pdf

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
using CairoMakie, LaTeXStrings, LinearAlgebra

const DATA = joinpath(@__DIR__, "..", "data")
const FIGS = joinpath(@__DIR__, "..", "figs")
mkpath(FIGS)

ref = load_ref(joinpath(DATA, "raw", "fig1_ref.jls"))
load_series(p) = [load_ref(joinpath(DATA, "raw", "fig1_p$(p)_r$(r).jls")) for r in 1:5]

const COL = Dict(2 => "#0072B2", 4 => "#D55E00", 6 => "#009E73")

fig = Figure(size = (820, 330))
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
    scatterlines!(ax1, Ns, ev; color = COL[p], marker = :circle, label = L"p = %$p")
    scatterlines!(ax2, Ns, [d.niter for d in s]; color = COL[p], marker = :circle, label = L"p = %$p")
end
hlines!(ax1, [1e-4]; color = :black, linestyle = :dot)
text!(ax1, 1.0e6, 1.15e-4; text = L"\varepsilon = 10^{-4}", fontsize = 12)
axislegend(ax1; position = :lb, framevisible = false, labelsize = 11)
ylims!(ax2, 0, 28)
axislegend(ax2; position = :lt, framevisible = false, labelsize = 11)

save(joinpath(FIGS, "fig61_fig1_convergence.pdf"), fig)
println("wrote figs/fig61_fig1_convergence.pdf")

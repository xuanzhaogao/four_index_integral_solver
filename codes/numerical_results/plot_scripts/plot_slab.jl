# 6.1 slab convergence figure: V-error vs DOF (three p-curves, with the
# archived no-edge-correction p=6 series as contrast) + N_iter companion.
# Output: figs/fig61_slab_convergence.pdf

include(joinpath(@__DIR__, "..", "common", "Lite.jl"))
using .Lite
using CairoMakie, LaTeXStrings, LinearAlgebra
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))

const DATA = joinpath(@__DIR__, "..", "exp61_convergence", "data")
const FIGS = joinpath(@__DIR__, "..", "figs")
mkpath(FIGS)

ref = load_ref(joinpath(DATA, "raw", "slab_ref.jls"))

load_series(dir, p) = [load_ref(joinpath(dir, "slab_p$(p)_r$(r).jls")) for r in 1:6]
load_series_uc(dir, p) = [load_ref(joinpath(dir, "slab_p$(p)_r$(r).jls")) for r in 1:5]

# Okabe-Ito
const _SC = sweep_colors(3)
const COL = Dict(2 => _SC[1], 4 => _SC[2], 6 => _SC[3])

fig = Figure(size = (FIG_W, FIG_H))
ax1 = Axis(fig[1, 1]; xscale = log10, yscale = log10,
    xlabel = L"\mathrm{DOF}\; N", ylabel = L"|V - V_\mathrm{ref}| / |V_\mathrm{ref}|",
    title = "(a) convergence, slab 10×10×1, ε = 10⁻⁴")
ax2 = Axis(fig[1, 2]; xscale = log10,
    xlabel = L"\mathrm{DOF}\; N", ylabel = L"N_\mathrm{iter}",
    title = "(b) GMRES iterations")

for p in (2, 4, 6)
    s = load_series(joinpath(DATA, "raw"), p)
    Ns = [d.N for d in s]
    ev = [abs(d.V - ref.V) / abs(ref.V) for d in s]
    its = [d.niter for d in s]
    scatterlines!(ax1, Ns, ev; color = COL[p], marker = :circle,
        linewidth = LW_DATA, markersize = MS, label = L"p = %$p")
    scatterlines!(ax2, Ns, its; color = COL[p], marker = :circle,
        linewidth = LW_DATA, markersize = MS, label = L"p = %$p")
end

# archived no-edge-correction series (errors recomputed vs the converged ref)
for p in (2, 4, 6)
    s = load_series_uc(joinpath(DATA, "raw", "noedges_run1"), p)
    Ns = [d.N for d in s]
    ev = [abs(d.V - ref.V) / abs(ref.V) for d in s]
    its = [d.niter for d in s]
    scatterlines!(ax1, Ns, ev; color = (COL[p], 0.55), marker = :utriangle,
        markersize = MS, linewidth = LW_GUIDE, linestyle = :dash)
    scatterlines!(ax2, Ns, its; color = (COL[p], 0.55), marker = :utriangle,
        markersize = MS, linewidth = LW_GUIDE, linestyle = :dash)
end
# legend proxy for the dashed family
scatterlines!(ax1, [NaN], [NaN]; color = :gray40, marker = :utriangle,
    linewidth = LW_GUIDE, linestyle = :dash, label = "no edge corr.")

axislegend(ax1; position = :lb, framevisible = false)
ylims!(ax2, 0, 28)
axislegend(ax2; position = :lt, framevisible = false)

save(joinpath(FIGS, "fig61_slab_convergence.pdf"), fig; px_per_unit = PX_PER_UNIT)
save(joinpath(FIGS, "fig61_slab_convergence.png"), fig; px_per_unit = 2)
println("wrote figs/fig61_slab_convergence.pdf")

# console: no-edges errors vs converged ref, for the record
println("\nno-edge-correction series vs converged reference:")
for p in (2, 4, 6)
    s = load_series_uc(joinpath(DATA, "raw", "noedges_run1"), p)
    println("  p=$p: ", join([string(round(abs(d.V - ref.V) / abs(ref.V); sigdigits = 3)) for d in s], "  "))
end

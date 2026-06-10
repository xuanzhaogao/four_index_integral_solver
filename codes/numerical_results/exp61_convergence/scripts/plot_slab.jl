# 6.1 slab convergence figure: V-error vs DOF (three p-curves, with the
# archived no-edge-correction p=6 series as contrast) + N_iter companion.
# Output: figs/fig61_slab_convergence.pdf

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
using CairoMakie, LaTeXStrings, LinearAlgebra

const DATA = joinpath(@__DIR__, "..", "data")
const FIGS = joinpath(@__DIR__, "..", "figs")
mkpath(FIGS)

ref = load_ref(joinpath(DATA, "raw", "slab_ref.jls"))

load_series(dir, p) = [load_ref(joinpath(dir, "slab_p$(p)_r$(r).jls")) for r in 1:5]

# Okabe-Ito
const COL = Dict(2 => "#0072B2", 4 => "#D55E00", 6 => "#009E73")

fig = Figure(size = (820, 330))
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
    scatterlines!(ax1, Ns, ev; color = COL[p], marker = :circle, label = L"p = %$p")
    scatterlines!(ax2, Ns, its; color = COL[p], marker = :circle, label = L"p = %$p")
end

# archived no-edge-correction series (errors recomputed vs the converged ref)
for p in (6,)
    s = load_series(joinpath(DATA, "raw", "noedges_run1"), p)
    Ns = [d.N for d in s]
    ev = [abs(d.V - ref.V) / abs(ref.V) for d in s]
    its = [d.niter for d in s]
    scatterlines!(ax1, Ns, ev; color = (:gray40, 0.9), marker = :utriangle,
        linestyle = :dash, label = L"p = 6\;\mathrm{(no\;edge\;corr.)}")
    scatterlines!(ax2, Ns, its; color = (:gray40, 0.9), marker = :utriangle, linestyle = :dash)
end

hlines!(ax1, [1e-4]; color = :black, linestyle = :dot)
text!(ax1, 3.0e5, 1.22e-4; text = L"\varepsilon = 10^{-4}", fontsize = 12)
axislegend(ax1; position = :lb, framevisible = false, labelsize = 11)
ylims!(ax2, 0, 28)
axislegend(ax2; position = :lt, framevisible = false, labelsize = 11)

save(joinpath(FIGS, "fig61_slab_convergence.pdf"), fig)
println("wrote figs/fig61_slab_convergence.pdf")

# console: no-edges errors vs converged ref, for the record
println("\nno-edge-correction series vs converged reference:")
for p in (2, 4, 6)
    s = load_series(joinpath(DATA, "raw", "noedges_run1"), p)
    println("  p=$p: ", join([string(round(abs(d.V - ref.V) / abs(ref.V); sigdigits = 3)) for d in s], "  "))
end

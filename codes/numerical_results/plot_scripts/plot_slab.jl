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
const _SM = sweep_markers(3)
const COL = Dict(2 => _SC[1], 4 => _SC[2], 6 => _SC[3])
const MK  = Dict(2 => _SM[1], 4 => _SM[2], 6 => _SM[3])

begin

    fig = Figure(size = (FIG_W, FIG_H))
    ax1 = Axis(fig[1, 1]; xscale = log10, yscale = log10,
        xlabel = L"N", ylabel = L"\mathcal{E}_{r}")
    ax2 = Axis(fig[1, 2]; xscale = log10,
        xlabel = L"N", ylabel = L"N_\mathrm{iter}")
    for (ax, lab) in ((ax1, "(a)"), (ax2, "(b)"))
        text!(ax, 0, 1; text = lab, space = :relative, align = (:left, :top),
            offset = (6, -6), font = :bold, fontsize = FS_BASE)
    end

    for p in (2, 4, 6)
        s = load_series(joinpath(DATA, "raw"), p)
        Ns = [d.N for d in s]
        ev = [abs(d.V - ref.V) / abs(ref.V) for d in s]
        its = [d.niter for d in s]
        scatterlines!(ax1, Ns, ev; color = COL[p], marker = MK[p],
            linewidth = LW_DATA, markersize = MS, label = L"p = %$p")
        scatterlines!(ax2, Ns, its; color = COL[p], marker = MK[p],
            linewidth = LW_DATA, markersize = MS, label = L"p = %$p")
    end

    # archived no-edge-correction series (errors recomputed vs the converged ref):
    # same color + marker per p, distinguished by the dashed (faded) line
    for p in (2, 4, 6)
        s = load_series_uc(joinpath(DATA, "raw", "noedges_run1"), p)
        Ns = [d.N for d in s]
        ev = [abs(d.V - ref.V) / abs(ref.V) for d in s]
        its = [d.niter for d in s]
        scatterlines!(ax1, Ns, ev; color = (COL[p], 0.55), marker = MK[p],
            markersize = MS, linewidth = LW_GUIDE, linestyle = :dash,
            label = L"p = %$p \text{ (no edge corr.)}")
        scatterlines!(ax2, Ns, its; color = (COL[p], 0.55), marker = MK[p],
            markersize = MS, linewidth = LW_GUIDE, linestyle = :dash)
    end

    axislegend(ax1; position = :lb)
    ylims!(ax2, 0, 30)
    # axislegend(ax2; position = :rb)

    save(joinpath(FIGS, "fig61_slab_convergence.pdf"), fig; px_per_unit = PX_PER_UNIT)
    
    fig
end

# console: no-edges errors vs converged ref, for the record
println("\nno-edge-correction series vs converged reference:")
for p in (2, 4, 6)
    s = load_series_uc(joinpath(DATA, "raw", "noedges_run1"), p)
    println("  p=$p: ", join([string(round(abs(d.V - ref.V) / abs(ref.V); sigdigits = 3)) for d in s], "  "))
end

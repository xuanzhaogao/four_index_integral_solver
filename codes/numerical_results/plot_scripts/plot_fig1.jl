# 6.1 Fig.-1 system convergence figure: V-error vs DOF + N_iter companion.
# Output: figs/fig61_fig1_convergence.pdf

include(joinpath(@__DIR__, "..", "common", "Lite.jl"))
using .Lite
using CairoMakie, LaTeXStrings, LinearAlgebra
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))

const DATA = joinpath(@__DIR__, "..", "exp61_convergence", "data")
const FIGS = joinpath(@__DIR__, "..", "figs")
mkpath(FIGS)

ref = load_ref(joinpath(DATA, "raw", "fig1_ref.jls"))
load_series(p) = [load_ref(joinpath(DATA, "raw", "fig1_p$(p)_r$(r).jls")) for r in 1:5]
# no-edge-correction contrast archive (run_fig1_noedges.jl); load whatever exists
const NOEDGE = joinpath(DATA, "raw", "fig1_noedges")
load_noedge(p) = [load_ref(joinpath(NOEDGE, "fig1_p$(p)_r$(r).jls"))
                  for r in 1:5 if isfile(joinpath(NOEDGE, "fig1_p$(p)_r$(r).jls"))]

const _SC = sweep_colors(3)
const _SM = sweep_markers(3)
const COL = Dict(2 => _SC[1], 4 => _SC[2], 6 => _SC[3])
const MK  = Dict(2 => _SM[1], 4 => _SM[2], 6 => _SM[3])

begin
    fig = Figure(size = (FIG_W, FIG_H), fontsize = FS_BASE)
    ax1 = Axis(fig[1, 1]; xscale = log10, yscale = log10,
        xlabel = L"N", ylabel = L"\mathcal{E}_{r}")
    ax2 = Axis(fig[1, 2]; xscale = log10,
        xlabel = L"N", ylabel = L"N_\mathrm{iter}")
    for (ax, lab) in ((ax1, "(a)"), (ax2, "(b)"))
        text!(ax, 0, 1; text = lab, space = :relative, align = (:left, :top),
            offset = (6, -6), font = :bold, fontsize = FS_BASE)
    end

    for p in (2, 4, 6)
        s = load_series(p)
        Ns = [d.N for d in s]
        ev = [abs(d.V - ref.V) / abs(ref.V) for d in s]
        scatterlines!(ax1, Ns, ev; color = COL[p], marker = MK[p],
            linewidth = LW_DATA, markersize = MS, label = L"p = %$p")
        scatterlines!(ax2, Ns, [d.niter for d in s]; color = COL[p], marker = MK[p],
            linewidth = LW_DATA, markersize = MS)
    end

    # no-edge-correction contrast (correct_edges=false): same color + marker per
    # p as the corrected curve, distinguished by the dashed (faded) line
    for p in (2, 4)
        s = load_noedge(p)
        isempty(s) && continue
        Ns = [d.N for d in s]
        ev = [abs(d.V - ref.V) / abs(ref.V) for d in s]
        scatterlines!(ax1, Ns, ev; color = (COL[p], 0.55), marker = MK[p],
            markersize = MS, linewidth = LW_GUIDE, linestyle = :dash,
            label = L"p = %$p \text{ (no edge corr.)}")

        scatterlines!(ax2, Ns, [d.niter for d in s]; color = (COL[p], 0.55),
            marker = MK[p], markersize = MS, linewidth = LW_GUIDE, linestyle = :dash)
    end

    axislegend(ax1; position = :lb)
    ylims!(ax1, 10^(-4.5), 10^(-1.5))
    ylims!(ax2, 15, 30)

    save(joinpath(FIGS, "fig61_fig1_convergence.pdf"), fig; px_per_unit = PX_PER_UNIT)

    fig
end

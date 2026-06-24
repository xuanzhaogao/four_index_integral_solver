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

const _SC = sweep_colors(3)
const COL = Dict(2 => _SC[1], 4 => _SC[2], 6 => _SC[3])

begin
    fig = Figure(size = (FIG_W, FIG_H), fontsize = FS_BASE)
    ax1 = Axis(fig[1, 1]; xscale = log10, yscale = log10,
        xlabel = "DOF", ylabel = L"\mathcal{E}_{r}")
    ax2 = Axis(fig[1, 2]; xscale = log10,
        xlabel = "DOF", ylabel = L"N_\mathrm{iter}")
    for (ax, lab) in ((ax1, "(a)"), (ax2, "(b)"))
        text!(ax, 0, 1; text = lab, space = :relative, align = (:left, :top),
            offset = (6, -6), font = :bold, fontsize = FS_BASE)
    end

    for p in (2, 4, 6)
        s = load_series(p)
        Ns = [d.N for d in s]
        ev = [abs(d.V - ref.V) / abs(ref.V) for d in s]
        scatterlines!(ax1, Ns, ev; color = COL[p], marker = :circle,
            linewidth = LW_DATA, markersize = MS, label = "p = $p")
        scatterlines!(ax2, Ns, [d.niter for d in s]; color = COL[p], marker = :circle,
            linewidth = LW_DATA, markersize = MS)
    end
    axislegend(ax1; position = :lb)
    ylims!(ax1, 10^(-4.5), 10^(-2))
    ylims!(ax2, 0, 30)
    # axislegend(ax2; position = :rb)

    save(joinpath(FIGS, "fig61_fig1_convergence.pdf"), fig; px_per_unit = PX_PER_UNIT)

    fig
end

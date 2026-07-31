# 6.3 figure: (a) V-error vs DOF, one curve per eps2 (p = 6, r = 1..5);
# (b) GMRES iterations vs DOF per eps2.
# Output: figs/fig63_contrast.pdf

include(joinpath(@__DIR__, "..", "common", "Lite.jl"))
using .Lite
using CairoMakie, LaTeXStrings
include(joinpath(@__DIR__, "..", "..", "fig_gen", "fig_style.jl"))

# RERUN_TAG: plot the rerun tree (data_v2/, *_v2 campaigns) instead of the published one.
const TAG = get(ENV, "RERUN_TAG", "")
const DATA = joinpath(@__DIR__, "..", "exp63_contrast", "data" * TAG)
const RAW = joinpath(DATA, "raw")
const FIGS = joinpath(@__DIR__, "..", "figs" * TAG); mkpath(FIGS)
mkpath(FIGS)

e2s = sort([parse(Float64, match(r"^ratio_ref_e2_(.+)\.jls$", f).captures[1])
            for f in readdir(RAW) if occursin(r"^ratio_ref_e2_.*\.jls$", f)])
e2s = filter(!=(2.0), e2s)        # eps2 = 2 excluded from the figure (kept in data/CSV)
rs_of(e2) = sort([parse(Int, match(r"_r(\d+)\.jls$", f).captures[1])
                  for f in readdir(RAW) if startswith(f, "ratio_test_e2_$(e2)_r")])

# Tol-bright sweep palette, one color + marker per eps2 (marker keeps the
# curves separable in grayscale / without color)
col(i) = sweep_colors(5)[mod1(i, 5)]
mk(i)  = sweep_markers(5)[mod1(i, 5)]

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

    for (i, e2) in enumerate(e2s)
        rf = load_ref(joinpath(RAW, "ratio_ref_e2_$(e2).jls"))
        ts = [load_ref(joinpath(RAW, "ratio_test_e2_$(e2)_r$(r).jls")) for r in rs_of(e2)]
        Ns = [t.N for t in ts]
        ev = [abs(t.V - rf.V) / abs(rf.V) for t in ts]
        g = (e2 - 4.0) / (e2 + 4.0)
        lab = L"\varepsilon_2 = %$(round(Int, e2))"
        scatterlines!(ax1, Ns, ev; color = col(i), marker = mk(i),
            linewidth = LW_DATA, markersize = MS, label = lab)
        scatterlines!(ax2, Ns, [t.niter for t in ts]; color = col(i), marker = mk(i),
            linewidth = LW_DATA, markersize = MS)
    end
    axislegend(ax1; position = :lb)

    ylims!(ax1, 10^(-4.5), 10^(-1.5))
    ylims!(ax2, 0, 50)
    fig
end

save(joinpath(FIGS, "fig63_contrast.pdf"), fig; px_per_unit = PX_PER_UNIT)
println("wrote figs/fig63_contrast.pdf")

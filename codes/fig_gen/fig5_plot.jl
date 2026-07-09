#=
Figure 5 plot script. Loads fig5_data.jls and writes fig5_tkm_validation.png.

Layout (1 x 2):
  (a) relative L2 error vs grid resolution n, comparing TKM at several
      tolerances eps and FMM (representative)
  (b) relative L2 error vs eta = min_alpha P_alpha / (l_alpha + L) at fixed
      n = 64, for the same tolerances eps
=#

using Serialization
using CairoMakie
using LaTeXStrings
using Printf

include(joinpath(@__DIR__, "fig_style.jl"))

const datapath = joinpath(@__DIR__, "fig5_data.jls")
const data = open(deserialize, datapath, "r")

const n_list      = data.n_list
const eta_list    = data.eta_list
const eps_list    = data.eps_list
const err_tkm_n   = data.err_tkm_n
const err_fmm_n   = data.err_fmm_n
const err_tkm_eta = data.err_tkm_eta
const n_for_eta   = data.n_for_eta

const eps_palette  = sweep_colors(length(eps_list))
const eps_colors   = Dict(eps_list[i] => eps_palette[i] for i in 1:length(eps_list))
const eps_mk_vec   = sweep_markers(length(eps_list))
const eps_markers  = Dict(eps_list[i] => eps_mk_vec[i] for i in 1:length(eps_list))

eps_label(eps) = L"\varepsilon = 10^{%$(round(Int, log10(eps)))}"

begin
    fig = Figure(size = (FIG_W, FIG_H))

    floor_y = 1e-16
    clip(y) = max(y, floor_y)

    # ---------------------------------------------------------------------
    # Panel (a): convergence vs n
    # ---------------------------------------------------------------------
    ax_a = Axis(fig[1, 1];
                yscale = log10,
                yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5),
                xlabel = L"n",
                ylabel = L"\mathcal{E}_u",
                xticks = (n_list, string.(n_list)),
                yticks = ([1e-13, 1e-11, 1e-9, 1e-7, 1e-5, 1e-3, 1e-1, 1e1],
                          [rich("10", superscript(string(e))) for e in (-13, -11, -9, -7, -5, -3, -1, 1)]))

    for eps in eps_list
        scatterlines!(ax_a, n_list, clip.(err_tkm_n[eps]);
                      color = eps_colors[eps], marker = eps_markers[eps],
                      markersize = MS, linewidth = LW_DATA,
                      label = L"TKM, $%$(eps_label(eps).s)$")
    end
    scatterlines!(ax_a, n_list, clip.(err_fmm_n[eps_list[end]]);
                  color = :black, marker = :star5,
                  markersize = MS, linewidth = LW_DATA,
                  linestyle = :dash,
                  label = L"\text{Direct Sum}")

    axislegend(ax_a; position = :lb)
    ylims!(ax_a, 1e-13, 1e1)

    # ---------------------------------------------------------------------
    # Panel (b): eta sweep at fixed n
    # ---------------------------------------------------------------------
    ax_b = Axis(fig[1, 2];
                # xscale = log10, 
                yscale = log10,
                # x has hand-picked, non-uniform major ticks (0.5,0.7,1.0,...) so
                # IntervalsBetween minor ticks land irregularly; keep only the y-minor grid.
                xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5),
                yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5),
                xlabel = L"\eta = \min_{\alpha}\, P_\alpha / (l_\alpha + L)",
                ylabel = L"\mathcal{E}_u",
                # xticks = ([0.5, 0.7, 1.0, 1.5, 2.0, 2.5],
                        #   ["0.5","0.7","1.0","1.5","2.0","2.5"]),
                yticks = ([1e-13, 1e-11, 1e-9, 1e-7, 1e-5, 1e-3, 1e-1, 1e1],
                          [rich("10", superscript(string(e))) for e in (-13, -11, -9, -7, -5, -3, -1, 1)]))

    for eps in eps_list
        scatterlines!(ax_b, eta_list[4:17], clip.(err_tkm_eta[eps])[4:17];
                      color = eps_colors[eps], marker = eps_markers[eps],
                      markersize = MS, linewidth = LW_DATA,
                      label = eps_label(eps))
    end

    vlines!(ax_b, [1.0]; color = (:black, 0.6),
            linestyle = :dash, linewidth = LW_GUIDE)
    text!(ax_b, L"\eta = 1"; position = (1.04, 1e-1),
          fontsize = FS_ANNOT, color = :black)

    axislegend(ax_b; position = :lb)
    xlims!(ax_b, 0.55, 1.5)
    ylims!(ax_b, 1e-13, 1e1)

    colgap!(fig.layout, 1, 30)

    for (ax, lab) in ((ax_a, "(a)"), (ax_b, "(b)"))
        text!(ax, 0, 1; text = lab, space = :relative, align = (:left, :top),
              offset = (6, -6), font = :bold, fontsize = FS_BASE)
    end

    outpath = joinpath(@__DIR__, "figs/fig5_tkm_validation.pdf")
    save(outpath, fig; px_per_unit = PX_PER_UNIT)
    @info "Saved figure" outpath

    fig
end

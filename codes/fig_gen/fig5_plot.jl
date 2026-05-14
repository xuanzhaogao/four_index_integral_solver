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

const datapath = joinpath(@__DIR__, "fig5_data.jls")
const data = open(deserialize, datapath, "r")

const n_list      = data.n_list
const eta_list    = data.eta_list
const eps_list    = data.eps_list
const err_tkm_n   = data.err_tkm_n
const err_fmm_n   = data.err_fmm_n
const err_tkm_eta = data.err_tkm_eta
const n_for_eta   = data.n_for_eta

const eps_colors = Dict(
    eps_list[1] => :crimson,
    eps_list[2] => :royalblue,
    eps_list[3] => :seagreen,
    eps_list[4] => :purple,
)

eps_label(eps) = L"\varepsilon = 10^{%$(round(Int, log10(eps)))}"

begin
    fig = Figure(size = (1000, 450), fontsize = 18)

    floor_y = 1e-16
    clip(y) = max(y, floor_y)

    # ---------------------------------------------------------------------
    # Panel (a): convergence vs n
    # ---------------------------------------------------------------------
    ax_a = Axis(fig[1, 1];
                yscale = log10,
                xlabel = L"n",
                ylabel = L"\mathcal{E}_u",
                xticks = (n_list, string.(n_list)))

    for eps in eps_list
        scatterlines!(ax_a, n_list, clip.(err_tkm_n[eps]);
                      color = eps_colors[eps], marker = :circle,
                      markersize = 12, linewidth = 2,
                      label = L"TKM, $%$(eps_label(eps).s)$")
    end
    scatterlines!(ax_a, n_list, clip.(err_fmm_n[eps_list[end]]);
                  color = :black, marker = :rect,
                  markersize = 12, linewidth = 2,
                  linestyle = :dash,
                  label = L"\text{FMM}")

    axislegend(ax_a; position = :lb)

    # ---------------------------------------------------------------------
    # Panel (b): eta sweep at fixed n
    # ---------------------------------------------------------------------
    ax_b = Axis(fig[1, 2];
                xscale = log10, yscale = log10,
                xlabel = L"\eta = \min_{\alpha}\, P_\alpha / (l_\alpha + L)",
                ylabel = L"\mathcal{E}_u",
                xticks = ([0.5, 0.7, 1.0, 1.5, 2.0, 2.5],
                          ["0.5","0.7","1.0","1.5","2.0","2.5"]))

    for eps in eps_list
        scatterlines!(ax_b, eta_list, clip.(err_tkm_eta[eps]);
                      color = eps_colors[eps], marker = :circle,
                      markersize = 12, linewidth = 2,
                      label = eps_label(eps))
    end

    vlines!(ax_b, [1.0]; color = (:black, 0.6),
            linestyle = :dash, linewidth = 1.5)
    text!(ax_b, L"\eta = 1"; position = (1.04, 1e-1),
          fontsize = 18, color = :black)

    axislegend(ax_b; position = :lb)

    colgap!(fig.layout, 1, 30)

    outpath = joinpath(@__DIR__, "fig5_tkm_validation.png")
    save(outpath, fig; px_per_unit = 4)
    @info "Saved figure" outpath

    fig
end

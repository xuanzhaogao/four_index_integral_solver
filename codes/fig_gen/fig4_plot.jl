#=
Figure 4 plot script. Loads fig4_data.jls and writes fig4_near_correction.png.

Layout:
  (a) E_near(d/h):  standard p x p GL  vs  near-field-corrected (c = 5, 8),
                    with the dynamically-upsampled curve as a floor reference
  (b) max residual error left in d/h > c by the corrected scheme, vs c
=#

using Serialization
using CairoMakie
using LaTeXStrings
using Printf

const datapath = joinpath(@__DIR__, "fig4_data.jls")
const data = open(deserialize, datapath, "r")

const p_quad   = data.p_quad
const eps_corr = data.eps_corr
const dh       = data.dh_list
const Estd     = data.E_std
const Eup      = data.E_up
const Ecorr    = data.E_corr
const c_scan   = data.c_scan
const Eresid_c = data.max_resid_far_c

begin
    fig = Figure(size = (1000, 450), fontsize = 18)

    # ---------------------------------------------------------------------
    # Panel (a): E_near(d/h)
    # ---------------------------------------------------------------------
    ax_a = Axis(fig[1, 1];
                xscale = log10, yscale = log10,
                xlabel = L"d / h",
                ylabel = L"E_{\mathrm{near}}")
                # xticks = ([0.2, 0.5, 1, 2, 5, 10, 20],
                          # ["0.2","0.5","1","2","5","10","20"]))

    floor_y = 1e-16
    clip(y) = max(y, floor_y)

    scatterlines!(ax_a, dh, clip.(Estd);
                  color = :crimson, marker = :rect,
                  markersize = 12, linewidth = 2,
                  label = L"standard $p \times p$ GL")

    scatterlines!(ax_a, dh, clip.(Ecorr[5]);
                  color = :royalblue, marker = :circle,
                  markersize = 12, linewidth = 2,
                  label = L"corrected, $c = 5$")

    scatterlines!(ax_a, dh, clip.(Ecorr[8]);
                  color = :seagreen, marker = :diamond,
                  markersize = 12, linewidth = 2,
                  label = L"corrected, $c = 8$")

    # eps_str = @sprintf("%.0e", eps_corr)
    # lines!(ax_a, dh, clip.(Eup);
    #        color = (:gray, 0.6), linewidth = 1.5, linestyle = :dash,
    #        label = L"dynamic $p_{\mathrm{up}}$ ($\varepsilon = %$(eps_str)$)")

    axislegend(ax_a; position = :rt) 

    # ---------------------------------------------------------------------
    # Panel (b): residual error in d/h > c (uncorrected far region)
    # ---------------------------------------------------------------------
    ax_b = Axis(fig[1, 2];
                xlabel = L"near-field threshold $c$",
                ylabel = L"\max_{d/h > c}\; E_{\mathrm{std}}",
                yscale = log10,
                xticks = (c_scan, string.(c_scan)))

    scatterlines!(ax_b, c_scan, max.(Eresid_c, 1e-16);
                  color = :crimson, marker = :circle,
                  markersize = 12, linewidth = 2)

    # for (cc, col, lab) in ((5, :royalblue, L"c = 5"),
    #                         (8, :seagreen,  L"c = 8"))
    #     vlines!(ax_b, [cc]; color = (col, 0.7),
    #             linestyle = :dash, linewidth = 1.5)
    #     text!(ax_b, lab; position = (cc + 0.15, 5e-2),
    #           fontsize = 16, color = col)
    # end

    colgap!(fig.layout, 1, 30)

    outpath = joinpath(@__DIR__, "fig4_near_correction.png")
    save(outpath, fig; px_per_unit = 4)
    @info "Saved figure" outpath

    fig
end

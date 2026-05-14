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
# Direct standard-GL error at d = c·h for each c. This is the residual the
# uncorrected far region contributes right at the near-field threshold,
# evaluated without log-grid binning bias.
const Eresid_c = data.E_std_at_c

begin
    fig = Figure(size = (1000, 450), fontsize = 18)

    # ---------------------------------------------------------------------
    # Panel (a): E_near(d/h)
    # ---------------------------------------------------------------------
    ax_a = Axis(fig[1, 1];
                xscale = log10, yscale = log10,
                xlabel = L"d / h",
                ylabel = L"\mathcal{E}_{\mathrm{near}}")
                # xticks = ([0.2, 0.5, 1, 2, 5, 10, 20],
                          # ["0.2","0.5","1","2","5","10","20"]))

    floor_y = 1e-16
    clip(y) = max(y, floor_y)

    scatterlines!(ax_a, dh, clip.(Estd);
                  color = :crimson, marker = :rect,
                  markersize = 12, linewidth = 2,
                  label = L"standard $8 \times 8$ GL")

    # Iterate over the c values present in the saved E_corr dict. Sorted so
    # color/marker assignment is stable across runs.
    corr_styles = [(:royalblue, :circle), (:seagreen, :diamond),
                   (:darkorange, :utriangle), (:purple, :star5)]
    for (i, cc) in enumerate(sort(collect(keys(Ecorr))))
        col, mk = corr_styles[mod1(i, length(corr_styles))]
        scatterlines!(ax_a, dh, clip.(Ecorr[cc]);
                      color = col, marker = mk,
                      markersize = 12, linewidth = 2,
                      label = L"corrected, $c = %$cc$")
    end

    # eps_str = @sprintf("%.0e", eps_corr)
    # lines!(ax_a, dh, clip.(Eup);
    #        color = (:gray, 0.6), linewidth = 1.5, linestyle = :dash,
    #        label = L"dynamic $p_{\mathrm{up}}$ ($\varepsilon = %$(eps_str)$)")

    axislegend(ax_a; position = :rt) 

    # ---------------------------------------------------------------------
    # Panel (b): residual error in d/h > c (uncorrected far region)
    # ---------------------------------------------------------------------
    ax_b = Axis(fig[1, 2];
                xlabel = L"$c$",
                ylabel = L"\mathcal{E}_{\mathrm{std}}",
                yscale = log10,
                xticks = (c_scan, string.(c_scan)))

    scatterlines!(ax_b, c_scan, max.(Eresid_c, 1e-16);
                  color = :crimson, marker = :circle,
                  markersize = 12, linewidth = 2,
                  label = L"E_{\mathrm{std}}")

    # Linear fit on the semi-log plane: log10(E) ≈ a + b·c. The slope b is
    # decades-per-unit-c, i.e. 10^b is the per-step decay ratio.
    # Restrict to the decay regime (drop points within one decade of the
    # observed floor) so the flattening at large c doesn't blunt the slope.
    fit_mask = isfinite.(Eresid_c) .& (Eresid_c .> 0)
    floor_obs = minimum(Eresid_c[fit_mask])
    fit_mask .&= Eresid_c .> 10 * floor_obs
    cs_fit   = Float64.(c_scan[fit_mask])
    ys_fit   = log10.(Eresid_c[fit_mask])
    n_fit    = length(cs_fit)
    c_mean   = sum(cs_fit) / n_fit
    y_mean   = sum(ys_fit) / n_fit
    b_slope  = sum((cs_fit .- c_mean) .* (ys_fit .- y_mean)) /
               sum((cs_fit .- c_mean).^2)
    a_int    = y_mean - b_slope * c_mean
    c_line   = range(minimum(cs_fit), maximum(cs_fit); length = 64)
    y_line   = 10 .^ (a_int .+ b_slope .* c_line)
    slope_str  = @sprintf("%.2f", b_slope)
    fit_label  = LaTeXString("fit: slope = " * slope_str * " per unit \$c\$")
    # lines!(ax_b, c_line, y_line;
    #        color = :black, linestyle = :dash, linewidth = 1.5,
    #        label = fit_label)
    # axislegend(ax_b; position = :rt, framevisible = false)
    # @info "Panel (b) linear fit" slope = b_slope intercept = a_int

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

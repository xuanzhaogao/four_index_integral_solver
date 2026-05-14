#=
Figure 6 plot script. Loads pre-computed data from fig6_data.jls.

Layout (1 x 2):
  (a) Relative error E_V vs post-refinement threshold h₀, one line per
      near-cutoff multiplier c.
  (b) N_FMM (solid, circle) and N_HCub (dashed, square) vs h₀, one colour
      per c.
=#

using Serialization
using CairoMakie
using LaTeXStrings
using Printf

const datapath = joinpath(@__DIR__, "fig6_data.jls")
const data = open(deserialize, datapath, "r")

const cs   = data.c_list
const h0s_all = data.h0_list
const c_colors = Dict(
    cs[1] => :crimson,
    cs[2] => :royalblue,
    cs[3] => :seagreen,
    cs[4] => :purple,
    cs[5] => :darkorange,
    cs[6] => :goldenrod,
)
const h0_colors = Dict(
    h0s_all[1] => :crimson,
    h0s_all[2] => :royalblue,
    h0s_all[3] => :seagreen,
    h0s_all[4] => :purple,
    h0s_all[5] => :darkorange,
)

c_label(c)   = L"c = %$(Int(c))"
h0_label(h0) = L"h_0 = 2^{%$(round(Int, log2(h0)))}"

begin
    fig = Figure(size = (1000, 450), fontsize = 18)

    floor_y = 1e-16
    clip(y) = max(y, floor_y)

    # ---------------------------------------------------------------------
    # Panel (a): E_V vs c, one line per h₀
    # ---------------------------------------------------------------------
    ax_a = Axis(fig[1, 1];
                yscale = log10,
                xlabel = L"c",
                ylabel = L"\mathcal{E}_V",
                xticks = (cs, [string(Int(c)) for c in cs]))

    for h0 in h0s_all
        Es = Float64[]
        for c in cs
            row = first(r for r in data.results[c] if r.h0 == h0)
            push!(Es, clip(row.E_V))
        end
        scatterlines!(ax_a, cs, Es;
                      color = h0_colors[h0], marker = :circle,
                      markersize = 12, linewidth = 2,
                      label = h0_label(h0))
    end
    axislegend(ax_a; position = :rt)

    # ---------------------------------------------------------------------
    # Panel (b): N_FMM (solid) and N_HCub (dashed) vs h₀
    # ---------------------------------------------------------------------
    ax_b = Axis(fig[1, 2];
                xscale = log10, yscale = log10,
                xlabel = L"h_0", ylabel = "count")
    ax_b.xreversed = true

    for c in cs
        rows = data.results[c]
        h0s = [r.h0 for r in rows]
        NFMM  = [r.N_FMM  for r in rows]
        NHcub = [r.N_HCub for r in rows]
        scatterlines!(ax_b, h0s, NFMM;
                      color = c_colors[c], marker = :circle,
                      markersize = 12, linewidth = 2)
        scatterlines!(ax_b, h0s, NHcub;
                      color = c_colors[c], marker = :rect,
                      markersize = 12, linewidth = 2, linestyle = :dash)
    end

    # Inline legend: marker shapes for FMM vs HCub, plus c colour swatches.
    style_elems = [
        [MarkerElement(color = :black, marker = :circle, markersize = 12),
         LineElement(color = :black, linewidth = 2, linestyle = :solid)],
        [MarkerElement(color = :black, marker = :rect, markersize = 12),
         LineElement(color = :black, linewidth = 2, linestyle = :dash)],
    ]
    color_elems = [LineElement(color = c_colors[c], linewidth = 2) for c in cs]
    axislegend(ax_b,
        [style_elems..., color_elems...],
        [L"N_{\mathrm{FMM}}", L"N_{\mathrm{HCub}}", [c_label(c) for c in cs]...];
        position = :rt, nbanks = 2, labelsize = 14)

    colgap!(fig.layout, 1, 30)

    outpath = joinpath(@__DIR__, "fig6_post_refinement.png")
    save(outpath, fig; px_per_unit = 4)
    @info "Saved figure" outpath

    fig
end

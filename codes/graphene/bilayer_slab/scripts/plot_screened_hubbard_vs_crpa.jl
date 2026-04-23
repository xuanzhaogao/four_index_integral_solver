using CairoMakie

include(joinpath(@__DIR__, "..", "src", "ScreenedHubbardComparison.jl"))
using .ScreenedHubbardComparison

const INPUT_FILE = DEFAULT_SCREENED_HUBBARD_CSV
const OUTPUT_FILE = joinpath(@__DIR__, "..", "figs", "screened_hubbard_vs_crpa.png")

function save_comparison_figure(curves)
    fig = Figure(size = (600, 400), fontsize = 18)
    ax = Axis(
        fig[1, 1],
        xlabel = "SoftMix bandwidth b",
        ylabel = "Interaction energy (eV)",
        title = "Screened Hubbard interactions vs graphene cRPA",
    )

    colors = Makie.wong_colors()
    # xticks = unique(vcat((curve.bandwidths for curve in curves)...))
    # ax.xticks = sort!(collect(xticks))

    xlims!(ax, 0.0, 1.1)

    for (i, curve) in enumerate(curves)
        color = colors[mod1(i, length(colors))]
        lines!(ax, curve.bandwidths, curve.values; color = color, linewidth = 3, label = "$(curve.pair)")
        scatter!(ax, curve.bandwidths, curve.values; color = color, markersize = 10)
        hlines!(ax, [curve.crpa_ev]; color = color, linestyle = :dash, linewidth = 2)
    end

    axislegend(ax; position = :lt, framevisible = true, nbanks = 2)
    save(OUTPUT_FILE, fig)
    return fig
end

df = load_softmix_hubbard_results(INPUT_FILE)
curves = pair_curve_data(df)
fig = save_comparison_figure(curves)

fig
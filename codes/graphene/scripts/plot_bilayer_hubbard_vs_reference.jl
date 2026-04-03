module PlotBilayerHubbardVsReference

using CairoMakie
using CSV
using DataFrames

const INPUT_FILE = normpath(joinpath(@__DIR__, "..", "data", "screened_hubbard_graphene.csv"))
const OUTPUT_FILE = normpath(joinpath(@__DIR__, "..", "figs", "bilayer_hubbard_vs_reference.png"))
const PAIR_ORDER = ["U_00", "U_01", "U_02", "U_03", "U_04", "U_05"]
const BLG_TABLE_I_CRPA_EV = Dict(
    "U_00" => 9.1,
    "U_01" => 4.7,
    "U_02" => 3.2,
    "U_03" => 2.9,
    "U_04" => 2.5,
    "U_05" => 2.3,
)

function load_softmix_results(path::AbstractString = INPUT_FILE)
    df = CSV.read(path, DataFrame)
    return df[df.mode .== "SoftMixInversePermittivity", :]
end

function pair_curve_data(df::DataFrame; value_col::Symbol = :u_total_ev)
    curves = NamedTuple[]
    for pair in PAIR_ORDER
        rows = df[df.pair .== pair, [:bandwidth, value_col]]
        nrow(rows) == 0 && continue
        sort!(rows, :bandwidth)
        push!(
            curves,
            (
                pair = pair,
                bandwidths = collect(rows.bandwidth),
                values = Float64.(rows[!, value_col]),
                reference_ev = BLG_TABLE_I_CRPA_EV[pair],
            ),
        )
    end
    return curves
end

function save_comparison_figure(curves; output_file::AbstractString = OUTPUT_FILE)
    fig = Figure(size = (700, 450), fontsize = 18)
    ax = Axis(
        fig[1, 1],
        xlabel = "SoftMix bandwidth b",
        ylabel = "Interaction energy (eV)",
        title = "Bilayer-slab softmix vs BLG Table I cRPA",
    )

    colors = Makie.wong_colors()
    xlims!(ax, 0.0, 1.1)

    for (i, curve) in enumerate(curves)
        color = colors[mod1(i, length(colors))]
        lines!(ax, curve.bandwidths, curve.values; color = color, linewidth = 3, label = curve.pair)
        scatter!(ax, curve.bandwidths, curve.values; color = color, markersize = 10)
        hlines!(ax, [curve.reference_ev]; color = color, linestyle = :dash, linewidth = 2)
    end

    axislegend(ax; position = :rt, framevisible = true, nbanks = 2)
    save(output_file, fig)
    return fig
end

function main()
    df = load_softmix_results(INPUT_FILE)
    curves = pair_curve_data(df)
    return save_comparison_figure(curves)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end

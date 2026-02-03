using CSV, DataFrames
using CairoMakie

df = CSV.read(joinpath(@__DIR__, "data/compare_correction.csv"), DataFrame)

begin

    fig = Figure(size = (1000, 400), fontsize = 16)
    ax_1 = Axis(
        fig[1, 1],
        xlabel = "l_min",
        ylabel = "relative error of potential",
        xscale = log10,
        xreversed = true,
        yscale = log10,
        title = "near target",
    )
    ax_2 = Axis(
        fig[1, 2],
        xlabel = "l_min",
        ylabel = "relative error of potential",
        xscale = log10,
        xreversed = true,
        yscale = log10,
        title = "far target",
    )

    markers = [:circle, :diamond, :xcross, :rect, :rtriangle, :star4]
    colors = [:blue, :red, :green, :orange, :purple, :brown]


    for (i, dz) in enumerate([0.1, 1.0, 10.0])
        df_dz = df[df.dz .== dz .&& df.p .== 4, :]
        df_dz_ref = df[(df.dz .== dz) .&& (df.p .== 5), :]

        scatterlines!(ax_1, 1.0 ./ (2 .^ df_dz.r), abs.((df_dz.phi_nn .- df_dz_ref.phi_cn[1]) ./ df_dz_ref.phi_cn[1]), label = "non-corrected, dz = $dz", markersize = 10, marker = markers[1], color = colors[2 * i - 1])
        scatterlines!(ax_1, 1.0 ./ (2 .^ df_dz.r), abs.((df_dz.phi_cn .- df_dz_ref.phi_cn[1]) ./ df_dz_ref.phi_cn[1]), label = "corrected, dz = $dz", markersize = 10, marker = markers[2], color = colors[2 * i])

        scatterlines!(ax_2, 1.0 ./ (2 .^ df_dz.r), abs.((df_dz.phi_nf .- df_dz_ref.phi_cf[1]) ./ df_dz_ref.phi_cf[1]), label = "non-corrected, dz = $dz", markersize = 10, marker = markers[1], color = colors[2 * i - 1])
        scatterlines!(ax_2, 1.0 ./ (2 .^ df_dz.r), abs.((df_dz.phi_cf .- df_dz_ref.phi_cf[1]) ./ df_dz_ref.phi_cf[1]), label = "corrected, dz = $dz", markersize = 10, marker = markers[2], color = colors[2 * i])
    end

    Legend(fig[1, 3], ax_1, orientation = :vertical)

    fig
end
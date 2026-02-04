using CSV, DataFrames
using CairoMakie

df = CSV.read(joinpath(@__DIR__, "data/near_source_corrected_p.csv"), DataFrame)

ps = sort(unique(df.p))
rs = sort(unique(df.r))

ref_pot = -0.00862773606428543

begin

    fig = Figure(size = (500, 400), fontsize = 16)

    ax = Axis(
        fig[1, 1],
        xlabel = "l_min",
        xscale = log10,
        xreversed = true,
        ylabel = "relative error",
        yscale = log10,
    )

    for p in ps
        df_p = df[df.p .== p, :]
        scatter!(ax, df_p.l_ec, abs.((df_p.potential_cff .- ref_pot) ./ ref_pot), label = "p = $p", markersize = 10)
    end

    axislegend(ax, position = :lb)

    fig
end
using CSV, DataFrames
using CairoMakie

df = CSV.read(joinpath(@__DIR__, "data/near_eval_accuracy.csv"), DataFrame)

trg_zs = sort(unique(df.trg_z))

begin
    fig = Figure(size = (500, 500), fontsize = 16)

    ax = Axis(
        fig[1, 1],
        xlabel = "desired accuracy",
        ylabel = "relative error",
        xscale = log10,
        xreversed = true,
        yscale = log10,
        title = "Near evaluation accuracy",
    )

    markers = [:xcross, :circle, :rect, :rtriangle, :diamond]
    colors = [:green, :red, :blue, :orange, :purple]
    ms = 10

    for (k, trg_z) in enumerate(trg_zs)
        df_z = df[df.trg_z .== trg_z, :]
        scatter!(ax, df_z.tol, df_z.err, label = "dz = $(round(trg_z, digits=2))", markersize = ms, marker = markers[k], color = colors[k])
    end

    # x = y line
    lines!(ax, [1e-10, 1e0], [1e-10, 1e0], color = :gray, linestyle = :dash)

    Legend(fig[0, :], ax, orientation = :horizontal, nbanks = 2)

    fig
end

save(joinpath(@__DIR__, "figs/near_eval_accuracy.svg"), fig)

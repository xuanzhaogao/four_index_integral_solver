using CSV, DataFrames
using CairoMakie, LaTeXStrings

df = CSV.read("data/batch_near_eval.csv", DataFrame)

begin
    fig = Figure(size = (500, 400), fontsize = 16)
    ax = Axis(fig[1, 1], xlabel = "tolerance", ylabel = "l2 relative error", yscale = log10, xscale = log10, xreversed = true)

    for range in [4.0, 6.0, 8.0]
        df_range = df[df.range .== range, :]
        scatterlines!(ax, df_range.tol, df_range.l2_rel_err, label = "range = $range")
    end

    xlims!(1e-1, 1e-13)

    axislegend(position = :rt)

    fig

end

save("figs/batch_near_eval.svg", fig)
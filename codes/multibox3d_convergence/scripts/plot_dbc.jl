using CSV, DataFrames
using CairoMakie

df = CSV.read("data/doublebox_convergence.csv", DataFrame)

ps = unique(df.n_quad)
rs = unique(df.edge_refine)

begin
    fig = Figure(size = (500, 400), fontsize = 20)
    ax = Axis(fig[1, 1], yscale = log10, xscale = log10, xreversed = true, xlabel = "l_edge", ylabel = "L2 relative error", title = "Double cubes", )
    for p in ps
        df_p = filter(row -> row.n_quad == p, df)
        scatter!(ax, 1 ./ 2 .^ rs, df_p.l2_error, label = "p = $p")
    end
    ylims!(ax, 1e-3, 1e-1)
    axislegend(ax)

    save("figs/doublebox_convergence.svg", fig)

    fig
end
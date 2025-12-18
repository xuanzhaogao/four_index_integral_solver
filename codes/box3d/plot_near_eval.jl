using CSV, DataFrames, CairoMakie

df_direct = CSV.read(joinpath(@__DIR__, "data/real_density_near_eval_direct.csv"), DataFrame)
df_upsampling = CSV.read(joinpath(@__DIR__, "data/real_density_near_eval_upsampling.csv"), DataFrame)

begin
    fig = Figure(size = (1200, 500), fontsize = 20)
    ax1 = Axis(fig[1, 1], xlabel = "n", ylabel = "Error", yscale = log10, title = "Direct Evaluation")
    scatterlines!(ax1, df_direct.n, df_direct.err, marker = :circle, markersize = 12)

    ax2 = Axis(fig[1, 2], xlabel = "n", ylabel = "Error", yscale = log10, title = "Upsampling Evaluation")
    for n in unique(df_upsampling.n)
        df_n = filter(row -> row.n == n, df_upsampling)
        scatterlines!(ax2, df_n.n_upsampling, df_n.err, marker = :circle, markersize = 12, label = "$(n)")
    end
    axislegend(ax2, position = :lb)


    ylims!(ax1, 1e-16, 1e1)
    ylims!(ax2, 1e-16, 1e1)
    save(joinpath(@__DIR__, "figs/real_density_near_eval.svg"), fig)
    fig
end
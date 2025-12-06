using CSV, DataFrames, CairoMakie

df = CSV.read(joinpath(@__DIR__, "data/single_box3d_gi.csv"), DataFrame)

df_non_reduce = filter(row -> row.reduce_quad == 0, df)
df_reduce = filter(row -> row.reduce_quad == 1, df)

n_quads = [1, 4, 8, 12, 16]
n_edges = [0, 2, 4, 6, 8, 10, 12]

begin
    fig = Figure(size = (1200, 1500), fontsize = 20)
    ax1 = Axis(fig[1, 1], xlabel = "n_quad", ylabel = "average L1 Error", title = "non-reduced", yscale = log10)
    ax2 = Axis(fig[1, 2], xlabel = "n_quad", ylabel = "average L1 Error", title = "reduced", yscale = log10)
    ax3 = Axis(fig[2, 1], xlabel = "n_quad", ylabel = "average L2 Error", title = "non-reduced", yscale = log10)
    ax4 = Axis(fig[2, 2], xlabel = "n_quad", ylabel = "average L2 Error", title = "reduced", yscale = log10)
    ax5 = Axis(fig[3, 1], xlabel = "n_quad", ylabel = "Linf Error", title = "non-reduced", yscale = log10)
    ax6 = Axis(fig[3, 2], xlabel = "n_quad", ylabel = "Linf Error", title = "reduced", yscale = log10)

    for (ax, df) in [(ax1, df_non_reduce), (ax2, df_reduce)]
        for n_quad in n_quads
            df_quad = filter(row -> row.n_quad == n_quad, df)
            scatterlines!(ax, df_quad.n_edge, df_quad.L1_err ./ df_quad.n_val, marker = :circle, markersize = 12, label = "n_quad = $(n_quad)")
        end
        axislegend(ax, position = :rt)
    end

    for (ax, df) in [(ax3, df_non_reduce), (ax4, df_reduce)]
        for n_quad in n_quads
            df_quad = filter(row -> row.n_quad == n_quad, df)
            scatterlines!(ax, df_quad.n_edge, df_quad.L2_err ./ df_quad.n_val, marker = :circle, markersize = 12, label = "n_quad = $(n_quad)")
        end
        axislegend(ax, position = :rt)
    end
    
    for (ax, df) in [(ax5, df_non_reduce), (ax6, df_reduce)]
        for n_quad in n_quads
            df_quad = filter(row -> row.n_quad == n_quad, df)
            scatterlines!(ax, df_quad.n_edge, df_quad.Linf_err, marker = :circle, markersize = 12, label = "n_quad = $(n_quad)")
        end
        axislegend(ax, position = :rt)
    end

    fig
end

save(joinpath(@__DIR__, "figs/single_box3d_gi.svg"), fig)
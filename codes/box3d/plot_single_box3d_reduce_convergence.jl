using CSV, DataFrames, CairoMakie, LaTeXStrings

df = CSV.read(joinpath(@__DIR__, "data/single_box3d_reduce_convergence.csv"), DataFrame)

df_non_reduce = CSV.read(joinpath(@__DIR__, "data/single_box3d_convergence.csv"), DataFrame)
df_non_reduce_p2 = filter(row -> row.n_quad == 2, df_non_reduce)
df_non_reduce_p4 = filter(row -> row.n_quad == 4, df_non_reduce)


begin
    fig = Figure(size = (1000, 800), fontsize = 20)

    ax1 = Axis(fig[1, 1], xlabel = L"r", ylabel = "Relative Error", yscale = log10, title = "Potential p_min = 2")
    ax2 = Axis(fig[1, 2], xlabel = L"p_{\text{max}}", ylabel = "Relative Error", yscale = log10, title = "Potential p_min = 2")
    ax3 = Axis(fig[2, 1], xlabel = L"r", ylabel = "Relative Error", yscale = log10, title = "Potential p_min = 4")
    ax4 = Axis(fig[2, 2], xlabel = L"p_{\text{max}}", ylabel = "Relative Error", yscale = log10, title = "Potential p_min = 4")

    n_quads = unique(df.n_quad_max)
    n_edges = unique(df.n_edge)

    for (n_quad_min, ax_r, ax_p, df_non_reduce) in [(2, ax1, ax2, df_non_reduce_p2), (4, ax3, ax4, df_non_reduce_p4)]

        for n_quad_max in n_quads
            df_quad = filter(row -> row.n_quad_max == n_quad_max && row.n_quad_min == n_quad_min, df)
            scatterlines!(ax_r, df_quad.n_edge, df_quad.pot_relerr, marker = :circle, markersize = 12, label = "$(n_quad_max)")
        end
        
        lines!(ax_r, df_non_reduce.n_edge, df_non_reduce.pot_relerr, linestyle = :dash, label = "nr", color = :black, linewidth = 2)

        axislegend(ax_r, L"p_{\text{max}}", position = :rt, nbanks = 3)
        xlims!(ax_r, -1, 11)
        ylims!(ax_r, 1e-4, 1e2)

        for n_edge in n_edges
            df_edge = filter(row -> row.n_edge == n_edge && row.n_quad_min == n_quad_min, df)
            scatterlines!(ax_p, df_edge.n_quad_max, df_edge.pot_relerr, marker = :circle, markersize = 12, label = "$(n_edge)")
        end

        axislegend(ax_p, L"r", position = :lt, nbanks = 3)
        xlims!(ax_p, 3, 15)
        ylims!(ax_p, 1e-4, 1e2)
    end

    save(joinpath(@__DIR__, "figs/single_box3d_reduce_convergence.svg"), fig)
    fig
end
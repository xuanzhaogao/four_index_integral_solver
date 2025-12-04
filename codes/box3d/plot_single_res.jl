using CSV, DataFrames, CairoMakie

df = CSV.read(joinpath(@__DIR__, "data/single_box3d.csv"), DataFrame)

begin
    fig = Figure(size = (1200, 1000), fontsize = 20)

    # 1. n_quad vs pot_relerr
    ax1 = Axis(fig[1, 1], xlabel = "n_quad", ylabel = "Relative Error", yscale = log10, title = "Convergence vs. n_quad")
    df1 = filter(row -> row.n_boxes == 1 && row.n_edge == 0 && row.n_corner == 0, df)
    scatterlines!(ax1, df1.n_quad, df1.pot_relerr, marker = :circle, markersize = 12, label = "n_box = 1", color = :blue)

    # 2. n_box vs pot_relerr
    ax2 = Axis(fig[1, 2], xlabel = "n_box", ylabel = "Relative Error", yscale = log10, title = "Convergence vs. n_box")
    df2 = filter(row -> row.n_quad == 12 && row.n_edge == 0 && row.n_corner == 0, df)
    scatterlines!(ax2, df2.n_boxes, df2.pot_relerr, marker = :rect, markersize = 12, label = "n_quad = 12", color = :blue)

    # 3. n_edge vs pot_relerr
    ax3 = Axis(fig[2, 1], xlabel = "n_edge", ylabel = "Relative Error", yscale = log10, title = "Convergence vs. n_edge")
    df3 = filter(row -> row.n_boxes == 4 && row.n_quad == 12 && row.n_corner == 0, df)
    scatterlines!(ax3, df3.n_edge, df3.pot_relerr, marker = :diamond, markersize = 12, label = "n_box=4, n_quad=12", color = :blue)

    # 4. n_corner vs pot_relerr
    ax4 = Axis(fig[2, 2], xlabel = "n_corner", ylabel = "Relative Error", yscale = log10, title = "Convergence vs. n_corner")
    df4 = filter(row -> row.n_boxes == 4 && row.n_quad == 12 && row.n_edge == 6, df)
    scatterlines!(ax4, df4.n_corner, df4.pot_relerr, marker = :star5, markersize = 12, label = "n_box=4, n_quad=12, n_edge=6", color = :blue)

    axislegend(ax1, position = :lb)
    axislegend(ax2, position = :lb)
    axislegend(ax3, position = :lb)
    axislegend(ax4, position = :lb)
    ylims!(ax1, 1e-4, 1e1)
    ylims!(ax2, 1e-4, 1e1)
    ylims!(ax3, 1e-4, 1e1)
    ylims!(ax4, 1e-4, 1e1)

    fig
end

save(joinpath(@__DIR__, "figs/single_box3d_convergence.svg"), fig)
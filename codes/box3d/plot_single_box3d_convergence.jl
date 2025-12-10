using CSV, DataFrames, CairoMakie, LaTeXStrings

df = CSV.read(joinpath(@__DIR__, "data/single_box3d_convergence.csv"), DataFrame)

n_quads = sort(unique(df.n_quad))
n_edges = sort(unique(df.n_edge))

begin
    fig = Figure(size = (500, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = L"N_{\mathrm{edge}}", ylabel = "Relative Error", yscale = log10, title = "Convergence vs. p")

    for n_quad in n_quads
        df_quad = filter(row -> row.n_quad == n_quad, df)
        scatterlines!(ax, n_edges, df_quad.pot_relerr, marker = :circle, markersize = 12, label = "$(n_quad)")
    end

    axislegend(ax, L"p", position = :lb, nbanks = 2)
    ylims!(ax, 1e-4, 1e1)

    save(joinpath(@__DIR__, "figs/single_box3d_convergence.svg"), fig)
    fig
end
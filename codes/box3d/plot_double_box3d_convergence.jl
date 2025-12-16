using CSV, DataFrames, CairoMakie, LaTeXStrings

df = CSV.read(joinpath(@__DIR__, "data/double_box3d_convergence.csv"), DataFrame)

begin
    fig = Figure(size = (1000, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = L"r", ylabel = "Relative Error", yscale = log10, title = "Potential")
    ax_flux = Axis(fig[1, 2], xlabel = L"r", ylabel = "Relative Error", yscale = log10, title = "Flux")

    n_quads = [1, 2, 3, 4]
    for n_quad in n_quads
        df_nquad = filter(row -> row.n_quad == n_quad, df)
        n_edges = sort(unique(df_nquad.n_edge))
        pot_relerrs = []
        flux_relerrs = []
        for n_edge in n_edges
            df_edge = filter(row -> row.n_edge == n_edge, df_nquad)
            push!(pot_relerrs, df_edge.pot_relerr[])
            push!(flux_relerrs, abs(1 - df_edge.gi[]))
        end
        scatterlines!(ax, n_edges, pot_relerrs, marker = :circle, markersize = 12, label = "$(n_quad)")
        scatterlines!(ax_flux, n_edges, flux_relerrs, marker = :circle, markersize = 12, label = "$(n_quad)")
    end

    axislegend(ax, L"p", position = :lb, nbanks = 2)
    axislegend(ax_flux, L"p", position = :lb, nbanks = 2)
    ylims!(ax, 1e-5, 1e1)
    ylims!(ax_flux, 1e-5, 1e1)

    # ax2 = Axis(fig[2, 1], xlabel = L"p", ylabel = "Relative Error", yscale = log10, title = "Potential")
    # ax2_flux = Axis(fig[2, 2], xlabel = L"p", ylabel = "Relative Error", yscale = log10, title = "Flux")
    # for n_edge in [3, 4]
    #     pot_relerrs = []
    #     flux_relerrs = []
    #     df_edge = filter(row -> row.n_edge == n_edge, df)
    #     n_quads = sort(unique(df_edge.n_quad))
    #     for n_quad in n_quads
    #         df_quad = filter(row -> row.n_quad == n_quad, df_edge)
    #         push!(pot_relerrs, df_quad.pot_relerr[])
    #         push!(flux_relerrs, abs(1 - df_quad.gi[]))
    #     end
    #     scatterlines!(ax2, n_quads, pot_relerrs, marker = :circle, markersize = 12, label = "$(n_edge)")
    #     scatterlines!(ax2_flux, n_quads, flux_relerrs, marker = :circle, markersize = 12, label = "$(n_edge)")
    # end

    # axislegend(ax2, L"r", position = :lb, nbanks = 2)
    # axislegend(ax2_flux, L"r", position = :lb, nbanks = 2)
    # ylims!(ax2, 1e-4, 1e1)
    # ylims!(ax2_flux, 1e-4, 1e1)

    save(joinpath(@__DIR__, "figs/double_box3d_convergence.svg"), fig)
    fig
end
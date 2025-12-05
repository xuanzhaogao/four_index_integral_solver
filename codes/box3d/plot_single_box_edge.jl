using CSV, DataFrames, CairoMakie

df = CSV.read(joinpath(@__DIR__, "data/single_box3d_edges.csv"), DataFrame)

df_non_reduce = filter(row -> row.reduce_quad == 0, df)

n_quads = [1, 4, 8, 12, 16]
n_edges = [0, 2, 4, 6, 8]

pot_relerrs_non_reduce = []
for n_quad in n_quads
    df_quad = filter(row -> row.n_quad == n_quad, df_non_reduce)
    pot_relerrs_quad = []
    for n_edge in n_edges
        df_edge = filter(row -> row.n_edge == n_edge, df_quad)
        pot_relerr = df_edge.pot_relerr
        push!(pot_relerrs_quad, pot_relerr...)
    end
    push!(pot_relerrs_non_reduce, pot_relerrs_quad)
end

df_reduce = filter(row -> row.reduce_quad == 1, df)
pot_relerrs_reduce = []
for n_quad in n_quads
    df_quad = filter(row -> row.n_quad == n_quad, df_reduce)
    pot_relerrs_quad = []
    for n_edge in n_edges
        df_edge = filter(row -> row.n_edge == n_edge, df_quad)
        pot_relerr = df_edge.pot_relerr
        push!(pot_relerrs_quad, pot_relerr...)
    end
    push!(pot_relerrs_reduce, pot_relerrs_quad)
end


# with the same n_quad, how error changes with n_edge?
begin
    fig = Figure(size = (1000, 400))
    ax = Axis(fig[1, 1], xlabel = "n_edge", ylabel = "Relative Error", yscale = log10, title = "Convergence vs. n_edge (non-reduced)")
    ax2 = Axis(fig[1, 2], xlabel = "n_edge", ylabel = "Relative Error", yscale = log10, title = "Convergence vs. n_edge (reduced)")

    for i in 1:length(n_quads)
        scatterlines!(ax, n_edges, pot_relerrs_non_reduce[i], markersize = 12, label = "Nq = $(n_quads[i])")
    end
    axislegend(ax, position = :lb, framevisible = false)

    # for i in 1:length(n_quads)
    #     scatterlines!(ax2, n_edges, pot_relerrs_reduce[i], markersize = 12, label = "Nq = $(n_quads[i])")
    # end
    # axislegend(ax2, position = :lb, framevisible = false)

    ylims!(ax, 1e-4, 1e2)
    ylims!(ax2, 1e-4, 1e2)

    fig
end

pot_relerrs_non_reduce_rev = []
for n_edge in n_edges
    pot_relerrs_edge = []
    for n_quad in n_quads
        df_quad = filter(row -> row.n_quad == n_quad, df_non_reduce)
        df_edge = filter(row -> row.n_edge == n_edge, df_quad)
        pot_relerr = df_edge.pot_relerr
        push!(pot_relerrs_edge, pot_relerr...)
    end
    push!(pot_relerrs_non_reduce_rev, pot_relerrs_edge)
end

begin
    fig = Figure(size = (600, 400))
    ax = Axis(fig[1, 1], xlabel = "N_quad", ylabel = "Relative Error", title = "Relative Error vs. order of quadrature (non-reduced)", yscale = log10, xticks = n_quads)
    for i in 1:length(n_edges)
        scatterlines!(ax, n_quads, pot_relerrs_non_reduce_rev[i], markersize = 12, label = "n_edge = $(n_edges[i])")
    end
    Legend(fig[1, 2], ax, "N_edge")
    save(joinpath(@__DIR__, "figs/single_box3d_edges_non_reduce_nquad.svg"), fig)
    fig
end

gi_non_reduce = []
for n_quad in n_quads
    df_quad = filter(row -> row.n_quad == n_quad, df_non_reduce)
    gi_quad = []
    for n_edge in n_edges
        df_edge = filter(row -> row.n_edge == n_edge, df_quad)
        gi = df_edge.gi
        push!(gi_quad, gi...)
    end
    push!(gi_non_reduce, gi_quad)
end

begin
    fig = Figure(size = (500, 400))
    ax = Axis(fig[1, 1], xlabel = "n_edge", ylabel = "Total Flux", title = "GI vs. n_edge (non-reduced)", yscale = log10)
    for i in 1:length(n_quads)
        scatterlines!(ax, n_edges, abs.(1 .- gi_non_reduce[i]), markersize = 12, label = "Nq = $(n_quads[i])")
    end
    axislegend(ax, position = :lb, framevisible = false)
    fig
end

gi_non_reduce_rev = []
for n_edge in n_edges
    gi_edge = []
    for n_quad in n_quads
        df_quad = filter(row -> row.n_quad == n_quad, df_non_reduce)
        df_edge = filter(row -> row.n_edge == n_edge, df_quad)
        gi = df_edge.gi
        push!(gi_edge, gi...)
    end
    push!(gi_non_reduce_rev, gi_edge)
end

begin
    fig = Figure(size = (500, 400))
    ax = Axis(fig[1, 1], xlabel = "Nq", ylabel = "Total Flux", title = "GI vs. Nq (non-reduced)", yscale = log10)
    for i in 1:length(n_quads)
        scatterlines!(ax, n_quads, abs.(1 .- gi_non_reduce_rev[i]), markersize = 12, label = "n_edge = $(n_edges[i])")
    end
    axislegend(ax, position = :lb, framevisible = false)
    fig
end
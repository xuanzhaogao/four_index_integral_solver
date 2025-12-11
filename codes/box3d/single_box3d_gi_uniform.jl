include(joinpath(@__DIR__, "single_box3d_utils.jl"))

n_edges = [2, 4, 8]
n_quads = [2, 4, 8]

n_boxes = 1
eps = 4.0

pts = range(0.9, 1.0, length = 1000)
trgs = zeros(3, length(pts))
for i in 1:length(pts)
    trgs[:, i] = [pts[i], 0.0, 1.0]
end

gi_trgs_list = []
for n_edge in n_edges
    gi_t = []
    for n_quad in n_quads
        dbox = BI.dielectric_box3d(eps, 1.0, n_boxes, n_quad, n_quad, n_edge, n_edge)

        @show BI.num_points(dbox)

        D_trgs = BI.laplace3d_D_trg_fmm3d(dbox, trgs, 1e-14)
        gi_trgs = D_trgs * ones(BI.num_points(dbox))

        @show n_quad, n_edge
        push!(gi_t, gi_trgs)
    end
    push!(gi_trgs_list, gi_t)
end

trgs_to_corners = zeros(3, length(pts))
for i in 1:length(pts)
    trgs_to_corners[:, i] = [pts[i], pts[i], 1.0]
end

gi_trgs_to_corners_list = []
for n_edge in n_edges
    gi_t = []
    for n_quad in n_quads
        dbox = BI.dielectric_box3d(eps, 1.0, n_boxes, n_quad, n_quad, n_edge, n_edge)

        @show BI.num_points(dbox)

        D_trgs = BI.laplace3d_D_trg_fmm3d(dbox, trgs_to_corners, 1e-14)
        gi_trgs = D_trgs * ones(BI.num_points(dbox))

        @show n_quad, n_edge
        push!(gi_t, gi_trgs)
    end
    push!(gi_trgs_to_corners_list, gi_t)
end


using CairoMakie

begin
    fig = Figure(size = (1000, 800), fontsize = 20)
    ax_1 = Axis(fig[1, 1], xlabel = "x", ylabel = "value", title = "p = 8", xticks = 0.95:0.01:1.0, yticks = -0.55:0.05:-0.2) # fixed n_quad to 8
    ax_2 = Axis(fig[1, 2], xlabel = "x", ylabel = "error", title = "p = 8", yscale = log10, xticks = 0.9:0.02:1.0) 

    ax_3 = Axis(fig[2, 1], xlabel = "x", ylabel = "value", title = "r = 8", xticks = 0.95:0.01:1.0, yticks = -0.55:0.05:-0.2) # fixed n_adapt to 8
    ax_4 = Axis(fig[2, 2], xlabel = "x", ylabel = "error", title = "r = 8", yscale = log10, xticks = 0.9:0.02:1.0) # fixed n_adapt to 8

    for j in 1:3
        lines!(ax_1, pts, gi_trgs_list[3][j], label = "r = $(n_edges[j])")
        lines!(ax_2, pts, abs.(0.5 .+ gi_trgs_list[3][j]))
    end

    for i in 1:3
        lines!(ax_3, pts, gi_trgs_list[i][3], label = "p = $(n_quads[i])")
        lines!(ax_4, pts, abs.(0.5 .+ gi_trgs_list[i][3]))
    end

    xmin = 0.95
    xmax = 1.005

    xlims!(ax_1, xmin, xmax)
    xlims!(ax_2, xmin, xmax)
    xlims!(ax_3, xmin, xmax)
    xlims!(ax_4, xmin, xmax)
    
    ylims!(ax_1, -0.54, -0.2)
    ylims!(ax_2, 1e-16, 1.0)
    ylims!(ax_3, -0.55, -0.2)
    ylims!(ax_4, 1e-16, 1.0)

    axislegend(ax_1, position = :lt)
    axislegend(ax_3, position = :lt)

    save(joinpath(@__DIR__, "figs/single_box3d_gi_uniform_to_edge.svg"), fig)

    fig
end


begin
    fig = Figure(size = (1000, 800), fontsize = 20)
    ax_1 = Axis(fig[1, 1], xlabel = "x", ylabel = "value", title = "p = 8", xticks = 0.95:0.01:1.0, yticks = -0.625:0.125:0.0) # fixed n_quad to 8
    ax_2 = Axis(fig[1, 2], xlabel = "x", ylabel = "error", title = "p = 8", yscale = log10, xticks = 0.9:0.02:1.0) 

    ax_3 = Axis(fig[2, 1], xlabel = "x", ylabel = "value", title = "r = 8", xticks = 0.95:0.01:1.0, yticks = -0.625:0.125:0.0) # fixed n_adapt to 8
    ax_4 = Axis(fig[2, 2], xlabel = "x", ylabel = "error", title = "r = 8", yscale = log10, xticks = 0.9:0.02:1.0) # fixed n_adapt to 8

    for j in 1:3
        lines!(ax_1, pts, gi_trgs_to_corners_list[3][j], label = "r = $(n_edges[j])")
        lines!(ax_2, pts, abs.(0.5 .+ gi_trgs_to_corners_list[3][j]))
    end

    for i in 1:3
        lines!(ax_3, pts, gi_trgs_to_corners_list[i][3], label = "p = $(n_quads[i])")
        lines!(ax_4, pts, abs.(0.5 .+ gi_trgs_to_corners_list[i][3]))
    end

    xmin = 0.95
    xmax = 1.005

    xlims!(ax_1, xmin, xmax)
    xlims!(ax_2, xmin, xmax)
    xlims!(ax_3, xmin, xmax)
    xlims!(ax_4, xmin, xmax)
    
    ylims!(ax_1, -0.625, 0.0)
    ylims!(ax_2, 1e-16, 1.0)
    ylims!(ax_3, -0.625, 0.0)
    ylims!(ax_4, 1e-16, 1.0)

    axislegend(ax_1, position = :lt)
    axislegend(ax_3, position = :lt)

    save(joinpath(@__DIR__, "figs/single_box3d_gi_uniform_to_corners.svg"), fig)

    fig
end
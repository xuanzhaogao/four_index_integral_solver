using BoundaryIntegral
import BoundaryIntegral as BI
using OMEinsum, LinearAlgebra
using CairoMakie
using LsqFit

# use the Green's identity as a bar for error
# number of quadrature points fixed to 16

function res_nonadaptive(n_panels::Int, eps_box::Float64, src::Tuple{Float64, Float64})
    box = BI.dielectric_box2d(n_panels, 16, adapt = false)
    lhs = BI.Lhs_dielectric_box2d(eps_box, box)
    rhs = BI.Rhs_dielectric_box2d(eps_box, box, src)
    sigma = BI.solve_lu(lhs, rhs)
    return (box, sigma)
end

function res_adaptive(n_panels::Int, eps_box::Float64, n_adapt::Int, src::Tuple{Float64, Float64}; n_quad::Int = 16)
    box = BI.dielectric_box2d(n_panels, n_quad, adapt = true, n_adapt = n_adapt)
    lhs = BI.Lhs_dielectric_box2d(eps_box, box)
    rhs = BI.Rhs_dielectric_box2d(eps_box, box, src)
    sigma = BI.solve_lu(lhs, rhs)
    return (box, sigma)
end

function res_adaptive_mbox(n_panels::Int, n_adapt::Int, n_quads::Int, eps_boxes::Vector{Float64}, rects::Vector{Vector{Tuple{Float64, Float64}}}, src::Tuple{Float64, Float64}, eps_src::Float64)
    mbox = BI.dielectric_mbox2d(eps_boxes, rects, n_panels, n_quads, n_adapt)
    lhs = BI.Lhs_dielectric_mbox2d(mbox)
    rhs = BI.Rhs_dielectric_mbox2d(mbox, src, eps_src)
    sigma = BI.solve_lu(lhs, rhs)
    return (mbox, sigma)
end


function plot_contourf(xs, ys, zs, src, title)
    non_zero_zs = filter(!iszero, zs)
    color_lim = if isempty(non_zero_zs)
        1.0
    else
        max(abs.(extrema(non_zero_zs))...)
    end

    sym_levels = range(-color_lim, stop = color_lim, length = 50)

    fig = Figure(size = (500, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y", title = title, aspect = DataAspect())
    co = contourf!(ax, xs, ys, zs, levels = sym_levels, colormap = :RdBu)
    scatter!(ax, [src[1]], [src[2]], color = :red, markersize = 10)
    tightlimits!(ax)

    Colorbar(fig[1, 2], co)
    return fig
end

function plot_contourf_error(xs, ys, zs, zs_ref, src, title)
    non_zero_zs = filter(!iszero, zs)
    color_lim = if isempty(non_zero_zs)
        1.0
    else
        max(abs.(extrema(non_zero_zs))...)
    end

    sym_levels = range(-color_lim, stop = color_lim, length = 50)

    fig = Figure(size = (1000, 400), fontsize = 20, title = title)
    ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y", title = "value", aspect = DataAspect())
    co = contourf!(ax, xs, ys, zs, levels = sym_levels, colormap = :RdBu)
    scatter!(ax, [src[1]], [src[2]], color = :red, markersize = 10)
    lines!(ax, [0.0, 1.0], [1.0, 1.0], color = :black)
    lines!(ax, [1.0, 1.0], [0.0, 1.0], color = :black)
    tightlimits!(ax)

    ax2 = Axis(fig[1, 3], xlabel = "x", ylabel = "y", title = "log10 error", aspect = DataAspect())
    co2 = contourf!(ax2, xs, ys, log10.(abs.(zs - zs_ref)), levels = 50)
    lines!(ax2, [0.0, 1.0], [1.0, 1.0], color = :black)
    lines!(ax2, [1.0, 1.0], [0.0, 1.0], color = :black)
    tightlimits!(ax2)

    Colorbar(fig[1, 2], co)
    Colorbar(fig[1, 4], co2)
    return fig
end

function plot_contourf_mbox_error(xs, ys, zs, zs_ref, src, box)
    non_zero_zs = filter(!iszero, zs)
    color_lim = if isempty(non_zero_zs)
        1.0
    else
        max(abs.(extrema(non_zero_zs))...)
    end

    sym_levels = range(-color_lim, stop = color_lim, length = 50)

    fig = Figure(size = (1000, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y", title = "value", aspect = DataAspect())
    co = contourf!(ax, xs, ys, zs, levels = sym_levels, colormap = :RdBu)
    scatter!(ax, [src[1]], [src[2]], color = :red, markersize = 10)

    for interface in box.interfaces
        for panel in interface[1].panels
            lines!(ax, [panel.points[1][1], panel.points[end][1]], [panel.points[1][2], panel.points[end][2]], color = :blue)
        end
    end

    tightlimits!(ax)

    ax2 = Axis(fig[1, 3], xlabel = "x", ylabel = "y", title = "log10 error", aspect = DataAspect())
    co2 = contourf!(ax2, xs, ys, log10.(abs.(zs - zs_ref)), levels = 50)
    tightlimits!(ax2)

    for interface in box.interfaces
        for panel in interface[1].panels
            lines!(ax2, [panel.points[1][1], panel.points[end][1]], [panel.points[1][2], panel.points[end][2]], color = :blue)
        end
    end

    Colorbar(fig[1, 2], co)
    Colorbar(fig[1, 4], co2)
    return fig
end

function heatmap_mbox_error(xs, ys, zs, zs_ref, src, box)

    fig = Figure(size = (1000, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y", title = "value", aspect = DataAspect())
    co = heatmap!(ax, xs, ys, zs)
    scatter!(ax, [src[1]], [src[2]], color = :red, markersize = 10)

    for interface in box.interfaces
        for panel in interface[1].panels
            lines!(ax, [panel.points[1][1], panel.points[end][1]], [panel.points[1][2], panel.points[end][2]], color = :blue)
        end
    end

    tightlimits!(ax)

    ax2 = Axis(fig[1, 3], xlabel = "x", ylabel = "y", title = "log10 error", aspect = DataAspect())
    co2 = heatmap!(ax2, xs, ys, log10.(abs.(zs - zs_ref)))
    tightlimits!(ax2)

    for interface in box.interfaces
        for panel in interface[1].panels
            lines!(ax2, [panel.points[1][1], panel.points[end][1]], [panel.points[1][2], panel.points[end][2]], color = :blue)
        end
    end

    Colorbar(fig[1, 2], co)
    Colorbar(fig[1, 4], co2)
    return fig
end
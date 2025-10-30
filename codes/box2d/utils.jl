using BoundaryIntegral
import BoundaryIntegral as BI
using OMEinsum, LinearAlgebra
using CairoMakie
using LsqFit

# use the Green's identity as a bar for error
# number of quadrature points fixed to 16

function res_nonadaptive(n_panels::Int, eps_box::Float64, src::Tuple{Float64, Float64})
    box = BI.dielectric_box2d(n_panels, 16, false)
    lhs = BI.Lhs_dielectric_box2d(eps_box, box)
    rhs = BI.Rhs_dielectric_box2d(eps_box, box, src)
    sigma = BI.solve_lu(lhs, rhs)
    return (box, sigma)
end

function res_adaptive(n_panels::Int, eps_box::Float64, n_adapt::Int, src::Tuple{Float64, Float64})
    box = BI.dielectric_box2d(n_panels, 16, true, n_adapt = n_adapt)
    lhs = BI.Lhs_dielectric_box2d(eps_box, box)
    rhs = BI.Rhs_dielectric_box2d(eps_box, box, src)
    sigma = BI.solve_lu(lhs, rhs)
    return (box, sigma)
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
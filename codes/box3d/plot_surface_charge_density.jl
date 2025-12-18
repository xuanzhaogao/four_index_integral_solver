using JLD2, CairoMakie
using BoundaryIntegral
import BoundaryIntegral as BI

begin
    data = load(joinpath(@__DIR__, "data/surface_charge_density.jld2"))
    sigma_f = data["sigma_f"]
    points = data["points"]

    xs = unique([point[1] for point in points])
    ys = unique([point[2] for point in points])

    fs = reshape(sigma_f, length(xs), length(ys))

    fig = Figure(size = (1000, 1000), fontsize = 20)
    ax = Axis3(fig[1, 1], xlabel = "x", ylabel = "y", zlabel = "z", title = "Cubic σ")
    surface!(ax, xs, ys, fs, colormap = :viridis, colorrange = (0.0, 0.1))
    xlims!(ax, -0.99, 0.99)
    ylims!(ax, -0.99, 0.99)
    zlims!(ax, 0.0, 0.1)

    save(joinpath(@__DIR__, "figs/cubic_surface_charge_density.png"), fig)

    fig
end


begin
    data = load(joinpath(@__DIR__, "data/surface_charge_density_thin_L10_a0.25_p6_r8.jld2"))
    sigma_f = data["sigma_f"]
    points = data["points"]

    xs = unique([point[1] for point in points])
    ys = unique([point[2] for point in points])

    fs = reshape(sigma_f, length(xs), length(ys))

    fig = Figure(size = (1000, 1000), fontsize = 20)
    ax = Axis3(fig[1, 1], xlabel = "x", ylabel = "y", zlabel = "log10(σ)", title = "L = 10, nxy = 40, p = 6, r = 8, log10(σ)")
    surface!(ax, xs, ys, log10.(fs), colormap = :viridis, colorrange = (-3, 0))
    xlims!(ax, -4.99, 4.99)
    ylims!(ax, -4.99, 4.99)
    # zlims!(ax, 0.0, 0.02)

    save(joinpath(@__DIR__, "figs/thin_surface_charge_density.png"), fig)

    fig
end
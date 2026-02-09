using CSV, DataFrames
import BoundaryIntegral as BI
using CairoMakie, LaTeXStrings
using Roots

df = CSV.read(joinpath(@__DIR__, "data/divergence_rate.csv"), DataFrame)

epss = unique(df.eps_in)
xs = unique(df.x)

epss_theory = range(minimum(epss), maximum(epss), length = 100)
g_theory = zeros(100)
for (i, eps) in enumerate(epss_theory)
    g_theory[i] = fzero(g -> BI.theta_shooting_even(π/2, eps, g), 1.0) - 1.0
end

g_x0 = [df[df.eps_in .== eps .&& df.x .≈ 0.0, :slope][1] for eps in epss]

begin
    fig = Figure(size = (500, 400), fontsize = 16)
    ax = Axis(fig[1, 1], xlabel = L"\epsilon", ylabel = L"\alpha", title = "Divergence Rate of Edge Singularities, x = 0.0")

    lines!(ax, epss_theory, g_theory, label = "theorical", color = :blue)
    scatter!(ax, epss, g_x0, label = "numerical", color = :red)

    axislegend(ax, position = :rt)

    save(joinpath(@__DIR__, "figs/divergence_rate_eps.svg"), fig)

    fig
end

begin
    fig = Figure(size = (500, 400), fontsize = 16)
    ax = Axis(fig[1, 1], xlabel = L"x", ylabel = L"\alpha", title = "Divergence Rate of Edge Singularities, ϵ = 10")

    eps = 10.0
    df_eps = df[df.eps_in .== eps, :]
    scatter!(ax, df_eps.x, df_eps.slope, label = "numerical", color = :red)
    
    g_4 = fzero(g -> BI.theta_shooting_even(π/2, eps, g), 1.0) - 1.0
    lines!(ax, [-0.5, 0], [g_4, g_4], label = "theorical", color = :blue)

    axislegend(ax, position = :rt)

    ylims!(ax, -0.3, -0.2)

    save(joinpath(@__DIR__, "figs/divergence_rate_x.svg"), fig)

    fig
end
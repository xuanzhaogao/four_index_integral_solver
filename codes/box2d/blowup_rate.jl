include("utils.jl")
using Roots

src = (0.4, 0.6)
n_panels = 8
n_adapt = 30

function catch_rt_corner(box, sigma, n_points)
    xs_top = Float64[]
    sigmas_top = Float64[]
    i = 0
    for panel in box.panels
        for point in panel.points
            i += 1
            if point[1] == 1.0
                push!(xs_top, point[2])
                push!(sigmas_top, sigma[i])
            end
        end
    end

    sp = sortperm(xs_top)
    xs_top = xs_top[sp]
    sigmas_top = sigmas_top[sp]

    return (1.0 .- xs_top[end - n_points:end], sigmas_top[end - n_points:end])
end

function theta_shooting_even(al, e, g)
	return sin(g * al / 2) * cos(g * (π - al / 2)) + cos(g * al / 2) * sin(g * (π - al/2)) / e
end

function first_root_even(al, e)
    return fzero(g -> theta_shooting_even(al, e, g), 1.0)
end

gammas1 = [-0.9, -0.6, -0.3, 0.3, 0.6, 0.9]
xss = []
sigmas_s = []
for gamma in gammas1
    eps_box = - (gamma + 1.0) / (gamma - 1.0)
    box, sigma = res_adaptive(n_panels, eps_box, n_adapt, src)
    xs, sigmas = catch_rt_corner(box, sigma, 16 * 5)
    push!(xss, xs)
    push!(sigmas_s, sigmas)
end

gammas = [-0.99:0.02:0.99...]
k_fits = Float64[]
k_theory = [first_root_even(pi / 2,  - (gamma + 1.0) / (gamma - 1.0)) for gamma in gammas]

for gamma in gammas
    eps_box = - (gamma + 1.0) / (gamma - 1.0)
    @show gamma, eps_box
    box, sigma = res_adaptive(n_panels, eps_box, n_adapt, src)

    xs, sigmas = catch_rt_corner(box, sigma, 100)

    @. model(x, p) = p[1] * x + p[2]
    fit = curve_fit(model, log10.(xs), log10.(abs.(sigmas)), [1.0, 0.0])
    @show gamma, fit.param[1]

    push!(k_fits, fit.param[1])
end

begin
    fig = Figure(size = (1000, 400), fontsize = 20)
    ax1 = Axis(fig[1, 1], xlabel = "x", ylabel = L"\sigma")
    ax2 = Axis(fig[1, 2], xlabel = "x", ylabel = L"\sigma")
    scatter!(ax1, xss[1], sigmas_s[1], label = L"\gamma = -0.9")
    scatter!(ax1, xss[2], sigmas_s[2], label = L"\gamma = -0.6")
    scatter!(ax1, xss[3], sigmas_s[3], label = L"\gamma = -0.3")
    axislegend(ax1, position = :rt)
    ylims!(ax1, -0.025, 0.01)

    scatter!(ax2, xss[6], sigmas_s[6], label = L"\gamma = 0.9")
    scatter!(ax2, xss[5], sigmas_s[5], label = L"\gamma = 0.6")
    scatter!(ax2, xss[4], sigmas_s[4], label = L"\gamma = 0.3")
    ylims!(ax2, -10.0, 100.0)
    axislegend(ax2, position = :rt)
    save("figs/density.svg", fig)
    fig
end

begin
    fig = Figure(size = (500, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = "(eps - 1) / (eps + 1)", ylabel = "k")
    scatter!(ax, gammas, k_fits, label = "fitted k")

    lines!(ax, gammas, k_theory .- 1.0, label = "theory (even parity)", color = :red, linewidth = 2)

    hlines!(ax, [1/3], color = :black, label = L"+1/3")
    hlines!(ax, [-1/3], color = :black, label = L"-1/3")
    axislegend(ax, position = :rt)
    save("figs/blowup_rate.svg", fig)
    fig
end
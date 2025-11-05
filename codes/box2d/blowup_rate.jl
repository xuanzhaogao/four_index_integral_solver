include("utils.jl")

src = (0.5, 0.5)
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

# begin
#     gamma = -0.99
#     eps_box = - (gamma + 1.0) / (gamma - 1.0)

#     box, sigma = res_adaptive(n_panels, eps_box, n_adapt, src)
#     xs, sigmas = catch_rt_corner(box, sigma, 16 * 5)

#     fig = Figure(size = (500, 400), fontsize = 20)
#     ax = Axis(fig[1, 1], xlabel = "1-x", ylabel = "sigma")
#     scatter!(ax, xs, sigmas)
#     save("blowup/0.99_$(n_adapt).svg", fig)
#     fig
# end

gammas = [-0.99:0.02:0.99...]
k_fits = Float64[]

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
    fig = Figure(size = (500, 400), fontsize = 20)
    ax = Axis(fig[1, 1], xlabel = "(eps - 1) / (eps + 1)", ylabel = "k")
    scatter!(ax, gammas, k_fits, label = "fitted k")
    hlines!(ax, [1/3], color = :red, label = L"+1/3")
    hlines!(ax, [-1/3], color = :red, label = L"-1/3")
    axislegend(ax, position = :rt)
    save("blowup_rate.svg", fig)
    fig
end
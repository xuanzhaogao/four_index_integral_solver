include("utils.jl")

eps_box = 2.0
src = (0.5, 0.5)
n_panels = 8
n_adapt = 45

box1, sigma1 = res_adaptive(n_panels, 2.0, n_adapt, (0.5, 0.5))
box2, sigma2 = res_adaptive(n_panels, 2.0, n_adapt, (0.9, 0.9))
box3, sigma3 = res_adaptive(n_panels, 1.1, n_adapt, (0.9, 0.9))
box4, sigma4 = res_adaptive(n_panels, 0.5, n_adapt, (0.9, 0.9))
box5, sigma5 = res_adaptive(n_panels, 1000000.0, n_adapt, (0.99, 0.99))

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


np = 16 * 15

xs1, sigmas1 = catch_rt_corner(box1, sigma1, np)
xs2, sigmas2 = catch_rt_corner(box2, sigma2, np)
xs3, sigmas3 = catch_rt_corner(box3, sigma3, np)
xs4, sigmas4 = catch_rt_corner(box4, sigma4, np)
xs5, sigmas5 = catch_rt_corner(box5, sigma5, np)

@. model(x, p) = p[1] * x + p[2]
fit1 = curve_fit(model, log10.(xs1), log10.(abs.(sigmas1)), [1.0, 0.0])
fit2 = curve_fit(model, log10.(xs2), log10.(abs.(sigmas2)), [1.0, 0.0])
fit3 = curve_fit(model, log10.(xs3), log10.(abs.(sigmas3)), [1.0, 0.0])
fit4 = curve_fit(model, log10.(xs4), log10.(abs.(sigmas4)), [1.0, 0.0])
fit5 = curve_fit(model, log10.(xs5), log10.(abs.(sigmas5)), [1.0, 0.0])

@show fit1.param[1], fit2.param[1], fit3.param[1], fit4.param[1], fit5.param[1]

begin
    fig = Figure(size = (1000, 400), fontsize = 20)
    ax1 = Axis(fig[1, 1], xlabel = "1-x", ylabel = "sigma", xscale = log10)
    ax2 = Axis(fig[1, 2], xlabel = "1-x", ylabel = "abs(sigma)", xscale = log10, yscale = log10)

    scatter!(ax1, xs1, (sigmas1), label = "eps = 2.0, d = 0.5")
    scatter!(ax1, xs2, (sigmas2), label = "eps = 2.0, d = 0.1")
    scatter!(ax1, xs3, (sigmas3), label = "eps = 1.1, d = 0.1")
    scatter!(ax1, xs4, (sigmas4), label = "eps = 0.5, d = 0.1")
    # scatter!(ax1, xs5, (sigmas5), label = "eps = 1e6, d = 0.1")
    axislegend(ax1, position = :rt)

    scatter!(ax2, xs1, abs.(sigmas1))
    scatter!(ax2, xs2, abs.(sigmas2))
    scatter!(ax2, xs3, abs.(sigmas3))
    scatter!(ax2, xs4, abs.(sigmas4))
    # scatter!(ax2, xs5, abs.(sigmas5))

    lines!(ax2, xs1, exp10.(model(log10.(xs1), fit1.param)), label = "k = $(round(fit1.param[1], digits = 2))")
    lines!(ax2, xs2, exp10.(model(log10.(xs2), fit2.param)), label = "k = $(round(fit2.param[1], digits = 2))")
    lines!(ax2, xs3, exp10.(model(log10.(xs3), fit3.param)), label = "k = $(round(fit3.param[1], digits = 2))")
    lines!(ax2, xs4, exp10.(model(log10.(xs4), fit4.param)), label = "k = $(round(fit4.param[1], digits = 2))")
    # lines!(ax2, xs5, exp10.(model(log10.(xs5), fit5.param)), label = "k = $(round(fit5.param[1], digits = 2))")
    # ylims!(ax2, 10^(-1.2), 10^(1))
    axislegend(ax2, position = :rt)

    save("blowup.svg", fig)
    fig
end

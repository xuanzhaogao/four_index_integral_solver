include("utils.jl")

eps_box = 2.0

srcs = [(0.5, 0.5), (0.9, 0.9), (0.99, 0.99)]
n_panels_nonadaptive = 1:30

g_nonadaptive = []
for src in srcs
    gt = Float64[]
    for n_panels in n_panels_nonadaptive
        box, sigma = res_nonadaptive(n_panels, eps_box, src)
        g = BI.l2d_singlelayer_gi(box, sigma, 2.0, 128) + 1.0 / eps_box
        push!(gt, g)
        @show n_panels, 1.0 - g
    end
    push!(g_nonadaptive, gt)
end

n_adapts = 1:40
g_adaptive = []
for src in srcs
    gt = Float64[]
    for n_adapt in n_adapts
        box, sigma = res_adaptive(8, eps_box, n_adapt, src)
        g = BI.l2d_singlelayer_gi(box, sigma, 2.0, 128) + 1.0 / eps_box
        push!(gt, g)
        @show n_adapt, 1.0 - g
    end
    push!(g_adaptive, gt)
end

n_quads = 2:2:32
g_adaptive_quad = []
for src in srcs
    gt = Float64[]
    for n_quad in n_quads
        box, sigma = res_adaptive(8, eps_box, 30, src, n_quad = n_quad)
        g = BI.l2d_singlelayer_gi(box, sigma, 2.0, 128) + 1.0 / eps_box
        push!(gt, g)
        @show n_quad, 1.0 - g
    end
    push!(g_adaptive_quad, gt)
end

begin
    fig = Figure(size = (500, 400), fontsize = 20)
    ax1 = Axis(fig[1, 1], yscale = log10, xlabel = "Number of panels", ylabel = "Error of Green's identity", title = "Non-adaptive")
    ds = [0.5, 0.1, 0.01]
    for i in eachindex(srcs)
        lines!(ax1, n_panels_nonadaptive .* 4, abs.(1.0 .- g_nonadaptive[i]), label = "d = $(ds[i])")
    end
    axislegend(ax1, position = :rt)
    save("convergence_nonadaptive.svg", fig)
end

begin
    fig = Figure(size = (1000, 400), fontsize = 20)
    ax1 = Axis(fig[1, 1], yscale = log10, xlabel = "Number of quadrature points", ylabel = "Error of Green's identity")
    ax2 = Axis(fig[1, 2], yscale = log10, xlabel = "Number of adaptions", ylabel = "Error of Green's identity")
    for i in eachindex(srcs)
        scatter!(ax1, n_quads, abs.(1.0 .- g_adaptive_quad[i]), label = "d = $(ds[i])")
        lines!(ax1, n_quads, abs.(1.0 .- g_adaptive_quad[i]))
        scatter!(ax2, n_adapts, abs.(1.0 .- g_adaptive[i]), label = "d = $(ds[i])")
        lines!(ax2, n_adapts, abs.(1.0 .- g_adaptive[i]))
    end
    xlims!(ax1, 0, 35)
    axislegend(ax2, position = :rt)
    save("convergence_adaptive.svg", fig)

    fig
end
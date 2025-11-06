include("utils.jl")

eps_box1 = 2.0
eps_box2 = 3.0
eps_box3 = 4.0
eps_boxes = [eps_box1, eps_box2, eps_box3]

srcs = [(0.0, 0.5), (0.4, 0.9), (0.49, 0.99)]
ds = [0.5, 0.1, 0.01]
rects = [BI.square(-1.0, -1.0), BI.square(0.0, -1.0), BI.square(-0.5, 0.0)]
eps_src = eps_box3

n_adapts = 1:40
g_adaptive = []
for src in srcs
    gt = zeros(length(n_adapts))
    Threads.@threads for i in eachindex(n_adapts)
        n_adapt = n_adapts[i]
        box, sigma = res_adaptive_mbox(8, n_adapt, 16, eps_boxes, rects, src, eps_src)
        g = BI.l2d_singlelayer_gi(box, sigma, 2.0, 128) + 1.0 / eps_src
        gt[i] = g
        @show n_adapt, 1.0 - g
    end
    push!(g_adaptive, gt)
end

n_quads = 2:2:32
g_adaptive_quad = []
for src in srcs
    gt = zeros(Float64, length(n_quads))
    Threads.@threads for i in eachindex(n_quads)
        n_quad = n_quads[i]
        box, sigma = res_adaptive_mbox(8, 20, n_quad, eps_boxes, rects, src, eps_src)
        g = BI.l2d_singlelayer_gi(box, sigma, 2.0, 128) + 1.0 / eps_src
        gt[i] = g
        @show n_quad, 1.0 - g
    end
    push!(g_adaptive_quad, gt)
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
    ylims!(ax1, 1e-16, 1.0)
    ylims!(ax2, 1e-16, 1.0)
    axislegend(ax2, position = :rt)
    save("convergence_mbox.svg", fig)

    fig
end
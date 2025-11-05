include("utils.jl")

function plot_contourf_error(src, gamma, n_panels, n_adapt)
    eps_box = - (gamma + 1.0) / (gamma - 1.0)

    xs = range(0.0, 1.5, 500)
    ys = range(0.0, 1.5, 500)
    zs_ref = zeros(length(xs), length(ys))
    zs = zeros(length(xs), length(ys))
    
    @info "Computing reference solution"
    box_ref, sigma_ref = res_adaptive(n_panels, eps_box, 30, src)
    @info "Reference solution computed"

    @info "Computing reference potential"
    Threads.@threads for i in eachindex(xs)
        for j in eachindex(ys)
            zs_ref[i, j] = BI.laplace2d_singlelayer_interface(box_ref, sigma_ref, (xs[i], ys[j]))
        end
    end
    @info "Reference potential computed"

    @info "Computing solution"
    box, sigma = res_adaptive(n_panels, eps_box, n_adapt, src)
    @info "Solution computed"

    @info "Computing solution potential"
    Threads.@threads for i in eachindex(xs)
        for j in eachindex(ys)
            zs[i, j] = BI.laplace2d_singlelayer_interface(box, sigma, (xs[i], ys[j]))
        end
    end
    @info "Solution potential computed"

    fig1 = plot_contourf_error(xs, ys, zs, zs_ref, src, "gamma = $(gamma), n_adapt = $(n_adapt)")

    save("contourf/contourf_$(gamma)_$(n_adapt).svg", fig1)

    return fig1
end

for gamma in -0.95:0.1:0.95
    plot_contourf_error((0.9, 0.9), gamma, 8, 20)
end

plot_contourf_error((0.9, 0.9), 0.999, 8, 20)
# plot_contourf_error((0.9, 0.9), -0.999, 8, 20)
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

    save("figs/contourf/contourf_$(gamma)_$(n_adapt).svg", fig1)

    return fig1
end

for gamma in -0.95:0.1:0.95
    plot_contourf_error((0.9, 0.9), gamma, 8, 20)
end

plot_contourf_error((0.9, 0.9), 0.999, 8, 20)
# plot_contourf_error((0.9, 0.9), -0.999, 8, 20)

function plot_contourf_error!(sf_s, sf_bar1, sf_e, sf_bar2, src, gamma, n_panels, n_adapt, title)
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

    non_zero_zs = filter(!iszero, zs)
    color_lim = if isempty(non_zero_zs)
        1.0
    else
        max(abs.(extrema(non_zero_zs))...)
    end

    sym_levels = range(-color_lim, stop = color_lim, length = 50)

    ax_solution = Axis(sf_s,  xlabel = "x", ylabel = "y", title = title, aspect = DataAspect())
    co = contourf!(ax_solution, xs, ys, zs, levels = sym_levels, colormap = :RdBu)
    scatter!(ax_solution, [src[1]], [src[2]], color = :red, markersize = 10)
    lines!(ax_solution, [0.0, 1.0], [1.0, 1.0], color = :black)
    lines!(ax_solution, [1.0, 1.0], [0.0, 1.0], color = :black)
    tightlimits!(ax_solution)
    Colorbar(sf_bar1, co)

    ax_error = Axis(sf_e,  xlabel = "x", ylabel = "y", title = title, aspect = DataAspect())
    co2 = contourf!(ax_error, xs, ys, log10.(abs.(zs - zs_ref)), levels = 50)
    scatter!(ax_error, [src[1]], [src[2]], color = :red, markersize = 10)
    lines!(ax_error, [0.0, 1.0], [1.0, 1.0], color = :black)
    lines!(ax_error, [1.0, 1.0], [0.0, 1.0], color = :black)
    tightlimits!(ax_error)
    Colorbar(sf_bar2, co2)

    return nothing
end

begin
    fig_solution = Figure(size = (1500, 800), fontsize = 20)
    fig_error = Figure(size = (1500, 800), fontsize = 20)

    plot_contourf_error!(fig_solution[1, 1], fig_solution[1, 2], fig_error[1, 1], fig_error[1, 2], (0.7, 0.8), -0.9, 8, 20, L"\gamma = -0.9, n_{adapt} = 20")
    plot_contourf_error!(fig_solution[1, 3], fig_solution[1, 4], fig_error[1, 3], fig_error[1, 4], (0.7, 0.8), -0.6, 8, 20, L"\gamma = -0.6, n_{adapt} = 20")
    plot_contourf_error!(fig_solution[1, 5], fig_solution[1, 6], fig_error[1, 5], fig_error[1, 6], (0.7, 0.8), -0.3, 8, 20, L"\gamma = -0.3, n_{adapt} = 20")
    plot_contourf_error!(fig_solution[2, 1], fig_solution[2, 2], fig_error[2, 1], fig_error[2, 2], (0.7, 0.8), 0.3, 8, 20, L"\gamma = 0.3, n_{adapt} = 20")
    plot_contourf_error!(fig_solution[2, 3], fig_solution[2, 4], fig_error[2, 3], fig_error[2, 4], (0.7, 0.8), 0.6, 8, 20, L"\gamma = 0.6, n_{adapt} = 20")
    plot_contourf_error!(fig_solution[2, 5], fig_solution[2, 6], fig_error[2, 5], fig_error[2, 6], (0.7, 0.8), 0.9, 8, 20, L"\gamma = 0.9, n_{adapt} = 20")
end

save("figs/contourf/contourf_solution.png", fig_solution, px_per_unit = 2)
save("figs/contourf/contourf_error.png", fig_error, px_per_unit = 2)
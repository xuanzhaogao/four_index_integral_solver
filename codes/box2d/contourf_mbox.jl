include("utils.jl")

function calc_contourf_mbox(src, n_panels, n_adapt, rects, eps_boxes, eps_src)

    xs = range(-1.5, 1.5, 500)
    ys = range(-1.5, 1.5, 500)
    zs_ref = zeros(length(xs), length(ys))
    zs = zeros(length(xs), length(ys))
    
    @info "Computing reference solution"
    box_ref, sigma_ref = res_adaptive_mbox(n_panels, 30, 16, eps_boxes, rects, src, eps_src)
    @info "Reference solution computed"

    @info "Computing reference potential"
    for i in eachindex(xs)
        Threads.@threads for j in eachindex(ys)
            zs_ref[i, j] = BI.laplace2d_singlelayer_interface(box_ref, sigma_ref, (xs[i], ys[j]))
        end
    end
    @info "Reference potential computed"

    @info "Computing solution"
    box, sigma = res_adaptive_mbox(n_panels, n_adapt, 16, eps_boxes, rects, src, eps_src)
    @info "Solution computed"

    @info "Computing solution potential"
    for i in eachindex(xs)
        Threads.@threads for j in eachindex(ys)
            zs[i, j] = BI.laplace2d_singlelayer_interface(box, sigma, (xs[i], ys[j]))
        end
    end
    @info "Solution potential computed"

    return (xs, ys, zs, zs_ref, box_ref)
end

rects = [BI.square(-1.0, -1.0), BI.square(0.0, -1.0), BI.square(-0.5, 0.0)]
eps_boxes = [2.0, 3.0, 4.0]
eps_src = 4.0

xs, ys, zs, zs_ref, box_ref = calc_contourf_mbox((0.4, 0.9), 8, 20, rects, eps_boxes, eps_src)
fig = plot_contourf_mbox_error(xs, ys, zs, zs_ref, (0.4, 0.9), box_ref)
save("contourf_mbox/contourf_1.svg", fig)


xs, ys, zs, zs_ref, box_ref = calc_contourf_mbox((0.4, 0.9), 8, 20, rects, [10000.0, 3.0, 4.0], eps_src)
fig = plot_contourf_mbox_error(xs, ys, zs, zs_ref, (0.4, 0.9), box_ref)
save("contourf_mbox/contourf_2.svg", fig)

fig = BI.viz_2d_dielectric_interfaces(box_ref)
save("contourf_mbox/box_1.svg", fig)

include("utils.jl")

using IterativeSolvers

function potential_mbox_fmm2d(src, n_panels, n_adapt, n_quads, rects, eps_boxes, eps_src)
    xs = range(-1.5, 1.5, 500)
    ys = range(-1.5, 1.5, 500)

    targets = zeros(2, length(xs) * length(ys))
    for i in eachindex(xs)
        for j in eachindex(ys)
            targets[1, (i - 1) * length(ys) + j] = xs[i]
            targets[2, (i - 1) * length(ys) + j] = ys[j]
        end
    end

    mbox = BI.dielectric_mbox2d(eps_boxes, rects, n_panels, n_quads, n_adapt)
    lhs = BI.Lhs_dielectric_mbox2d_fmm2d(mbox)
    rhs = BI.Rhs_dielectric_mbox2d(mbox, src, eps_src)

    sigma, history = gmres(lhs, rhs, reltol = 1e-8, log = true, verbose = true, restart = 100)

    poteval = BI.laplace2d_pottarg_fmm2d(mbox, targets, 1e-8)
    potentials = poteval * sigma

    return transpose(reshape(potentials, length(xs), length(ys))), mbox
end

xs = range(-1.5, 1.5, 500)
ys = range(-1.5, 1.5, 500)
rects = [BI.square(-1.0, -1.0), BI.square(0.0, -1.0), BI.square(-0.5, 0.0)]
eps_src = 4.0

zs, mbox = potential_mbox_fmm2d((0.4, 0.9), 8, 20, 16, rects, [2.0, 3.0, 4.0], eps_src)
zs_ref, mbox_ref = potential_mbox_fmm2d((0.4, 0.9), 8, 30, 16, rects, [2.0, 3.0, 4.0], eps_src)

fig = heatmap_mbox_error(xs, ys, zs, zs_ref, (0.4, 0.9), mbox)
save("figs/mbox_fmm2d_1.png", fig)


zs, mbox = potential_mbox_fmm2d((0.4, 0.9), 8, 20, 16, rects, [20.0, 3.0, 4.0], eps_src)
zs_ref, mbox_ref = potential_mbox_fmm2d((0.4, 0.9), 8, 30, 16, rects, [20.0, 3.0, 4.0], eps_src)

fig = heatmap_mbox_error(xs, ys, zs, zs_ref, (0.4, 0.9), mbox)
save("figs/mbox_fmm2d_2.png", fig)


zs, mbox = potential_mbox_fmm2d((0.4, 0.9), 8, 20, 16, rects, [100.0, 3.0, 4.0], eps_src)
zs_ref, mbox_ref = potential_mbox_fmm2d((0.4, 0.9), 8, 30, 16, rects, [100.0, 3.0, 4.0], eps_src)

fig = heatmap_mbox_error(xs, ys, zs, zs_ref, (0.4, 0.9), mbox)
save("figs/mbox_fmm2d_3.png", fig)

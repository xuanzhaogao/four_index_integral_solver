include("utils.jl")

eps_box = 2.0

src_1 = (0.5, 0.5)
n_panels = 8
n_adapt = 20

box, sigma = res_adaptive(n_panels, eps_box, n_adapt, src_1)

xs = [0.0:0.005:1.5...]
ys = [0.0:0.005:1.5...]
zs = zeros(length(xs), length(ys))
for i in eachindex(xs)
    for j in eachindex(ys)
        xs[i] < 1.0 && ys[j] < 1.0 && continue
        zs[i, j] = BI.laplace2d_singlelayer_surface(box, sigma, (xs[i], ys[j]))
    end
end
fig1 = plot_contourf(xs, ys, zs, src_1, "single layer d = 0.5")
save("contourf_1.svg", fig1)

src_2 = (0.9, 0.9)
box2, sigma2 = res_adaptive(n_panels, eps_box, n_adapt, src_2)
zs2 = zeros(length(xs), length(ys))
for i in eachindex(xs)
    for j in eachindex(ys)
        xs[i] < 1.0 && ys[j] < 1.0 && continue
        zs2[i, j] = BI.laplace2d_singlelayer_surface(box2, sigma2, (xs[i], ys[j]))
    end
end
fig2 = plot_contourf(xs, ys, zs2, src_2, "single layer d = 0.9")
save("contourf_2.svg", fig2)
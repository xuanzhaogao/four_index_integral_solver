include("transform.jl")
include("laplace.jl")

using FastGaussQuadrature
using Plots

function integrand_3d(density, xt, yt, zt, nxt, nyt, nzt)
    f = x -> density(x[1], x[2]) * laplace3d_DT(xt - x[1], yt - x[2], zt, nxt, nyt, nzt)
    return f
end

function integrate_3d(fi, nx, ny)
    x, wx = gausslegendre(nx)
    y, wy = gausslegendre(ny)
    t = 0.0
    for (xi, wi) in zip(x, wx)
        for (yi, wyi) in zip(y, wy)
            t += wi * wyi * fi((xi, yi))
        end
    end
    return t
end

Nx, Ny = 16, 16
Nxu, Nyu = 128, 128

x, wx = gausslegendre(Nx)
y, wy = gausslegendre(Ny)
xu, wxu = gausslegendre(Nxu)
yu, wyu = gausslegendre(Nyu)

density = (x, y) -> y * cos(-x) + x * y^2 * exp(- (y-0.5)^2)

fvals = reshape([density(xi, yi) for xi in x, yi in y], Nx * Ny)

C = legendre_transform_matrix_2d(x, wx, y, wy; Mx=Nx, My=Ny)
cvals = C * fvals

E = legendre_eval_matrix_2d(xu, yu; Mx=Nx, My=Ny)
res = E * cvals

fvals_u = reshape([density(xi, yi) for xi in xu, yi in yu], Nxu * Nyu)

heatmap(xu, yu, reshape(res .- fvals_u, Nxu, Nyu))

heatmap(xu, yu, reshape(log10.(abs.(res .- fvals_u)), Nxu, Nyu))

z0 = 0.1
nx, ny, nz = 0.0, 0.0, 1.0

N_trgx, N_trgy = 100, 100
x_trgs = range(-2.0, 2.0, length=N_trgx)
y_trgs = range(-2.0, 2.0, length=N_trgy)

f_mat = zeros(N_trgx * N_trgy, Nxu * Nyu)
for i in 1:N_trgx, j in 1:N_trgy
    for k in 1:Nxu, l in 1:Nyu
        f_mat[i + (j-1)*N_trgx, k + (l-1)*Nxu] = laplace3d_DT(x_trgs[i] - xu[k], y_trgs[j] - yu[l], z0, nx, ny, nz) * wyu[l] * wxu[k]
    end
end

ops_mat = f_mat * E * C
res_trgs = ops_mat * fvals

res_mat = reshape(res_trgs, N_trgx, N_trgy)

res_low = zeros(N_trgx, N_trgy)
for i in 1:N_trgx, j in 1:N_trgy
    fi = integrand_3d(density, x_trgs[i], y_trgs[j], z0, nx, ny, nz)
    res_low[i, j] = integrate_3d(fi, 32, 32)
end

res_high_128 = zeros(N_trgx, N_trgy)
for i in 1:N_trgx, j in 1:N_trgy
    fi = integrand_3d(density, x_trgs[i], y_trgs[j], z0, nx, ny, nz)
    res_high_128[i, j] = integrate_3d(fi, 128, 128)
end

heatmap(x_trgs, y_trgs, log10.(abs.(res_low .- res_mat)))
fig_2d_err = heatmap(x_trgs, y_trgs, log10.(abs.(res_high_128 .- res_mat)), aspect_ratio=:equal, dpi=300)
xlims!(fig_2d_err, -2.0, 2.0)
ylims!(fig_2d_err, -2.0, 2.0)

savefig("2d_err.png")

# shifted evaluation on another grid

Ns = 16
Ne = 100
x, w = gausslegendre(Ns)

Lxa, Lxb = 1.1, 1.5
Lya, Lyb = 2.2, 2.3

xs = (Lxa + Lxb) / 2 .+ ((Lxb - Lxa) / 2) .* x
ys = (Lya + Lyb) / 2 .+ ((Lyb - Lya) / 2) .* y

wx = (Lxb - Lxa) / 2 .* w
wy = (Lyb - Lya) / 2 .* w

xe = range(Lxa, Lxb, length=Ne)
ye = range(Lya, Lyb, length=Ne)

xeo = (xe .- (Lxa + Lxb) / 2) ./ (Lxb - Lxa) .* 2
yeo = (ye .- (Lya + Lyb) / 2) ./ (Lyb - Lya) .* 2

fvals = reshape([density(xi, yi) for xi in xs, yi in ys], Ns * Ns)

C = legendre_transform_matrix_2d(x, w, y, w; Mx=Ns, My=Ns)
E = legendre_eval_matrix_2d(xeo, yeo; Mx=Ns, My=Ns)

eval_vals = (E * C) * fvals
f_real = reshape([density(xi, yi) for xi in xe, yi in ye], Ne * Ne)

heatmap(xe, ye, reshape(log10.(abs.(eval_vals .- f_real)), Ne, Ne), aspect_ratio=:equal)
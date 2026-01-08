include("transform.jl")
include("laplace.jl")

using FastGaussQuadrature
using LegendrePolynomials
using LinearAlgebra

using Plots

N = 32
x, w = gausslegendre(N)

density = x -> x^3 * exp(-x)

fvals = @. density(x)   # example function values

M = N    # max degree of Legendre expansion
c = values_to_legendre_coeffs(x, w, fvals; M=M)

eval_legendre_series(x, c) .- fvals

xs = -1:0.001:1

plot(xs, eval_legendre_series(xs, c))
plot!(xs, density.(xs))

plot(xs, log10.(abs.((eval_legendre_series(xs, c) .- density.(xs)) ./ density.(xs))))


P_mat = legendre_eval_matrix(xs, N)
mat_vals = P_mat * c

plot(xs, mat_vals)
plot!(xs, density.(xs))

plot(log10.(abs.((mat_vals .- density.(xs)))))

function integrand(density, xt, yt, nxt, nyt)
    f = x -> density(x[1]) * laplace2d_DT(xt - x[1], yt, nxt, nyt)
    return f
end

fi = integrand(density, 0.5, 0.01, 0.0, 1.0)

function integrate(fi, n)
    x, w = gausslegendre(n)
    t = 0.0
    for (xi, wi) in zip(x, w)
        t += wi * fi(xi)
    end
    return t
end

ints = Float64[]
for n in 1:2048
    push!(ints, integrate(fi, n))
end

plot(ints)
plot(log10.(abs.(ints .- ints[end])))

truth = ints[end]

function fake_density(density, N)
    x, w = gausslegendre(N)
    fvals = @. density(x)   # example function values

    c = values_to_legendre_coeffs(x, w, fvals; M=N)

    return x -> eval_legendre_series(x, c)
end

fdi = fake_density(density, 16)

fints = Float64[]
for n in 1:2048
    push!(fints, integrate(integrand(fdi, 0.5, 0.01, 0.0, 1.0), n))
end

plot(fints)
plot(log10.(abs.(fints .- truth)))

x, w = gausslegendre(16)
xl, wl = gausslegendre(2048)

T = legendre_transform_matrix(x, w; M=16)

P1 = legendre_eval_matrix(xl, 16)
x0, y0 = 0.5, 0.01
nx, ny = 0.0, 1.0

f_vec = [laplace2d_DT(x0 - xi, y0, nx, ny) * wi for (xi, wi) in zip(xl, wl)]

ops = f_vec' * P1 * T

res = ops * density.(x)

res - truth

trgs_x = range(-2.0, 2.0, length=1024)
trgs_y = 0.01

f_mat = zeros(length(trgs_x), length(xl))
for i in 1:length(trgs_x)
    for j in 1:length(xl)
        f_mat[i, j] = laplace2d_DT(trgs_x[i] - xl[j], trgs_y, nx, ny) * wl[j]
    end
end

ops_mat = f_mat * P1 * T
res_vec = ops_mat * density.(x)

truth_vec = [integrate(integrand(density, trgs_x[i], trgs_y, nx, ny), 2048) for i in 1:length(trgs_x)]

truth_vec - res_vec



plot(trgs_x, truth_vec)
plot!(trgs_x, res_vec)
plot(trgs_x, log10.(abs.(truth_vec .- res_vec)))
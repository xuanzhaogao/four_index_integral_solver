#!/usr/bin/env julia

using Random
using LinearAlgebra
using Statistics
using BoxDMK
using FBCPoisson

Random.seed!(20260302)

function gauss_legendre_nodes_weights(n::Int, a::Float64, b::Float64)
    n >= 2 || throw(ArgumentError("n must be >= 2"))
    beta = [k / sqrt(4k^2 - 1) for k in 1:(n - 1)]
    eig = eigen(SymTridiagonal(zeros(n), beta))
    x = ((b - a) / 2) .* collect(eig.values) .+ (a + b) / 2
    w = ((b - a) / 2) .* collect(2 .* (eig.vectors[1, :]).^2)
    return x, w
end

function make_targets_10x10x10(a::Float64, b::Float64)
    xs = collect(range(a, b; length = 10))
    nt = 10^3
    targets = Matrix{Float64}(undef, 3, nt)
    idx = 1
    for x in xs, y in xs, z in xs
        targets[:, idx] .= (x, y, z)
        idx += 1
    end
    return targets
end

function main()
    boxlen = 1.18
    sigma = 0.12
    center = (0.12, -0.08, 0.17)

    norder = 6
    eps = 1e-5

    nq = 24
    nfft = 48
    tol = 1e-8

    normc = 1 / ((2 * pi)^(3 / 2) * sigma^3)
    gaussian_density(x, y, z) = normc * exp(-((x - center[1])^2 + (y - center[2])^2 + (z - center[3])^2) / (2 * sigma^2))
    rnorm(x, y, z) = sqrt(x^2 + y^2 + z^2)
    rho(x, y, z) = gaussian_density(x, y, z) * sin(rnorm(x, y, z))

    density_cb(x, _) = rho(x[1], x[2], x[3])
    targets = make_targets_10x10x10(-0.5, 0.5)

    prob = LaplaceProblem(density = density_cb, nd = 1, ndim = 3, boxlen = boxlen)
    opts = BDMKOptions(eps = eps, norder = norder)
    tree, _ = solve_problem(prob; compute = :potential, opts = opts)
    p_box = vec(evaluate_targets(prob, tree, targets; compute = :potential, eps = eps).pote)

    xs, wx = gauss_legendre_nodes_weights(nq, -boxlen / 2, boxlen / 2)
    ys, wy = gauss_legendre_nodes_weights(nq, -boxlen / 2, boxlen / 2)
    zs, wz = gauss_legendre_nodes_weights(nq, -boxlen / 2, boxlen / 2)

    ns = nq^3
    sources = Matrix{Float64}(undef, 3, ns)
    charges = Vector{Float64}(undef, ns)
    idx = 1
    for i in eachindex(xs), j in eachindex(ys), k in eachindex(zs)
        x = xs[i]
        y = ys[j]
        z = zs[k]
        sources[:, idx] .= (x, y, z)
        charges[idx] = rho(x, y, z) * wx[i] * wy[j] * wz[k]
        idx += 1
    end

    p_fbc = FBCPoisson.lfbc3d(nfft, sources, charges, targets, tol, 1)
    p_fbc_4pi = 4 * pi .* p_fbc

    diff_raw = p_box - p_fbc
    diff_scaled = p_box - p_fbc_4pi

    println("density = gaussian_density * sin(sqrt(x^2 + y^2 + z^2))")
    println("targets=$(size(targets, 2)) sources=$ns")
    println("scaled4pi: rel_l2=$(norm(diff_scaled) / norm(p_fbc_4pi)) max_abs=$(maximum(abs.(diff_scaled))) rmse=$(sqrt(mean(abs2, diff_scaled)))")
end

main()

# first I want see how close will make the direct evaluation of rhs fails

using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
using CairoMakie
using BoundaryIntegral.FMM3D

function Rhs_grid_fmm3d(
    pts,
    normal,
    vs::VolumeSource{Float64, 3},
    eps_src::Float64,
    thresh::Float64,)
    xs, ys, zs = vs.axes
    weights = vs.weights
    density = vs.density
    nx, ny, nz = length(xs), length(ys), length(zs)
    n_sources = nx * ny * nz
    sources = zeros(Float64, 3, n_sources)
    charges = zeros(Float64, n_sources)
    idx = 0
    for ix in 1:nx, iy in 1:ny, iz in 1:nz
        idx += 1
        sources[1, idx] = xs[ix]
        sources[2, idx] = ys[iy]
        sources[3, idx] = zs[iz]
        charges[idx] = weights[ix, iy, iz] * density[ix, iy, iz]
    end

    n_points = length(pts)
    targets = zeros(Float64, 3, n_points)
    normals = zeros(Float64, 3, n_points)
    for (i, pt) in enumerate(pts)
        targets[1, i] = pt[1]
        targets[2, i] = pt[2]
        targets[3, i] = pt[3]
        normals[1, i] = normal[1]
        normals[2, i] = normal[2]
        normals[3, i] = normal[3]
    end

    vals = lfmm3d(thresh, sources, charges = charges, targets = targets, pgt = 2)
    grad = vals.gradtarg
    Rhs = zeros(Float64, n_points)
    for i in 1:n_points
        Rhs[i] = dot(normals[:, i], grad[:, i]) / (4π * eps_src)
    end
    return Rhs
end


function different_rhs(z_g, N)
    L = 10.0

    xs = range(-L / 2, stop = L / 2, length = N)
    ys = range(-L / 2, stop = L / 2, length = N)
    z = 0.0

    normal = (0.0, 0.0, 1.0)

    sigma = 1.0

    center = (0.0, 0.0, z_g)
    vs = BI.GaussianVolumeSource(center, sigma, 1e-8)

    rhs_exact = [dot(BI.gaussian_laplace3d_grad(center, (x, y, z), sigma), normal) for x in xs, y in ys]

    pts = Vector{NTuple{3, Float64}}()
    for x in xs, y in ys
        push!(pts, (x, y, 0.0))
    end

    rhs_fmm = Rhs_grid_fmm3d(pts, normal, vs, 1.0, 1e-12)
    rhs_fmm = reshape(rhs_fmm, N, N)

    return rhs_exact, rhs_fmm
end

begin
    fig = Figure(size = (1000, 400), fontsize = 20)
    ax = Axis(fig[1, 1], title = "exact value, z = 0.1, σ = 1.0", xlabel = "x", ylabel = "y", aspect = DataAspect())

    rhs_exact, rhs_fmm = different_rhs(0.1, 400)

    xs = range(-10.0 / 2, stop = 10.0 / 2, length = 400)
    ys = range(-10.0 / 2, stop = 10.0 / 2, length = 400)

    hm_1 = heatmap!(ax, xs, ys, log10.(abs.(rhs_exact)))
    Colorbar(fig[1, 2], hm_1)

    ax2 = Axis(fig[1, 3], title = "error", xlabel = "x", ylabel = "y", aspect = DataAspect())
    hm_2 = heatmap!(ax2, xs, ys, log10.(abs.(rhs_fmm .- rhs_exact)))
    Colorbar(fig[1, 4], hm_2)

    save("figs/rhs_error_z_0.1.png", fig)

    fig
end

for z in 1.0:3.0:10.0
    fig = Figure(size = (1000, 400), fontsize = 20)
    ax = Axis(fig[1, 1], title = "exact value, z = $(z), σ = 1.0", xlabel = "x", ylabel = "y", aspect = DataAspect())

    rhs_exact, rhs_fmm = different_rhs(z, 200)
    xs = range(-10.0 / 2, stop = 10.0 / 2, length = 200)
    ys = range(-10.0 / 2, stop = 10.0 / 2, length = 200)

    hm_1 = heatmap!(ax, xs, ys, log10.(abs.(rhs_exact)))
    Colorbar(fig[1, 2], hm_1)

    ax2 = Axis(fig[1, 3], title = "error", xlabel = "x", ylabel = "y", aspect = DataAspect())
    hm_2 = heatmap!(ax2, xs, ys, log10.(abs.(rhs_fmm .- rhs_exact)))
    Colorbar(fig[1, 4], hm_2)

    save("figs/rhs_error_z_$(z).png", fig)

    fig
end

errors = []
for z in 0.1:0.1:10.0
    rhs_exact, rhs_fmm = different_rhs(z, 10)
    error = norm(rhs_fmm .- rhs_exact, Inf) / maximum(abs.(rhs_exact))
    @show z, error
    push!(errors, error)
end

begin
    fig = Figure(size = (500, 400), fontsize = 20)
    ax = Axis(fig[1, 1], title = "σ = 1.0, fmm_tol = 1e-8", xlabel = "z", ylabel = "relative error", yscale = log10)
    scatterlines!(ax, 0.1:0.1:10.0, errors, color = :blue)
    save("figs/rhs_error_vs_z.png", fig)
    fig
end
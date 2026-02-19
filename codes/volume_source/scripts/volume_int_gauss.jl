# first I want see how close will make the direct evaluation of rhs fails

using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
using CairoMakie
using MPI
ENV["PVFMM_ALLOW_REVISE"] = "1"
push!(LOAD_PATH, expanduser("~/codes/pvfmm/julia"))
using PVFMM

function _pvfmm_comm_self()
    if !MPI.Initialized()
        MPI.Init(finalize_atexit = false)
        MPI.add_finalize_hook!(() -> GC.gc(true))
    end
    return MPI.COMM_SELF
end

const _PVFMM_VS_REF = Ref{Any}(nothing)
const _PVFMM_ORIGIN_REF = Ref{NTuple{3, Float64}}((0.0, 0.0, 0.0))
const _PVFMM_SIDE_REF = Ref{Float64}(1.0)

function _pvfmm_density_callback(coord::Ptr{Cdouble}, n::Clong, out::Ptr{Cdouble}, ::Ptr{Cvoid})::Cvoid
    coords = unsafe_wrap(Vector{Float64}, coord, 3 * Int(n))
    vals = unsafe_wrap(Vector{Float64}, out, Int(n)) # data_dim = 1
    vs = _PVFMM_VS_REF[]::VolumeSource{Float64, 3}
    origin = _PVFMM_ORIGIN_REF[]
    side = _PVFMM_SIDE_REF[]
    scale = side^2

    @inbounds for i in 1:Int(n)
        base = 3 * (i - 1)
        ux = coords[base + 1]
        uy = coords[base + 2]
        uz = coords[base + 3]
        x = origin[1] + side * ux
        y = origin[2] + side * uy
        z = origin[3] + side * uz
        vals[i] = scale * _volume_density_trilinear(vs, x, y, z)
    end
    return
end

const _PVFMM_DENSITY_CFUN = @cfunction(
    _pvfmm_density_callback,
    Cvoid,
    (Ptr{Cdouble}, Clong, Ptr{Cdouble}, Ptr{Cvoid}),
)

@inline function _flatten_targets(pts::Vector{NTuple{3, Float64}})
    trg = Vector{Float64}(undef, 3 * length(pts))
    for (i, p) in enumerate(pts)
        base = 3 * (i - 1)
        trg[base + 1] = p[1]
        trg[base + 2] = p[2]
        trg[base + 3] = p[3]
    end
    return trg
end

@inline function _volume_density_trilinear(vs::VolumeSource{Float64, 3}, x::Float64, y::Float64, z::Float64)
    xs, ys, zs = vs.axes
    nx, ny, nz = length(xs), length(ys), length(zs)

    if x < xs[1] || x > xs[end] || y < ys[1] || y > ys[end] || z < zs[1] || z > zs[end]
        return 0.0
    end

    ix = clamp(searchsortedlast(xs, x), 1, nx - 1)
    iy = clamp(searchsortedlast(ys, y), 1, ny - 1)
    iz = clamp(searchsortedlast(zs, z), 1, nz - 1)

    x0, x1 = xs[ix], xs[ix + 1]
    y0, y1 = ys[iy], ys[iy + 1]
    z0, z1 = zs[iz], zs[iz + 1]
    tx = (x - x0) / (x1 - x0)
    ty = (y - y0) / (y1 - y0)
    tz = (z - z0) / (z1 - z0)

    d = vs.density
    c000 = d[ix, iy, iz]
    c100 = d[ix + 1, iy, iz]
    c010 = d[ix, iy + 1, iz]
    c110 = d[ix + 1, iy + 1, iz]
    c001 = d[ix, iy, iz + 1]
    c101 = d[ix + 1, iy, iz + 1]
    c011 = d[ix, iy + 1, iz + 1]
    c111 = d[ix + 1, iy + 1, iz + 1]

    c00 = (1 - tx) * c000 + tx * c100
    c10 = (1 - tx) * c010 + tx * c110
    c01 = (1 - tx) * c001 + tx * c101
    c11 = (1 - tx) * c011 + tx * c111
    c0 = (1 - ty) * c00 + ty * c10
    c1 = (1 - ty) * c01 + ty * c11
    return (1 - tz) * c0 + tz * c1
end

function _tree_region(pts::Vector{NTuple{3, Float64}}, vs::VolumeSource{Float64, 3})
    xs, ys, zs = vs.axes
    src_min = (minimum(xs), minimum(ys), minimum(zs))
    src_max = (maximum(xs), maximum(ys), maximum(zs))

    trg_x = map(p -> p[1], pts)
    trg_y = map(p -> p[2], pts)
    trg_z = map(p -> p[3], pts)
    trg_min = (minimum(trg_x), minimum(trg_y), minimum(trg_z))
    trg_max = (maximum(trg_x), maximum(trg_y), maximum(trg_z))

    lower = (
        min(src_min[1], trg_min[1]),
        min(src_min[2], trg_min[2]),
        min(src_min[3], trg_min[3]),
    )
    upper = (
        max(src_max[1], trg_max[1]),
        max(src_max[2], trg_max[2]),
        max(src_max[3], trg_max[3]),
    )

    span = (
        upper[1] - lower[1],
        upper[2] - lower[2],
        upper[3] - lower[3],
    )
    side = 1.05 * maximum(span)
    cx = 0.5 * (lower[1] + upper[1])
    cy = 0.5 * (lower[2] + upper[2])
    cz = 0.5 * (lower[3] + upper[3])
    origin = (cx - 0.5 * side, cy - 0.5 * side, cz - 0.5 * side)
    return origin, side
end

function _uniform_leaves(depth::Int)
    leaves = Vector{NTuple{4, Float64}}()
    n = 1 << depth
    h = 1.0 / n
    for iz in 0:(n - 1), iy in 0:(n - 1), ix in 0:(n - 1)
        push!(leaves, (ix * h, iy * h, iz * h, h))
    end
    return leaves
end

function _select_uniform_depth(vs::VolumeSource{Float64, 3}, origin::NTuple{3, Float64}, side::Float64, tol::Float64)
    for depth in 2:8
        n = 1 << depth
        h = 1.0 / n
        max_rel_err = 0.0
        for kz in 0:(n - 1), ky in 0:(n - 1), kx in 0:(n - 1)
            ux = (kx + 0.5) * h
            uy = (ky + 0.5) * h
            uz = (kz + 0.5) * h
            x = origin[1] + side * ux
            y = origin[2] + side * uy
            z = origin[3] + side * uz
            c = _volume_density_trilinear(vs, x, y, z)

            c8 = 0.0
            for oz in (0.0, h), oy in (0.0, h), ox in (0.0, h)
                xc = origin[1] + side * (kx * h + ox)
                yc = origin[2] + side * (ky * h + oy)
                zc = origin[3] + side * (kz * h + oz)
                c8 += _volume_density_trilinear(vs, xc, yc, zc)
            end
            c8 *= 0.125
            rel = abs(c - c8) / max(abs(c), 1e-14)
            max_rel_err = max(max_rel_err, rel)
        end
        if max_rel_err <= tol
            return depth
        end
    end
    return 8
end

function _cheb_nodes_1d(cheb_deg::Int)
    n = cheb_deg + 1
    ξ = Vector{Float64}(undef, n)
    for j in 0:cheb_deg
        ξ[j + 1] = 0.5 * (1.0 - cos(pi * (2j + 1) / (2n)))
    end
    return ξ
end

function _leaf_coefficients(
    vs::VolumeSource{Float64, 3},
    leaves::Vector{NTuple{4, Float64}},
    cheb_deg::Int,
    origin::NTuple{3, Float64},
    side::Float64,
)
    ξ = _cheb_nodes_1d(cheb_deg)
    nleaf = length(leaves)
    ncheb = (cheb_deg + 1)^3
    node_val = Vector{Float64}(undef, nleaf * ncheb) # dof = 1
    leaf_coord = Vector{Float64}(undef, 3 * nleaf)

    for (leaf_i, leaf) in enumerate(leaves)
        lx, ly, lz, h = leaf
        leaf_base = 3 * (leaf_i - 1)
        leaf_coord[leaf_base + 1] = lx
        leaf_coord[leaf_base + 2] = ly
        leaf_coord[leaf_base + 3] = lz

        local_idx = 0
        for jz in 1:length(ξ), jy in 1:length(ξ), jx in 1:length(ξ)
            local_idx += 1
            ux = lx + h * ξ[jx]
            uy = ly + h * ξ[jy]
            uz = lz + h * ξ[jz]
            x = origin[1] + side * ux
            y = origin[2] + side * uy
            z = origin[3] + side * uz
            node_val[(leaf_i - 1) * ncheb + local_idx] = _volume_density_trilinear(vs, x, y, z)
        end
    end

    coeff = PVFMM.nodes_to_coeff(nleaf, cheb_deg, 1, node_val)
    return leaf_coord, coeff
end

function Rhs_grid_pvfmm(
    pts,
    normal,
    vs::VolumeSource{Float64, 3},
    eps_src::Float64,
    tol::Float64,
)
    origin, side = _tree_region(pts, vs)

    targets_u = Vector{Float64}(undef, 3 * length(pts))
    for (i, pt) in enumerate(pts)
        base = 3 * (i - 1)
        targets_u[base + 1] = (pt[1] - origin[1]) / side
        targets_u[base + 2] = (pt[2] - origin[2]) / side
        targets_u[base + 3] = (pt[3] - origin[3]) / side
    end

    cheb_deg = 4
    multipole_order = 8
    comm = _pvfmm_comm_self()
    _PVFMM_VS_REF[] = vs
    _PVFMM_ORIGIN_REF[] = origin
    _PVFMM_SIDE_REF[] = side

    tree = PVFMM.from_function(
        PVFMM.FMMVolumeTree{Float64},
        cheb_deg,
        1,
        _PVFMM_DENSITY_CFUN,
        C_NULL,
        targets_u,
        comm,
        tol,
        150,
        false,
        0,
    )
    ctx = PVFMM.FMMVolumeContext(multipole_order, cheb_deg, PVFMM.LaplaceGradient, comm; T = Float64)
    grad_u = PVFMM.evaluate(tree, ctx, length(pts))

    Rhs = zeros(Float64, length(pts))
    for i in eachindex(pts)
        base = 3 * (i - 1)
        grad_x = (
            grad_u[base + 1] / side,
            grad_u[base + 2] / side,
            grad_u[base + 3] / side,
        )
        Rhs[i] = (normal[1] * grad_x[1] + normal[2] * grad_x[2] + normal[3] * grad_x[3]) / (4π * eps_src)
    end
    return - 4π .* Rhs
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

    rhs_fmm = Rhs_grid_pvfmm(pts, normal, vs, 1.0, 1e-6)
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

# errors = []
# for z in 0.1:0.1:10.0
#     rhs_exact, rhs_fmm = different_rhs(z, 10)
#     error = norm(rhs_fmm .- rhs_exact, Inf) / maximum(abs.(rhs_exact))
#     @show z, error
#     push!(errors, error)
# end

# begin
#     fig = Figure(size = (500, 400), fontsize = 20)
#     ax = Axis(fig[1, 1], title = "σ = 1.0, fmm_tol = 1e-8", xlabel = "z", ylabel = "relative error", yscale = log10)
#     scatterlines!(ax, 0.1:0.1:10.0, errors, color = :blue)
#     save("figs/rhs_error_vs_z.png", fig)
#     fig
# end

using BoundaryIntegral
import BoundaryIntegral as BI
using FastGaussQuadrature
using JLD2

using BoundaryIntegral: legendre_transform_matrix_2d, legendre_eval_matrix_2d

function laplace3d_DT(x, y, z, nx, ny, nz)
    return (nx * x + ny * y + nz * z) / (sqrt(x^2 + y^2 + z^2))^3
end

# a linear operator that integrates the density on the original tgl grid at z = 0 to another tgl grid on the opposite side at z = 2
function laplace3d_f2f_op(n_origin, n_upsampling)
    x, w = gausslegendre(n_origin)
    xu, wu = gausslegendre(n_upsampling)
    C = legendre_transform_matrix_2d(x, w, x, w; Mx=n_origin, My=n_origin)
    E = legendre_eval_matrix_2d(xu, xu; Mx=n_origin, My=n_origin)

    f_eval = zeros(n_origin^2, n_upsampling^2)
    for i in 1:n_origin, j in 1:n_origin
        for k in 1:n_upsampling, l in 1:n_upsampling
            f_eval[i + (j-1)*n_origin, k + (l-1)*n_upsampling] = laplace3d_DT(xu[k] - x[i], xu[l] - x[j], 2.0, 0.0, 0.0, 1.0) * wu[l] * wu[k]
        end
    end

    return f_eval * E * C
end

function target_interface(n_target)
    x, w = gausslegendre(n_target)
    points = Vector{NTuple{3, Float64}}(undef, n_target^2)
    x_start = 1.0
    x_end = 2.0
    x_scaled = (x_end + x_start) / 2 .+ (x_end - x_start) / 2 .* x
    for i in 1:n_target, j in 1:n_target
        points[(i - 1) * n_target + j] = (x_scaled[i], x_scaled[j], - 0.5)
    end
    norm = (0.0, 0.0, 1.0)
    weights = ones(n_target^2)
    corners = [(-1.0, -1.0, -0.5), (1.0, -1.0, -0.5), (1.0, 1.0, -0.5), (-1.0, 1.0, -0.5)]
    panels = BI.Panel(n_target^2, points, norm, weights, corners)
    interface = BI.Interface(1, [panels])
    target_interface = BI.DielectricInterfaces(1, [(interface, 4.0, 1.0)])
    return target_interface
end

res_ref = load(joinpath(@__DIR__, "cache/single_thin_box_convergence_L10_r8_p6_r8.jld2"))
tbox = res_ref["tbox"]
sigma = res_ref["sigma"]

direct_eval_res = []
upsampling_res = []

x_opposite = range(-1.0, 1.0, length = 20)

for nt in [4, 8, 16, 32, 64, 128]
    ti = target_interface(nt)
    sigma_f = BI.nystrom_interpolation_dielectric_box3d(tbox, ti, (0.2, 0.3, 0.4), 4.0, sigma, 1e-6)
    sigma_f = reshape(sigma_f, nt, nt)
    @show "sigma_f is calculated"

    res = zeros(length(x_opposite), length(x_opposite))
    x, w = gausslegendre(nt)
    for i in 1:length(x_opposite)
        Threads.@threads for j in 1:length(x_opposite)
            for k in 1:nt, l in 1:nt
                res[i, j] += w[k] * w[l] * sigma_f[k, l] * laplace3d_DT(x_opposite[i] - x[k], x_opposite[j] - x[l], 2.0, 0.0, 0.0, 1.0)
            end
        end
    end

    push!(direct_eval_res, res)
end
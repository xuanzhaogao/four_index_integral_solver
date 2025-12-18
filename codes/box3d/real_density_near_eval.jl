using BoundaryIntegral
import BoundaryIntegral as BI
using FastGaussQuadrature
using LinearAlgebra
using JLD2

using BoundaryIntegral: legendre_transform_matrix_2d, legendre_eval_matrix_2d

function laplace3d_DT(x, y, z, nx, ny, nz)
    return (nx * x + ny * y + nz * z) / (sqrt(x^2 + y^2 + z^2))^3
end

# a linear operator that integrates the density on the original tgl grid at z = 0 to another tgl grid on the opposite side at z = 2
function laplace3d_f2f_op(n_origin, n_upsampling, xt, yt)
    x, w = gausslegendre(n_origin)
    xu, wu = gausslegendre(n_upsampling)
    C = legendre_transform_matrix_2d(x, w, x, w; Mx=n_origin, My=n_origin)
    E = legendre_eval_matrix_2d(xu, xu; Mx=n_origin, My=n_origin)

    f_eval = zeros(length(xt) * length(yt), n_upsampling^2)
    for i in 1:length(xt), j in 1:length(yt)
        for k in 1:n_upsampling, l in 1:n_upsampling
            f_eval[i + (j-1)*length(xt), k + (l-1)*n_upsampling] = laplace3d_DT(xt[i] - xu[k], yt[j] - xu[l], 0.2, 0.0, 0.0, 1.0) * wu[l] * wu[k]
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

x_opposite = range(-1.0, 1.0, length = 20)


sigma_f_list = []

for nt in [4, 8, 16, 32, 64, 128]
    ti = target_interface(nt)
    sigma_f = BI.nystrom_interpolation_dielectric_box3d(tbox, ti, (0.2, 0.3, 0.4), 4.0, sigma, 1e-6)
    sigma_f = reshape(sigma_f, nt, nt)
    @show "sigma_f is calculated"

    push!(sigma_f_list, sigma_f)
end

direct_eval_res = []
for (i, nt) in enumerate([4, 8, 16, 32, 64, 128])
    sigma_f = sigma_f_list[i]
    res = zeros(length(x_opposite), length(x_opposite))
    x, w = gausslegendre(nt)
    for i in 1:length(x_opposite)
        Threads.@threads for j in 1:length(x_opposite)
            for k in 1:nt, l in 1:nt
                res[i, j] += w[k] * w[l] * sigma_f[k, l] * laplace3d_DT(x_opposite[i] - x[k], x_opposite[j] - x[l], 0.2, 0.0, 0.0, 1.0)
            end
        end
    end

    push!(direct_eval_res, res)
end

direct_err = [LinearAlgebra.norm(direct_eval_res[end] - direct_eval_res[i], Inf) for i in 1:length(direct_eval_res) - 1]

ref = direct_eval_res[end]

upsampling_res = []
for (i, nt) in enumerate([4, 8, 16])
    ress = []
    for n_upsampling in [16, 32, 64, 128]

        sigma_f = sigma_f_list[i]
        ops = laplace3d_f2f_op(nt, n_upsampling, x_opposite, x_opposite)
        res = ops * reshape(sigma_f, nt^2)

        res = reshape(res, length(x_opposite), length(x_opposite))
        err = LinearAlgebra.norm(res - ref, Inf)
        @show "error = $err"

        push!(ress, res)
    end

    push!(upsampling_res, ress)
end

using CSV, DataFrames

df_direct = joinpath(@__DIR__, "data/real_density_near_eval_direct.csv")
df_upsampling = joinpath(@__DIR__, "data/real_density_near_eval_upsampling.csv")

CSV.write(df_direct, DataFrame(n = [], err = []))
CSV.write(df_upsampling, DataFrame(n = [], n_upsampling = [], err = []))

for (i, nt) in enumerate([4, 8, 16, 32, 64, 128])
    CSV.write(df_direct, DataFrame(n = [nt], err = [LinearAlgebra.norm(direct_eval_res[i] - ref, Inf)]), append = true)
end

for (i, nt) in enumerate([4, 8, 16])
    for (j, nu) in enumerate([16, 32, 64, 128])
        CSV.write(df_upsampling, DataFrame(n = [nt], n_upsampling = [nu], err = [LinearAlgebra.norm(upsampling_res[i][j] - ref, Inf)]), append = true)
    end
end
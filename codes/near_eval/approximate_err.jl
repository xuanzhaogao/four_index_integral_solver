# check the error induced by approximating the potential by a point charge near the surface

using BoundaryIntegral
import BoundaryIntegral as BI
using FastGaussQuadrature
using LinearAlgebra
using CSV, DataFrames

function laplace3d_doublelayer(x, y, z, ptx, pty, ptz)
    r2 = (x - ptx)^2 + (y - pty)^2 + (z - ptz)^2
    return (z - ptz) / (r2 * sqrt(r2))
end

CSV.write("data/doublelayer_represent_error.csv", DataFrame(pt_z = Float64[], n = Int[], err = Float64[]))

pt_x = 0.2
pt_y = 0.3

trg_x = range(-1.0, 1.0, length = 100)
trg_y = range(-1.0, 1.0, length = 100)

for pt_z in [0.1, 0.2, 0.5, 1.0]

    ref_val = zeros(length(trg_x), length(trg_y))
    for i in 1:length(trg_x), j in 1:length(trg_y)
        ref_val[i, j] = laplace3d_doublelayer(trg_x[i], trg_y[j], 0.0, pt_x, pt_y, pt_z)
    end
    ref_val = reshape(ref_val, length(trg_x) * length(trg_y))

    for n in 2:2:64
        x, w = gausslegendre(n)
        val = zeros(n, n)
        for i in 1:n, j in 1:n
            val[i, j] = laplace3d_doublelayer(x[i], x[j], 0.0, pt_x, pt_y, pt_z)
        end
        vec_val = reshape(val, n^2)
        U = BI.interp_matrix_2d_gl_tensor(x, w, x, w, trg_x, trg_y)

        upsampled_val = U * vec_val

        err = norm((upsampled_val .- ref_val), Inf)
        CSV.write("data/doublelayer_represent_error.csv", append=true, DataFrame(pt_z = pt_z, n = n, err = err))
    end
end
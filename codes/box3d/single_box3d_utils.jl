using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra, Krylov
using CSV, DataFrames, JLD2

function solve_single_box3d(eps, n_box, n_quad, reduce_quad, n_edge, n_corner, src; rtol = 1e-6)
    dbox = BI.dielectric_box3d(eps, 1.0, n_box, n_quad, reduce_quad, n_edge, n_corner)
    lhs = BI.Lhs_dielectric_mbox3d_fmm3d(dbox, rtol)
    rhs =  BI.Rhs_dielectric_mbox3d(dbox, src, eps)
    sigma, status = Krylov.gmres(lhs, rhs, rtol=rtol, verbose = 1)
    gi = dot(sigma, BI.all_weights(dbox)) + 1 / eps

    n_val = length(sigma)
    n_iter = status.niter

    println("gi = $gi, n_val = $n_val, n_iter = $n_iter")
    return dbox, sigma, gi, n_val, n_iter
end
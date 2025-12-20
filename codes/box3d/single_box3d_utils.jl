using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra, Krylov
using CSV, DataFrames, JLD2

function solve_single_box3d(eps, n_box, n_quad_max, n_quad_min, n_edge, n_corner, src; rtol = 1e-6)
    dbox = BI.dielectric_box3d(eps, 1.0, n_box, n_quad_max, n_quad_min, n_edge, n_corner)
    lhs = BI.Lhs_dielectric_mbox3d_fmm3d(dbox, rtol)
    rhs =  BI.Rhs_dielectric_mbox3d(dbox, src, eps)
    sigma, status = Krylov.gmres(lhs, rhs, rtol=rtol, verbose = 1)
    gi = dot(sigma, BI.all_weights(dbox)) + 1 / eps

    n_val = length(sigma)
    n_iter = status.niter

    println("gi = $gi, n_val = $n_val, n_iter = $n_iter")
    return dbox, sigma, gi, n_val, n_iter
end

function solve_double_box3d(eps_1, eps_2, n_quad_max, n_quad_min, n_edge, n_corner, src; rtol = 1e-6)
    dbox = BI.dielectric_double_box3d(eps_1, eps_2, 1.0, n_quad_max, n_quad_min, n_edge, n_corner)
    lhs = BI.Lhs_dielectric_mbox3d_fmm3d(dbox, rtol)
    rhs =  BI.Rhs_dielectric_mbox3d(dbox, src, eps_2)
    sigma, status = Krylov.gmres(lhs, rhs, rtol=rtol, verbose = 1)
    gi = dot(sigma, BI.all_weights(dbox)) + 1 / eps_2

    n_val = length(sigma)
    n_iter = status.niter

    println("gi = $gi, n_val = $n_val, n_iter = $n_iter")
    return dbox, sigma, gi, n_val, n_iter
end

function solve_single_thin_box3d(eps, Lx, Ly, Lz, nx, ny, nz, n_quad_max, n_quad_min, n_edge, n_corner, src; rtol = 1e-6)
    tbox = BI.dielectric_arbitrary_box3d(eps, 1.0, Lx, Ly, Lz, nx, ny, nz, n_quad_max, n_quad_min, n_edge, n_corner)
    lhs = BI.Lhs_dielectric_mbox3d_fmm3d(tbox, rtol)
    rhs =  BI.Rhs_dielectric_mbox3d(tbox, src, eps)
    sigma, status = Krylov.gmres(lhs, rhs, rtol=rtol, verbose = 1)
    gi = dot(sigma, BI.all_weights(tbox)) + 1 / eps

    n_val = length(sigma)
    n_iter = status.niter

    println("gi = $gi, n_val = $n_val, n_iter = $n_iter")
    return tbox, sigma, gi, n_val, n_iter
end

function solve_single_thin_box3d_far_src(eps, Lx, Ly, Lz, nx, ny, nz, n_quad_max, n_quad_min, n_edge, n_corner, src; rtol = 1e-6)
    tbox = BI.dielectric_arbitrary_box3d(eps, 1.0, Lx, Ly, Lz, nx, ny, nz, n_quad_max, n_quad_min, n_edge, n_corner)
    lhs = BI.Lhs_dielectric_mbox3d_fmm3d(tbox, rtol)
    rhs =  BI.Rhs_dielectric_mbox3d(tbox, src, 1.0)
    sigma, status = Krylov.gmres(lhs, rhs, rtol=rtol, verbose = 1)
    gi = dot(sigma, BI.all_weights(tbox))

    n_val = length(sigma)
    n_iter = status.niter

    println("gi = $gi, n_val = $n_val, n_iter = $n_iter")
    return tbox, sigma, gi, n_val, n_iter
end

function solve_single_box3d_weighted(eps, n_box, n_quad_max, n_quad_min, n_edge, n_corner, src; rtol = 1e-6)
    dbox = BI.dielectric_box3d(eps, 1.0, n_box, n_quad_max, n_quad_min, n_edge, n_corner)
    lhs = BI.Lhs_dielectric_mbox3d_fmm3d_weighted(dbox, rtol)
    rhs =  BI.Rhs_dielectric_mbox3d_weighted(dbox, src, eps)
    sigma, status = Krylov.gmres(lhs, rhs, rtol=rtol, verbose = 1)
    gi = dot(sigma, BI.all_weights(dbox)) + 1 / eps

    n_val = length(sigma)
    n_iter = status.niter

    println("gi = $gi, n_val = $n_val, n_iter = $n_iter")
    return dbox, sigma, gi, n_val, n_iter
end
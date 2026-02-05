# this script is used to test the convergence of DT_corrected against the original one
using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra, Krylov
using CSV, DataFrames, JLD2
using Random


function solve_single_box3d_adaptive_mesh(Lx, Ly, Lz, n_quad, l_ec, eps_in, eps_out, fmm_tol, rhs_tol, ps)
    tbox = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, n_quad, ps, 1.0, l_ec, rhs_tol, eps_in, eps_out, max_depth = 1000)

    lhs = BI.Lhs_dielectric_box3d_fmm3d(tbox, fmm_tol)

    rhs = BI.Rhs_dielectric_box3d(tbox, ps, eps_in)
    sigma, status = Krylov.gmres(lhs, rhs, atol=fmm_tol, verbose = 1)
    total_flux = dot(sigma, BI.all_weights(tbox))

    n_val = length(sigma)
    n_iter = status.niter

    println("total_flux = $total_flux, n_val = $n_val, n_iter = $n_iter")
    return tbox, sigma, total_flux, n_val, n_iter
end

function solve_single_box3d_nonadaptive_mesh(Lx, Ly, Lz, n_quad, l_ec, eps_in, eps_out, fmm_tol, ps)
    tbox = BI.single_dielectric_box3d(Lx, Ly, Lz, n_quad, l_ec, eps_in, eps_out)

    lhs = BI.Lhs_dielectric_box3d_fmm3d(tbox, fmm_tol)

    rhs = BI.Rhs_dielectric_box3d(tbox, ps, eps_in)
    sigma, status = Krylov.gmres(lhs, rhs, atol=fmm_tol, verbose = 1)
    total_flux = dot(sigma, BI.all_weights(tbox))

    n_val = length(sigma)
    n_iter = status.niter

    println("total_flux = $total_flux, n_val = $n_val, n_iter = $n_iter")
    return tbox, sigma, total_flux, n_val, n_iter
end

function solve_single_box3d_adaptive_mesh_corrected(Lx, Ly, Lz, n_quad, l_ec, eps_in, eps_out, fmm_tol, rhs_tol, up_tol, max_order, ps; include_edges_src = false, include_edges_trg = false)
    tbox = BI.single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, n_quad, ps, 1.0, l_ec, rhs_tol, eps_in, eps_out, max_depth = 1000)

    lhs = BI.Lhs_dielectric_box3d_fmm3d_corrected(tbox, fmm_tol, up_tol, max_order, include_edges_src = include_edges_src, include_edges_trg = include_edges_trg)
    rhs = BI.Rhs_dielectric_box3d(tbox, ps, eps_in)

    sigma, status = Krylov.gmres(lhs, rhs, atol=fmm_tol, verbose = 1)
    total_flux = dot(sigma, BI.all_weights(tbox))

    n_val = length(sigma)
    n_iter = status.niter

    println("total_flux = $total_flux, n_val = $n_val, n_iter = $n_iter")
    return tbox, sigma, total_flux, n_val, n_iter
end


function solve_single_box3d_adaptive_mesh_varquad_corrected(Lx, Ly, Lz, n_quad, l_ec, eps_in, eps_out, fmm_tol, rhs_tol, up_tol, max_order, ps; include_edges_src = false, include_edges_trg = false)
    tbox = BI.single_dielectric_box3d_rhs_adaptive_varquad(Lx, Ly, Lz, n_quad, ps, 1.0, l_ec, rhs_tol, eps_in, eps_out, max_depth = 1000)

    lhs = BI.Lhs_dielectric_box3d_fmm3d_corrected(tbox, fmm_tol, up_tol, max_order, include_edges_src = include_edges_src, include_edges_trg = include_edges_trg)
    rhs = BI.Rhs_dielectric_box3d(tbox, ps, eps_in)

    sigma, status = Krylov.gmres(lhs, rhs, atol=fmm_tol, verbose = 1)
    total_flux = dot(sigma, BI.all_weights(tbox))

    n_val = length(sigma)
    n_iter = status.niter

    println("total_flux = $total_flux, n_val = $n_val, n_iter = $n_iter")
    return tbox, sigma, total_flux, n_val, n_iter
end
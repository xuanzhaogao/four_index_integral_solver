using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra, Krylov
using CSV, DataFrames, JLD2
using Random

function solve_single_thin_box3d(Lx, Ly, Lz, n_quad, l_ec, eps_in, eps_out, fmm_tol)
    tbox = BI.single_dielectric_box3d(Lx, Ly, Lz, n_quad, l_ec, eps_in, eps_out)

    lhs = BI.Lhs_dielectric_box3d_fmm3d(tbox, fmm_tol)

    rhs = 3e4 .* BI.Rhs_dielectric_box3d(tbox, PointSource((101.0, 102.0, 103.0), 1.0), eps_in)
    sigma, status = Krylov.gmres(lhs, rhs, rtol=fmm_tol, verbose = 1)
    total_flux = dot(sigma, BI.all_weights(tbox))

    n_val = length(sigma)
    n_iter = status.niter

    println("total_flux = $total_flux, n_val = $n_val, n_iter = $n_iter")
    return tbox, sigma, total_flux, n_val, n_iter
end

function solve_single_thin_box3d_corrected(Lx, Ly, Lz, n_quad, l_ec, eps_in, eps_out, fmm_tol, up_tol, max_order, include_edges_src, include_edges_trg)
    tbox = BI.single_dielectric_box3d(Lx, Ly, Lz, n_quad, l_ec, eps_in, eps_out)

    lhs = BI.Lhs_dielectric_box3d_fmm3d_corrected(tbox, fmm_tol, up_tol, max_order, include_edges_src = include_edges_src, include_edges_trg = include_edges_trg)
    rhs = 3e4 .* BI.Rhs_dielectric_box3d(tbox, PointSource((101.0, 102.0, 103.0), 1.0), eps_out)

    sigma, status = Krylov.gmres(lhs, rhs, rtol=fmm_tol, verbose = 1)
    total_flux = dot(sigma, BI.all_weights(tbox))

    n_val = length(sigma)
    n_iter = status.niter

    println("total_flux = $total_flux, n_val = $n_val, n_iter = $n_iter")
    return tbox, sigma, total_flux, n_val, n_iter
end

p = 4
Lx = 1.0
Ly = 1.0
Lz = 1.0
l_ec = 0.049
eps_in = 4.0
eps_out = 1.0
fmm_tol = 1e-4
up_tol = 1e-5
max_order = 128

target = (0.4, 0.4, 0.4)

tbox_u, sigma_u, total_flux_u, n_val_u, n_iter_u = solve_single_thin_box3d(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol)
potential_u = BI.laplace3d_pottrg_near(tbox_u, target, sigma_u, up_tol)

@show potential_u

# corrected case, edge excluded
tbox_cff, sigma_cff, total_flux_cff, n_val_cff, n_iter_cff = solve_single_thin_box3d_corrected(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, up_tol, max_order, false, false)
potential_cff = BI.laplace3d_pottrg_near(tbox_cff, target, sigma_cff, 1e-6)

@show potential_cff

tbox_ctf, sigma_ctf, total_flux_ctf, n_val_ctf, n_iter_ctf = solve_single_thin_box3d_corrected(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, up_tol, max_order, true, false)
potential_ctf = BI.laplace3d_pottrg_near(tbox_ctf, target, sigma_ctf, 1e-6)

@show potential_ctf

tbox_cft, sigma_cft, total_flux_cft, n_val_cft, n_iter_cft = solve_single_thin_box3d_corrected(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, up_tol, max_order, false, true)
potential_cft = BI.laplace3d_pottrg_near(tbox_cft, target, sigma_cft, 1e-6)

@show potential_cft

# corrected case, edge included
tbox_ctt, sigma_ctt, total_flux_ctt, n_val_ctt, n_iter_ctt = solve_single_thin_box3d_corrected(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, up_tol, max_order, true, true)
potential_ctt = BI.laplace3d_pottrg_near(tbox_ctt, target, sigma_ctt, 1e-6)

@show potential_u, potential_cff, potential_ctf, potential_cft, potential_ctt
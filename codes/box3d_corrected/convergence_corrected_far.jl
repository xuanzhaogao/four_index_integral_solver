# this script is used to test the convergence of DT_corrected against the original one
using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra, Krylov
using CSV, DataFrames, JLD2
using Random

function solve_single_thin_box3d(Lx, Ly, Lz, n_quad, l_ec, eps_in, eps_out, fmm_tol)
    tbox = BI.single_dielectric_box3d(Lx, Ly, Lz, n_quad, l_ec, eps_in, eps_out)

    lhs = BI.Lhs_dielectric_box3d_fmm3d(tbox, fmm_tol)

    rhs =  100.0 .* BI.Rhs_dielectric_box3d(tbox, PointSource((21.0, 22.0, 11.0), 1.0), eps_in)
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
    rhs =  100.0 .* BI.Rhs_dielectric_box3d(tbox, PointSource((21.0, 22.0, 11.0), 1.0), eps_out)

    sigma, status = Krylov.gmres(lhs, rhs, rtol=fmm_tol, verbose = 1)
    total_flux = dot(sigma, BI.all_weights(tbox))

    n_val = length(sigma)
    n_iter = status.niter

    println("total_flux = $total_flux, n_val = $n_val, n_iter = $n_iter")
    return tbox, sigma, total_flux, n_val, n_iter
end

df = joinpath(@__DIR__, "data/convergence_corrected_far.csv")
CSV.write(df, DataFrame(p = Int[], r = Int[], L = Float64[], l_ec = Float64[], total_flux_u = Float64[], total_flux_cff = Float64[], total_flux_ctf = Float64[], total_flux_cft = Float64[], total_flux_ctt = Float64[], n_val_u = Int[], n_val_cff = Int[], n_val_ctf = Int[], n_val_cft = Int[], n_val_ctt = Int[], n_iter_u = Int[], n_iter_cff = Int[], n_iter_ctf = Int[], n_iter_cft = Int[], n_iter_ctt = Int[]))

begin
    for L in [5.0, 10.0, 20.0]
        Lx = L
        Ly = L
        Lz = 1.0

        ps = [2, 4, 6]
        rs = 0:2:6
        eps_in = 4.0
        eps_out = 1.0

        fmm_tol = 1e-4
        up_tol = 1e-5
        max_order = 48

        for p in ps
            for r in rs
                l_ec = 1.0 / 2^r * 1.01

                # uncorrected case
                tbox_u, sigma_u, total_flux_u, n_val_u, n_iter_u = solve_single_thin_box3d(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol)

                # corrected case, edge excluded
                tbox_cff, sigma_cff, total_flux_cff, n_val_cff, n_iter_cff = solve_single_thin_box3d_corrected(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, up_tol, max_order, false, false)

                tbox_ctf, sigma_ctf, total_flux_ctf, n_val_ctf, n_iter_ctf = solve_single_thin_box3d_corrected(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, up_tol, max_order, true, false)

                tbox_cft, sigma_cft, total_flux_cft, n_val_cft, n_iter_cft = solve_single_thin_box3d_corrected(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, up_tol, max_order, false, true)

                # corrected case, edge included
                tbox_ctt, sigma_ctt, total_flux_ctt, n_val_ctt, n_iter_ctt = solve_single_thin_box3d_corrected(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, up_tol, max_order, true, true)

                @show p, r, L, l_ec, total_flux_u, total_flux_cff, total_flux_ctf, total_flux_cft, total_flux_ctt

                res = Dict("tbox_u" => tbox_u, "sigma_u" => sigma_u, "tbox_cff" => tbox_cff, "sigma_cff" => sigma_cff, "tbox_ctf" => tbox_ctf, "sigma_ctf" => sigma_ctf, "tbox_cft" => tbox_cft, "sigma_cft" => sigma_cft, "tbox_ctt" => tbox_ctt, "sigma_ctt" => sigma_ctt)

                save(joinpath(@__DIR__, "cache/convergence_corrected_far_L$(L)_p$(p)_r$(r).jld2"), "res", res)

                CSV.write(df, DataFrame(p = p, r = r, L = L, l_ec = l_ec, total_flux_u = total_flux_u, total_flux_cff = total_flux_cff, total_flux_ctf = total_flux_ctf, total_flux_cft = total_flux_cft, total_flux_ctt = total_flux_ctt, n_val_u = n_val_u, n_val_cff = n_val_cff, n_val_ctf = n_val_ctf, n_val_cft = n_val_cft, n_val_ctt = n_val_ctt, n_iter_u = n_iter_u, n_iter_cff = n_iter_cff, n_iter_ctf = n_iter_ctf, n_iter_cft = n_iter_cft, n_iter_ctt = n_iter_ctt), append = true)
            end
        end
    end
end
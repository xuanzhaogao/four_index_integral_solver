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

function solve_single_thin_box3d_corrected(Lx, Ly, Lz, n_quad, l_ec, eps_in, eps_out, fmm_tol, up_tol, max_order, include_edges)
    tbox = BI.single_dielectric_box3d(Lx, Ly, Lz, n_quad, l_ec, eps_in, eps_out)

    lhs = BI.Lhs_dielectric_box3d_fmm3d_corrected(tbox, fmm_tol, up_tol, max_order, include_edges = include_edges)
    rhs =  100.0 .* BI.Rhs_dielectric_box3d(tbox, PointSource((21.0, 22.0, 11.0), 1.0), eps_out)

    sigma, status = Krylov.gmres(lhs, rhs, rtol=fmm_tol, verbose = 1)
    total_flux = dot(sigma, BI.all_weights(tbox))

    n_val = length(sigma)
    n_iter = status.niter

    println("total_flux = $total_flux, n_val = $n_val, n_iter = $n_iter")
    return tbox, sigma, total_flux, n_val, n_iter
end

df = joinpath(@__DIR__, "data/convergence_corrected_far.csv")
# CSV.write(df, DataFrame(p = Int[], r = Int[], L = Float64[], l_ec = Float64[], total_flux_u = Float64[], total_flux_cf = Float64[], total_flux_ct = Float64[], n_val_u = Int[], n_val_cf = Int[], n_val_ct = Int[], n_iter_u = Int[], n_iter_cf = Int[], n_iter_ct = Int[]))

begin
    # for L in [5.0, 10.0, 20.0]
    for L in [20.0]
        Lx = L
        Ly = L
        Lz = 1.0

        l_panel = 1.0
        ps = [2, 4, 6]
        rs = 0:2:8
        eps_in = 4.0
        eps_out = 1.0

        fmm_tol = 1e-4
        up_tol = 1e-5
        max_order = 12

        for p in ps
            for r in rs
                l_ec = l_panel / 2^r * 1.01

                # uncorrected case
                tbox_u, sigma_u, total_flux_u, n_val_u, n_iter_u = solve_single_thin_box3d(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol)

                # corrected case, edge excluded
                tbox_cf, sigma_cf, total_flux_cf, n_val_cf, n_iter_cf = solve_single_thin_box3d_corrected(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, up_tol, max_order, false)

                # corrected case, edge included
                tbox_ct, sigma_ct, total_flux_ct, n_val_ct, n_iter_ct = solve_single_thin_box3d_corrected(Lx, Ly, Lz, p, l_ec, eps_in, eps_out, fmm_tol, up_tol, max_order, true)

                @show p, r, L, l_ec, total_flux_u, total_flux_cf, total_flux_ct

                res = Dict("tbox_u" => tbox_u, "sigma_u" => sigma_u, "tbox_cf" => tbox_cf, "sigma_cf" => sigma_cf, "tbox_ct" => tbox_ct, "sigma_ct" => sigma_ct)

                save(joinpath(@__DIR__, "cache/convergence_corrected_far_L$(L)_p$(p)_r$(r).jld2"), "res", res)

                CSV.write(df, DataFrame(p = p, r = r, L = L, l_ec = l_ec, total_flux_u = total_flux_u, total_flux_cf = total_flux_cf, total_flux_ct = total_flux_ct, n_val_u = n_val_u, n_val_cf = n_val_cf, n_val_ct = n_val_ct, n_iter_u = n_iter_u, n_iter_cf = n_iter_cf, n_iter_ct = n_iter_ct), append = true)
            end
        end
    end
end
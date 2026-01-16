# this script is used to test the convergence of DT_corrected against the original one
using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra, Krylov
using CSV, DataFrames, JLD2
using Random

function solve_single_thin_box3d(Lx, Ly, Lz, n_quad, l_panel, l_ec, eps_in, eps_out, fmm_tol)
    tbox = BI.single_dielectric_box3d(Lx, Ly, Lz, n_quad, l_panel, l_ec, eps_in, eps_out)

    lhs = BI.Lhs_dielectric_box3d_fmm3d(tbox, fmm_tol)

    rhs =  BI.Rhs_dielectric_box3d(tbox, PointSource((0.2, 0.3, 0.1), 1.0), eps_in)
    sigma, status = Krylov.gmres(lhs, rhs, rtol=fmm_tol, verbose = 1)
    total_flux = dot(sigma, BI.all_weights(tbox)) + 1 / eps_in

    n_val = length(sigma)
    n_iter = status.niter

    println("total_flux = $total_flux, n_val = $n_val, n_iter = $n_iter")
    return tbox, sigma, total_flux, n_val, n_iter
end

function solve_single_thin_box3d_corrected(Lx, Ly, Lz, n_quad, l_panel, l_ec, eps_in, eps_out, fmm_tol, up_tol, max_order, l_min)
    tbox = BI.single_dielectric_box3d(Lx, Ly, Lz, n_quad, l_panel, l_ec, eps_in, eps_out)

    lhs = BI.Lhs_dielectric_box3d_fmm3d_corrected(tbox, fmm_tol, up_tol, max_order, l_min)
    rhs =  BI.Rhs_dielectric_box3d(tbox, PointSource((0.2, 0.3, 0.1), 1.0), eps_in)

    sigma, status = Krylov.gmres(lhs, rhs, rtol=fmm_tol, verbose = 1)
    total_flux = dot(sigma, BI.all_weights(tbox)) + 1 / eps_in

    n_val = length(sigma)
    n_iter = status.niter

    println("total_flux = $total_flux, n_val = $n_val, n_iter = $n_iter")
    return tbox, sigma, total_flux, n_val, n_iter
end

df = joinpath(@__DIR__, "data/convergence_corrected.csv")
CSV.write(df, DataFrame(p = Int[], r = Int[], L = Float64[], l_ec = Float64[], total_flux_u = Float64[], total_flux_c = Float64[], n_val_u = Int[], n_val_c = Int[], n_iter_u = Int[], n_iter_c = Int[]))

begin
    for L in [5.0, 10.0, 20.0]
        Lx = L
        Ly = L
        Lz = 1.0

        l_panel = 1.0
        # ps = [2, 3, 4]
        ps = [6]
        rs = 0:2:6
        eps_in = 4.0
        eps_out = 1.0

        fmm_tol = 1e-4
        up_tol = 1e-6
        max_order = 12

        for p in ps
            for r in rs
                l_ec = l_panel / 2^r * 1.01

                # uncorrected case
                tbox_u, sigma_u, total_flux_u, n_val_u, n_iter_u = solve_single_thin_box3d(Lx, Ly, Lz, p, l_panel, l_ec, eps_in, eps_out, fmm_tol)

                # corrected case
                tbox_c, sigma_c, total_flux_c, n_val_c, n_iter_c = solve_single_thin_box3d_corrected(Lx, Ly, Lz, p, l_panel, l_ec, eps_in, eps_out, fmm_tol, up_tol, max_order, l_ec / 2)

                @show p, r, L, l_ec, total_flux_u, total_flux_c

                res = Dict("tbox_u" => tbox_u, "sigma_u" => sigma_u, "tbox_c" => tbox_c, "sigma_c" => sigma_c)

                save(joinpath(@__DIR__, "cache/convergence_corrected_L$(L)_p$(p)_r$(r).jld2"), "res", res)

                CSV.write(df, DataFrame(p = p, r = r, L = L, l_ec = l_ec, total_flux_u = total_flux_u, total_flux_c = total_flux_c, n_val_u = n_val_u, n_val_c = n_val_c, n_iter_u = n_iter_u, n_iter_c = n_iter_c), append = true)
            end
        end
    end
end
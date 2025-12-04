# test the convergence of a single box in 3d space
using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra, Krylov
using CSV, DataFrames, JLD2

function solve_single_box3d(eps, n_box, n_quad, n_edge, n_corner, src)
    dbox = BI.dielectric_box3d(eps, 1.0, n_box, n_quad, n_edge, n_corner)
    lhs = BI.Lhs_dielectric_mbox3d_fmm3d(dbox, 1e-6)
    rhs =  BI.Rhs_dielectric_mbox3d(dbox, src, eps)
    sigma, _ = Krylov.gmres(lhs, rhs, rtol=1e-6, verbose = 1)
    gi = dot(sigma, BI.all_weights(dbox)) + 1 / eps
    println("gi = $gi")
    return dbox, sigma, gi
end

trgs = load(joinpath(@__DIR__, "data/sphere_trgs.jld"))["trgs"]

# res_ref = solve_single_box3d(4.0, 6, 12, 8, 10, (0.2, 0.3, 0.4)), gi = 1.0005676883216563
res_ref = load(joinpath(@__DIR__, "data/single_box3d_ref.jld"))["res"]

df = joinpath(@__DIR__, "data/single_box3d.csv")
CSV.write(df, DataFrame(n_boxes = [], n_quad = [], n_edge = [], n_corner = [], gi = [], pot_abserr = [], pot_relerr = []))

# 1. fix n_box = 1, n_edge = 0, n_corner = 0, only vary n_quad
for n_quad in [2, 4, 8, 16, 32, 64]
    dbox, sigma, gi = solve_single_box3d(4.0, 1, n_quad, 0, 0, (0.2, 0.3, 0.4))
    S_trg = BI.laplace3d_pottarg_fmm3d(dbox, trgs, 1e-6)
    pot_trg = S_trg * sigma
    pot_abserr = maximum(abs.(pot_trg - res_ref[4]))
    pot_relerr = maximum(abs.((pot_trg - res_ref[4]) ./ res_ref[4]))
    CSV.write(df, DataFrame(n_boxes = [1], n_quad = [n_quad], n_edge = [0], n_corner = [0], gi = [gi], pot_abserr = [pot_abserr], pot_relerr = [pot_relerr]), append = true)
end

# 2. fix n_quad = 12, vary n_box (1, 2, 3, 4, 5, 6)
# non_adaptive refinement
for n_box in [1, 2, 3, 4, 5, 6]
    dbox, sigma, gi = solve_single_box3d(4.0, n_box, 12, 0, 0, (0.2, 0.3, 0.4))
    S_trg = BI.laplace3d_pottarg_fmm3d(dbox, trgs, 1e-6)
    pot_trg = S_trg * sigma
    pot_abserr = maximum(abs.(pot_trg - res_ref[4]))
    pot_relerr = maximum(abs.((pot_trg - res_ref[4]) ./ res_ref[4]))
    CSV.write(df, DataFrame(n_boxes = [n_box], n_quad = [12], n_edge = [0], n_corner = [0], gi = [gi], pot_abserr = [pot_abserr], pot_relerr = [pot_relerr]), append = true)
end

# adaptive refinement for edges
# 3. fix n_quad = 12, n_box = 4, vary_n_edge (2, 4, 6, 8)
for n_edge in [2, 4, 6, 8]
    dbox, sigma, gi = solve_single_box3d(4.0, 4, 12, n_edge, n_edge, (0.2, 0.3, 0.4))
    S_trg = BI.laplace3d_pottarg_fmm3d(dbox, trgs, 1e-6)
    pot_trg = S_trg * sigma
    pot_abserr = maximum(abs.(pot_trg - res_ref[4]))
    pot_relerr = maximum(abs.((pot_trg - res_ref[4]) ./ res_ref[4]))
    CSV.write(df, DataFrame(n_boxes = [4], n_quad = [12], n_edge = [n_edge], n_corner = [0], gi = [gi], pot_abserr = [pot_abserr], pot_relerr = [pot_relerr]), append = true)
end

# adaptive refinement for corners
# 4. fix n_quad = 12, n_box = 4, n_edges = 6, n_corners = [7, 8, 9, 10]
for n_corner in [7, 8, 9, 10]
    dbox, sigma, gi = solve_single_box3d(4.0, 4, 12, 6, n_corner, (0.2, 0.3, 0.4))
    S_trg = BI.laplace3d_pottarg_fmm3d(dbox, trgs, 1e-6)
    pot_trg = S_trg * sigma
    pot_abserr = maximum(abs.(pot_trg - res_ref[4]))
    pot_relerr = maximum(abs.((pot_trg - res_ref[4]) ./ res_ref[4]))
    CSV.write(df, DataFrame(n_boxes = [4], n_quad = [12], n_edge = [6], n_corner = [n_corner], gi = [gi], pot_abserr = [pot_abserr], pot_relerr = [pot_relerr]), append = true)
end
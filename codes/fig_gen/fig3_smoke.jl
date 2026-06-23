#=
Smoke test for the fig3 redo with proper edge correction.

Goal: validate lhs_dielectric_box3d_fmm3d_corrected(...; correct_edges=true)
in the standalone fig3 setting and measure (a) per-solve cost vs the old
uncorrected operator and (b) how the far-field potential differs / converges.

Cheap configs only (shallow k). Prints timing, N, niter, gmres_res, and the
relative far-field difference corrected-vs-uncorrected at each level.
=#

using LinearAlgebra
using FastGaussQuadrature
using BoundaryIntegral
const BI = BoundaryIntegral

const Lx, Ly, Lz = 2.0, 2.0, 2.0
const eps_d, eps_0 = 10.0, 1.0
const eps_src = eps_0
const ps = BI.PointSource((0.5, 0.6, 100.0), 1.0e4)
const fmm_tol   = 1e-6
const gmres_tol = 1e-6
const up_tol    = 1e-6
const max_order = 64

# far-field sphere (same as fig3_data_pscan.jl)
const R_far = 5.0
const n_theta = 12
const n_phi   = 24
far_targets = Matrix{Float64}(undef, 3, n_theta * n_phi)
let kk = 1
    for i in 1:n_theta
        θ = π * (i - 0.5) / n_theta
        for j in 1:n_phi
            φ = 2π * (j - 1) / n_phi
            far_targets[1, kk] = R_far * sin(θ) * cos(φ)
            far_targets[2, kk] = R_far * sin(θ) * sin(φ)
            far_targets[3, kk] = R_far * cos(θ)
            kk += 1
        end
    end
end

function solve_far(p_quad, l_ec; corrected::Bool)
    iface = single_dielectric_box3d(Lx, Ly, Lz, p_quad, l_ec, eps_d, eps_0)
    npan = length(iface.panels)
    N = npan * p_quad^2
    t_lhs = @elapsed Lhs = corrected ?
        lhs_dielectric_box3d_fmm3d_corrected(iface, fmm_tol, up_tol, max_order; correct_edges = true) :
        lhs_dielectric_box3d_fmm3d(iface, fmm_tol)
    rhs = rhs_dielectric_box3d(iface, ps, eps_src)
    t_solve = @elapsed sigma = BI.solve_gmres(Lhs, rhs, gmres_tol, gmres_tol)
    gres = norm(Lhs * sigma - rhs) / norm(rhs)
    S = BI.laplace3d_pottrg(iface, far_targets)
    u_far = S * sigma
    return (; npan, N, t_lhs, t_solve, gres, u_far)
end

p_quad = 3
println("# fig3 smoke: p=$p_quad  (corrected = correct_edges=true vs uncorrected pure FMM)")
println("# k   l_ec        npan      N         t_lhs[s]  t_solve[s]  gmres_res   rel||u_c-u_u||")
u_c_prev = nothing
u_u_prev = nothing
for k in 3:5
    l_ec = 1.01 / 2.0^k
    ru = solve_far(p_quad, l_ec; corrected = false)
    rc = solve_far(p_quad, l_ec; corrected = true)
    reldiff = norm(rc.u_far .- ru.u_far) / norm(rc.u_far)
    println(rpad("  $k", 5), rpad(round(l_ec; sigdigits=3), 12),
            rpad(rc.npan, 10), rpad(rc.N, 10),
            "corr[lhs=$(round(rc.t_lhs; digits=1)) solve=$(round(rc.t_solve; digits=1)) res=$(round(rc.gres; sigdigits=2))] ",
            "unc[lhs=$(round(ru.t_lhs; digits=1)) solve=$(round(ru.t_solve; digits=1)) res=$(round(ru.gres; sigdigits=2))] ",
            "reldiff=$(round(reldiff; sigdigits=4))")
    # self-convergence of each operator across k
    if u_c_prev !== nothing
        dc = norm(rc.u_far .- u_c_prev) / norm(u_c_prev)
        du = norm(ru.u_far .- u_u_prev) / norm(u_u_prev)
        println("       self-conv k-1->k:  corrected Δ=$(round(dc; sigdigits=4))   uncorrected Δ=$(round(du; sigdigits=4))")
    end
    global u_c_prev = rc.u_far
    global u_u_prev = ru.u_far
end

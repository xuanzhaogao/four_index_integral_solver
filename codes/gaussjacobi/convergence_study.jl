import BoundaryIntegral as BI
using LinearAlgebra, Printf

eps_in  = 1.0; eps_out = 200.0
ps      = BI.PointSource((0.001, 0.001), 1.0)
exact_flux = 1.0/eps_out - 1.0/eps_in

println("eps_in=$eps_in, eps_out=$eps_out, gamma=$(round(BI.corner_singularity_power(pi/2, eps_in, eps_out), sigdigits=4))")
println("exact_flux = $exact_flux")
println()

# -------------------------------------------------------
# Table 1: sweep n_quad, fixed l_corner=0.05
# -------------------------------------------------------
println("=== Table 1: n_quad sweep  (l_corner=0.05, l_panel=0.2) ===")
println(@sprintf("%-8s  %-6s  %-12s  %-12s  %-8s", "n_quad", "n_pts", "GL_flux_err", "GJ_flux_err", "ratio"))

for nq in [4, 6, 8, 10, 12, 16, 20]
    box_gl = BI.single_dielectric_box2d(1.0, 1.0, nq, 0.2, 0.05, eps_in, eps_out, Float64; use_singular=false)
    box_gj = BI.single_dielectric_box2d(1.0, 1.0, nq, 0.2, 0.05, eps_in, eps_out, Float64; use_singular=true)

    lhs_gl = BI.lhs_dielectric_box2d(box_gl)
    rhs_gl = BI.rhs_dielectric_box2d(box_gl, ps, eps_in)
    x_gl   = BI.solve_lu(lhs_gl, rhs_gl)

    lhs_gj = BI.lhs_dielectric_box2d(box_gj)
    rhs_gj = BI.rhs_dielectric_box2d(box_gj, ps, eps_in)
    x_gj   = BI.solve_lu(lhs_gj, rhs_gj)

    e_gl = abs(dot(BI.all_weights(box_gl), x_gl) - exact_flux)
    e_gj = abs(dot(BI.all_weights(box_gj), x_gj) - exact_flux)
    npts = BI.num_points(box_gl)
    println(@sprintf("%-8d  %-6d  %-12.3e  %-12.3e  %-8.2f", nq, npts, e_gl, e_gj, e_gl/e_gj))
    flush(stdout)
end

println()

# -------------------------------------------------------
# Table 2: sweep l_corner, fixed n_quad=12
# -------------------------------------------------------
println("=== Table 2: l_corner sweep  (n_quad=12, l_panel=0.2) ===")
println(@sprintf("%-10s  %-6s  %-12s  %-12s  %-8s", "l_corner", "n_pts", "GL_flux_err", "GJ_flux_err", "ratio"))

for lc in [0.2, 0.1, 0.05, 0.02, 0.01]
    box_gl = BI.single_dielectric_box2d(1.0, 1.0, 12, 0.2, lc, eps_in, eps_out, Float64; use_singular=false)
    box_gj = BI.single_dielectric_box2d(1.0, 1.0, 12, 0.2, lc, eps_in, eps_out, Float64; use_singular=true)

    lhs_gl = BI.lhs_dielectric_box2d(box_gl)
    rhs_gl = BI.rhs_dielectric_box2d(box_gl, ps, eps_in)
    x_gl   = BI.solve_lu(lhs_gl, rhs_gl)

    lhs_gj = BI.lhs_dielectric_box2d(box_gj)
    rhs_gj = BI.rhs_dielectric_box2d(box_gj, ps, eps_in)
    x_gj   = BI.solve_lu(lhs_gj, rhs_gj)

    e_gl = abs(dot(BI.all_weights(box_gl), x_gl) - exact_flux)
    e_gj = abs(dot(BI.all_weights(box_gj), x_gj) - exact_flux)
    npts = BI.num_points(box_gl)
    println(@sprintf("%-10.3f  %-6d  %-12.3e  %-12.3e  %-8.2f", lc, npts, e_gl, e_gj, e_gl/e_gj))
    flush(stdout)
end

println()

# -------------------------------------------------------
# Table 3: n_total sweep — equal work comparison
# increase n_quad for GJ, refine l_corner for GL to match n_pts
# -------------------------------------------------------
println("=== Table 3: equal n_pts comparison ===")
println(@sprintf("%-8s  %-6s  %-12s  %-8s  %-6s  %-12s", "GL_nq", "n_pts_GL", "GL_flux_err", "GJ_nq", "n_pts_GJ", "GJ_flux_err"))

gl_configs = [(4,0.05),(6,0.05),(8,0.05),(12,0.05),(16,0.05),(20,0.05)]
gj_configs = [(4,0.05),(6,0.05),(8,0.05),(12,0.05),(16,0.05),(20,0.05)]

for ((nq_gl,lc_gl),(nq_gj,lc_gj)) in zip(gl_configs, gj_configs)
    box_gl = BI.single_dielectric_box2d(1.0,1.0,nq_gl,0.2,lc_gl,eps_in,eps_out,Float64; use_singular=false)
    box_gj = BI.single_dielectric_box2d(1.0,1.0,nq_gj,0.2,lc_gj,eps_in,eps_out,Float64; use_singular=true)

    x_gl = BI.solve_lu(BI.lhs_dielectric_box2d(box_gl), BI.rhs_dielectric_box2d(box_gl,ps,eps_in))
    x_gj = BI.solve_lu(BI.lhs_dielectric_box2d(box_gj), BI.rhs_dielectric_box2d(box_gj,ps,eps_in))

    e_gl = abs(dot(BI.all_weights(box_gl), x_gl) - exact_flux)
    e_gj = abs(dot(BI.all_weights(box_gj), x_gj) - exact_flux)

    println(@sprintf("%-8d  %-6d  %-12.3e  %-8d  %-6d  %-12.3e",
        nq_gl, BI.num_points(box_gl), e_gl, nq_gj, BI.num_points(box_gj), e_gj))
    flush(stdout)
end

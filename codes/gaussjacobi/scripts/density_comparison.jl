"""
Compare the density σ(s) near box corners for three discretizations:
  1. GL-ref  : highly refined GL (reference solution)
  2. GL      : moderate GL
  3. GL+GJ   : GJ on corner panels, GL elsewhere

We extract σ on each edge approaching each corner separately.
The density should behave as σ(s) ~ C · s^α  where α = γ-1, s = dist from corner.
"""

import BoundaryIntegral as BI
using LinearAlgebra, Printf

# ── problem parameters ───────────────────────────────────────────────────────
eps_in  = 1.0
eps_out = 200.0
ps      = BI.PointSource((0.33, 0.44), 1.0)   # asymmetric to avoid cancellation
exact_flux = 1.0/eps_out - 1.0/eps_in

gamma = BI.corner_singularity_power(pi/2, eps_in, eps_out)
alpha = gamma - 1   # density power: σ ~ s^α near corner
@printf("γ = %.8f,  α = γ-1 = %.8f\n", gamma, alpha)
@printf("exact_flux = %.12e\n\n", exact_flux)

# ── solve helper ─────────────────────────────────────────────────────────────
function solve_box(box)
    lhs = BI.lhs_dielectric_box2d(box)
    rhs = BI.rhs_dielectric_box2d(box, ps, eps_in)
    x   = BI.solve_lu(lhs, rhs)
    flux_err = abs(dot(BI.all_weights(box), x) - exact_flux)
    return x, flux_err
end

"""
Extract (s, σ) on each panel individually. Returns a vector of
(panel_corners, is_singular, s_vals, sigma_vals) tuples.
"""
function extract_panel_data(box, x)
    result = []
    offset = 0
    for p in box.panels
        np = BI.num_points(p)
        x_panel = x[offset+1:offset+np]
        sigma_phys = copy(x_panel)
        if p.is_singular
            for j in 1:np
                sigma_phys[j] = x_panel[j] * (1 + p.gl_xs[j])^p.singular_exponent
            end
        end
        push!(result, (corners=p.corners, is_singular=p.is_singular,
                       points=copy(p.points), sigma=sigma_phys,
                       sigma_bare=copy(x_panel), gl_xs=copy(p.gl_xs)))
        offset += np
    end
    return result
end

"""
For a given corner, find the two adjacent panels (one on each edge)
and return (s, σ) data sorted by s = distance from corner.
"""
function corner_edge_data(panel_data, corner_pt; tol=1e-8)
    edges = []
    for pd in panel_data
        d_a = norm(pd.corners[1] .- corner_pt)
        d_b = norm(pd.corners[2] .- corner_pt)
        if d_a < tol || d_b < tol
            s_vals = [norm(pt .- corner_pt) for pt in pd.points]
            perm = sortperm(s_vals)
            push!(edges, (corners=pd.corners, is_singular=pd.is_singular,
                          s=s_vals[perm], sigma=pd.sigma[perm],
                          sigma_bare=pd.sigma_bare[perm]))
        end
    end
    return edges
end

# ── 1. Reference: very refined GL ────────────────────────────────────────────
@printf("Building reference GL (n=32, l_corner=0.001)...\n")
box_ref = BI.single_dielectric_box2d(1.0, 1.0, 32, 0.05, 0.001,
                                      eps_in, eps_out; use_singular=false)
x_ref, err_ref = solve_box(box_ref)
@printf("  N=%d, flux_err=%.3e\n", BI.num_points(box_ref), err_ref)
pd_ref = extract_panel_data(box_ref, x_ref)

# ── 2. Moderate GL ───────────────────────────────────────────────────────────
n_quad = 16; l_panel = 0.2; l_corner = 0.05
@printf("Building GL (n=%d, lc=%.3f)...\n", n_quad, l_corner)
box_gl = BI.single_dielectric_box2d(1.0, 1.0, n_quad, l_panel, l_corner,
                                     eps_in, eps_out; use_singular=false)
x_gl, err_gl = solve_box(box_gl)
@printf("  N=%d, flux_err=%.3e\n", BI.num_points(box_gl), err_gl)
pd_gl = extract_panel_data(box_gl, x_gl)

# ── 3. GL+GJ ────────────────────────────────────────────────────────────────
@printf("Building GL+GJ (n=%d, lc=%.3f)...\n", n_quad, l_corner)
box_gj = BI.single_dielectric_box2d(1.0, 1.0, n_quad, l_panel, l_corner,
                                     eps_in, eps_out; use_singular=true)
x_gj, err_gj = solve_box(box_gj)
@printf("  N=%d, flux_err=%.3e\n", BI.num_points(box_gj), err_gj)
pd_gj = extract_panel_data(box_gj, x_gj)

# ── power-law fit ────────────────────────────────────────────────────────────
function fit_power_law(s, sig; s_max=Inf)
    mask = (s .> 0) .& (s .< s_max) .& (abs.(sig) .> 0)
    sum(mask) < 3 && return NaN, NaN, sum(mask)
    ls = log.(s[mask])
    lsig = log.(abs.(sig[mask]))
    n = length(ls)
    b = (n * dot(ls, lsig) - sum(ls)*sum(lsig)) / (n*dot(ls,ls) - sum(ls)^2)
    a = (sum(lsig) - b*sum(ls)) / n
    return b, exp(a), n
end

# ── print per-corner, per-edge density comparison ────────────────────────────
corners = [(-0.5, 0.5), (0.5, 0.5), (0.5, -0.5), (-0.5, -0.5)]

for corner in corners
    println("\n" * "="^90)
    @printf("Corner %s\n", corner)
    println("="^90)

    edges_ref = corner_edge_data(pd_ref, corner)
    edges_gl  = corner_edge_data(pd_gl, corner)
    edges_gj  = corner_edge_data(pd_gj, corner)

    # group edges by direction (horizontal or vertical)
    for (ie, (er, egl, egj)) in enumerate(zip(edges_ref, edges_gl, edges_gj))
        @printf("\n  Edge %d: ref corners=%s, GL corners=%s, GJ corners=%s\n",
                ie, er.corners, egl.corners, egj.corners)
        @printf("          GJ is_singular=%s\n", egj.is_singular)

        # power-law fit for this edge
        a_ref, C_ref, _ = fit_power_law(er.s, er.sigma)
        a_gl,  C_gl,  _ = fit_power_law(egl.s, egl.sigma)
        a_gj,  C_gj,  _ = fit_power_law(egj.s, egj.sigma)

        @printf("  Power-law fit (all pts):  ref α=%.4f  GL α=%.4f  GJ α=%.4f  (expected %.4f)\n",
                a_ref, a_gl, a_gj, alpha)

        # print first few data points for each
        println("\n  --- Reference (first 10) ---")
        @printf("  %-14s  %-16s  %-14s\n", "s", "σ(s)", "|σ|/s^α")
        for i in 1:min(10, length(er.s))
            @printf("  %.8e  %+.10e  %.8e\n",
                    er.s[i], er.sigma[i], abs(er.sigma[i]) / er.s[i]^alpha)
        end

        println("\n  --- GL ---")
        @printf("  %-14s  %-16s  %-14s\n", "s", "σ(s)", "|σ|/s^α")
        for i in 1:min(10, length(egl.s))
            @printf("  %.8e  %+.10e  %.8e\n",
                    egl.s[i], egl.sigma[i], abs(egl.sigma[i]) / egl.s[i]^alpha)
        end

        println("\n  --- GL+GJ ---")
        @printf("  %-14s  %-16s  %-14s  %-14s\n", "s", "σ(s)", "|σ|/s^α", "σ_bare")
        for i in 1:min(10, length(egj.s))
            @printf("  %.8e  %+.10e  %.8e  %+.10e\n",
                    egj.s[i], egj.sigma[i], abs(egj.sigma[i]) / egj.s[i]^alpha,
                    egj.sigma_bare[i])
        end
    end
end

# ── summary: per-edge power-law fits ─────────────────────────────────────────
println("\n" * "="^90)
@printf("Summary of power-law fits: σ ~ C·s^α  (expected α = %.6f)\n", alpha)
println("="^90)
@printf("%-20s  %-8s  %-10s  %-10s  %-10s\n", "corner/edge", "sing?", "α_ref", "α_GL", "α_GJ")
println("-"^65)

for corner in corners
    edges_ref = corner_edge_data(pd_ref, corner)
    edges_gl  = corner_edge_data(pd_gl, corner)
    edges_gj  = corner_edge_data(pd_gj, corner)
    for (ie, (er, egl, egj)) in enumerate(zip(edges_ref, edges_gl, edges_gj))
        a_ref, _, _ = fit_power_law(er.s, er.sigma)
        a_gl,  _, _ = fit_power_law(egl.s, egl.sigma)
        a_gj,  _, _ = fit_power_law(egj.s, egj.sigma)
        label = @sprintf("%s e%d", corner, ie)
        @printf("%-20s  %-8s  %-10.4f  %-10.4f  %-10.4f\n",
                label, egj.is_singular, a_ref, a_gl, a_gj)
    end
end

# ── write data for the first corner, first edge (for plotting) ───────────────
dir = joinpath(@__DIR__, "..", "data")
corner = corners[1]
edges_ref = corner_edge_data(pd_ref, corner)
edges_gl  = corner_edge_data(pd_gl, corner)
edges_gj  = corner_edge_data(pd_gj, corner)

for (label, edges) in [("ref", edges_ref), ("gl", edges_gl), ("gj", edges_gj)]
    open(joinpath(dir, "density_$(label).dat"), "w") do io
        @printf(io, "# corner=%s\n", corner)
        for (ie, e) in enumerate(edges)
            @printf(io, "# edge %d: corners=%s  is_singular=%s\n", ie, e.corners, e.is_singular)
            for (s, sig) in zip(e.s, e.sigma)
                @printf(io, "%.15e  %.15e  %d\n", s, sig, ie)
            end
        end
    end
end
# also write sigma_bare for GJ
open(joinpath(dir, "density_gj_bare.dat"), "w") do io
    @printf(io, "# corner=%s\n", corner)
    for (ie, e) in enumerate(edges_gj)
        e.is_singular || continue
        for (s, sb) in zip(e.s, e.sigma_bare)
            @printf(io, "%.15e  %.15e  %d\n", s, sb, ie)
        end
    end
end

println("\nDone. Data files written to $dir")

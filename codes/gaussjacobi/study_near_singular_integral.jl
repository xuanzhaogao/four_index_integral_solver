"""
Study: accurate evaluation of ∫ K(s(t), xₖ) · σ(t) dt
where σ(t) = (1+t)^β φ(t) and φ is given by nodal values at GJ nodes.

Three methods compared:
  (A) Direct GJ sum:   Σⱼ K(s(tⱼ), xₖ) · wⱼ · φⱼ
  (B) hcub + barycentric: ∫ K · (1+t)^β · [Σⱼ φⱼ Lⱼ(t)] dt   (barycentric_row! inside integrand)
  (C) hcub + Jacobi:   same integral but φ(t) evaluated via Jacobi expansion
  (REF) hcub with very tight tolerance on method (B) as the reference

We test for varying:
  - n_quad (number of GJ nodes)
  - distance d of target xₖ from the corner (singularity endpoint of the source panel)
  - φ prescribed as a smooth function (constant, polynomial, trigonometric)
"""

import BoundaryIntegral as BI
using LinearAlgebra, FastGaussQuadrature, HCubature, Printf, StaticArrays

# ── helpers ────────────────────────────────────────────────────────────────

function laplace2d_DT_kernel(src, trg, normal_trg)
    # ∂G/∂n_trg = -1/(2π) · (trg - src)·n / |trg - src|²
    d = trg .- src
    r2 = d[1]^2 + d[2]^2
    return -(d[1]*normal_trg[1] + d[2]*normal_trg[2]) / (2π * r2)
end

function jacobi_poly_vals!(v::AbstractVector{T}, beta::T, t::T) where T
    n = length(v)
    n == 0 && return v
    v[1] = one(T)
    n == 1 && return v
    v[2] = ((beta + 2) * t - beta) / 2
    for m in 1:(n - 2)
        abm  = T(2m) + beta
        mp1  = T(m + 1)
        mbp1 = T(m) + beta + 1
        a_c  = (abm + 1) * (abm + 2) / (2 * mp1 * mbp1)
        b_c  = -(beta^2) * (abm + 1) / (2 * mp1 * mbp1 * abm)
        c_c  = T(m) * (T(m) + beta) * (abm + 2) / (mp1 * mbp1 * abm)
        v[m + 2] = (a_c * t + b_c) * v[m + 1] - c_c * v[m]
    end
    return v
end

function jacobi_norms_sq(beta::T, n::Int) where T
    [T(2)^(beta + 1) / (2 * T(m) + beta + 1) for m in 0:(n - 1)]
end

# Barycentric weights for arbitrary nodes
function bary_weights(x::Vector{T}) where T
    n = length(x)
    w = ones(T, n)
    for j in 1:n
        for k in 1:n
            k == j && continue
            w[j] *= x[j] - x[k]
        end
        w[j] = one(T) / w[j]
    end
    w ./= maximum(abs, w)
    return w
end

function bary_eval(x::Vector{T}, w::Vector{T}, phi::Vector{T}, t::T) where T
    # Barycentric interpolation of phi at t
    for i in 1:length(x)
        isapprox(t, x[i]; atol=1e-14) && return phi[i]
    end
    num = zero(T); den = zero(T)
    for i in 1:length(x)
        wi_t = w[i] / (t - x[i])
        num += wi_t * phi[i]
        den += wi_t
    end
    return num / den
end

# ── main study ──────────────────────────────────────────────────────────────

beta = 0.6703              # singularity exponent (γ - 1 for ε_in=1, ε_out=200, 90° corner)
# Source panel: bottom edge from corner (0,0) to (L_panel, 0)
L_panel = 0.05
a_panel = (0.0, 0.0)
b_panel = (L_panel, 0.0)
mid     = ((a_panel[1] + b_panel[1]) / 2, (a_panel[2] + b_panel[2]) / 2)
half    = ((b_panel[1] - a_panel[1]) / 2, (b_panel[2] - a_panel[2]) / 2)
Lhalf   = L_panel / 2

# Smooth test function φ (the unknown's smooth part)
phi_func(t) = 1.0 + 0.5*t + 0.3*t^2   # smooth polynomial

# Reference tolerance (used as "exact")
REF_RTOL = 1e-14
REF_ATOL = 1e-14

println("="^70)
println("Source panel: [0, $(L_panel)] × {0}, β=$(round(beta, digits=4))")
println("φ(t) = 1 + 0.5t + 0.3t²  (smooth polynomial)")
println("="^70)

# ── Table 1: varying n_quad, fixed distance ──────────────────────────────

println("\n--- Table 1: n_quad sweep, fixed target distance from corner ---")
println("Target: left-side panel node at (0, d) with d = 0.001  (normal = (1,0))")
d_trg = 0.001
xk    = (0.0, d_trg)
nk    = (1.0, 0.0)   # outward normal on left panel

println(@sprintf("%-8s  %-14s  %-14s  %-14s  %-14s", "n_quad", "DirectGJ", "hcub+bary", "hcub+Jacobi", "ref_rtol=1e-14"))

for nq in [4, 6, 8, 10, 12, 16, 20, 24]
    gj_xs_raw, gj_ws_raw = gaussjacobi(nq, 0.0, beta)
    gj_xs = Float64.(gj_xs_raw)
    gj_ws = Float64.(gj_ws_raw)

    # φ values at GJ nodes
    phi_vals = phi_func.(gj_xs)

    # Physical source positions
    s_pts = [(mid[1] + t * half[1], mid[2] + t * half[2]) for t in gj_xs]

    # (A) Direct GJ sum: Σⱼ K(sⱼ, xₖ) · wⱼ · L/2 · φⱼ
    val_direct = sum(laplace2d_DT_kernel(s_pts[j], xk, nk) * gj_ws[j] * Lhalf * phi_vals[j]
                     for j in 1:nq)

    # Barycentric weights
    bw = bary_weights(gj_xs)

    # (B) hcub + barycentric
    pv_buf = Vector{Float64}(undef, nq)
    val_bary, _ = hquadrature(
        t -> begin
            s_pt = (mid[1] + t * half[1], mid[2] + t * half[2])
            K = laplace2d_DT_kernel(s_pt, xk, nk)
            phi_t = bary_eval(gj_xs, bw, phi_vals, t)
            K * (1 + t)^beta * phi_t * Lhalf
        end,
        -1.0, 1.0; rtol=1e-10, atol=1e-14)

    # (C) hcub + Jacobi
    pv = Vector{Float64}(undef, nq)
    h  = jacobi_norms_sq(beta, nq)
    # Jacobi coefficients: aₘ = Σⱼ φⱼ Pₘ(tⱼ) wⱼ / hₘ
    P_at_nodes = [jacobi_poly_vals!(copy(pv), beta, t) for t in gj_xs]  # nq vectors of length nq
    a_coeffs = zeros(nq)
    for m in 1:nq
        for j in 1:nq
            a_coeffs[m] += phi_vals[j] * P_at_nodes[j][m] * gj_ws[j]
        end
        a_coeffs[m] /= h[m]
    end
    pv2 = Vector{Float64}(undef, nq)
    val_jacobi, _ = hquadrature(
        t -> begin
            s_pt = (mid[1] + t * half[1], mid[2] + t * half[2])
            K = laplace2d_DT_kernel(s_pt, xk, nk)
            jacobi_poly_vals!(pv2, beta, t)
            phi_t = dot(a_coeffs, pv2)
            K * (1 + t)^beta * phi_t * Lhalf
        end,
        -1.0, 1.0; rtol=1e-10, atol=1e-14)

    # (REF) high-accuracy reference
    val_ref, _ = hquadrature(
        t -> begin
            s_pt = (mid[1] + t * half[1], mid[2] + t * half[2])
            K = laplace2d_DT_kernel(s_pt, xk, nk)
            phi_t = bary_eval(gj_xs, bw, phi_vals, t)
            K * (1 + t)^beta * phi_t * Lhalf
        end,
        -1.0, 1.0; rtol=REF_RTOL, atol=REF_ATOL)

    err_direct  = abs(val_direct  - val_ref)
    err_bary    = abs(val_bary    - val_ref)
    err_jacobi  = abs(val_jacobi  - val_ref)

    println(@sprintf("%-8d  %-14.4e  %-14.4e  %-14.4e  %-14.4e",
        nq, err_direct, err_bary, err_jacobi, abs(val_ref)))
end

# ── Table 2: varying target distance from corner, fixed n_quad ───────────

println("\n--- Table 2: target distance sweep, n_quad=12 ---")
println("Target: left-side panel node at (0, d), normal = (1,0)")
nq = 12
gj_xs_raw, gj_ws_raw = gaussjacobi(nq, 0.0, beta)
gj_xs    = Float64.(gj_xs_raw)
gj_ws    = Float64.(gj_ws_raw)
phi_vals = phi_func.(gj_xs)
s_pts    = [(mid[1] + t * half[1], mid[2] + t * half[2]) for t in gj_xs]
bw       = bary_weights(gj_xs)
h        = jacobi_norms_sq(beta, nq)
pv_node  = [jacobi_poly_vals!(zeros(nq), beta, t) for t in gj_xs]
a_coeffs = [sum(phi_vals[j] * pv_node[j][m] * gj_ws[j] for j in 1:nq) / h[m] for m in 1:nq]
pv2      = Vector{Float64}(undef, nq)

println(@sprintf("%-12s  %-14s  %-14s  %-14s", "d_target", "DirectGJ", "hcub+bary", "hcub+Jacobi"))

for d_trg in [0.1, 0.05, 0.02, 0.01, 0.005, 0.001, 0.0005, 0.0001]
    xk = (0.0, d_trg)
    nk = (1.0, 0.0)

    val_direct = sum(laplace2d_DT_kernel(s_pts[j], xk, nk) * gj_ws[j] * Lhalf * phi_vals[j]
                     for j in 1:nq)

    val_bary, _ = hquadrature(
        t -> begin
            s_pt = (mid[1] + t * half[1], mid[2] + t * half[2])
            K = laplace2d_DT_kernel(s_pt, xk, nk)
            K * (1 + t)^beta * bary_eval(gj_xs, bw, phi_vals, t) * Lhalf
        end,
        -1.0, 1.0; rtol=1e-10, atol=1e-14)

    val_jacobi, _ = hquadrature(
        t -> begin
            s_pt = (mid[1] + t * half[1], mid[2] + t * half[2])
            K = laplace2d_DT_kernel(s_pt, xk, nk)
            jacobi_poly_vals!(pv2, beta, t)
            K * (1 + t)^beta * dot(a_coeffs, pv2) * Lhalf
        end,
        -1.0, 1.0; rtol=1e-10, atol=1e-14)

    val_ref, _ = hquadrature(
        t -> begin
            s_pt = (mid[1] + t * half[1], mid[2] + t * half[2])
            K = laplace2d_DT_kernel(s_pt, xk, nk)
            K * (1 + t)^beta * bary_eval(gj_xs, bw, phi_vals, t) * Lhalf
        end,
        -1.0, 1.0; rtol=REF_RTOL, atol=REF_ATOL)

    err_direct = abs(val_direct - val_ref)
    err_bary   = abs(val_bary   - val_ref)
    err_jacobi = abs(val_jacobi - val_ref)

    println(@sprintf("%-12.4e  %-14.4e  %-14.4e  %-14.4e",
        d_trg, err_direct, err_bary, err_jacobi))
end

# ── Table 3: n_quad sweep at very small distance (cross-corner stress test) ─

println("\n--- Table 3: n_quad sweep at d=0.0001 (near-corner stress test) ---")
d_trg = 0.0001
xk = (0.0, d_trg)
nk = (1.0, 0.0)
println(@sprintf("%-8s  %-14s  %-14s  %-14s", "n_quad", "DirectGJ", "hcub+bary", "hcub+Jacobi"))

for nq in [4, 6, 8, 10, 12, 16, 20]
    gj_xs_raw, gj_ws_raw = gaussjacobi(nq, 0.0, beta)
    gj_xs    = Float64.(gj_xs_raw)
    gj_ws    = Float64.(gj_ws_raw)
    phi_vals = phi_func.(gj_xs)
    s_pts    = [(mid[1] + t * half[1], mid[2] + t * half[2]) for t in gj_xs]
    bw       = bary_weights(gj_xs)
    h        = jacobi_norms_sq(beta, nq)
    pv_node  = [jacobi_poly_vals!(zeros(nq), beta, t) for t in gj_xs]
    a_coeffs = [sum(phi_vals[j] * pv_node[j][m] * gj_ws[j] for j in 1:nq) / h[m] for m in 1:nq]
    pv2      = Vector{Float64}(undef, nq)

    val_direct = sum(laplace2d_DT_kernel(s_pts[j], xk, nk) * gj_ws[j] * Lhalf * phi_vals[j]
                     for j in 1:nq)

    val_bary, _ = hquadrature(
        t -> begin
            s_pt = (mid[1] + t * half[1], mid[2] + t * half[2])
            K = laplace2d_DT_kernel(s_pt, xk, nk)
            K * (1 + t)^beta * bary_eval(gj_xs, bw, phi_vals, t) * Lhalf
        end,
        -1.0, 1.0; rtol=1e-10, atol=1e-14)

    val_jacobi, _ = hquadrature(
        t -> begin
            s_pt = (mid[1] + t * half[1], mid[2] + t * half[2])
            K = laplace2d_DT_kernel(s_pt, xk, nk)
            jacobi_poly_vals!(pv2, beta, t)
            K * (1 + t)^beta * dot(a_coeffs, pv2) * Lhalf
        end,
        -1.0, 1.0; rtol=1e-10, atol=1e-14)

    val_ref, _ = hquadrature(
        t -> begin
            s_pt = (mid[1] + t * half[1], mid[2] + t * half[2])
            K = laplace2d_DT_kernel(s_pt, xk, nk)
            K * (1 + t)^beta * bary_eval(gj_xs, bw, phi_vals, t) * Lhalf
        end,
        -1.0, 1.0; rtol=REF_RTOL, atol=REF_ATOL)

    err_direct = abs(val_direct - val_ref)
    err_bary   = abs(val_bary   - val_ref)
    err_jacobi = abs(val_jacobi - val_ref)

    println(@sprintf("%-8d  %-14.4e  %-14.4e  %-14.4e",
        nq, err_direct, err_bary, err_jacobi))
end

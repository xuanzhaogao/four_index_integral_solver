#=
Data generation for Figure 4 (near-field quadrature correction at the
panel-pair level).

Two parallel unit panels:
    P = [-1/2, 1/2]^2 x {0},  n_P = e_z   (target)
    Q = [-1/2, 1/2]^2 x {d},              (source, density sigma_Q)

For each d, evaluate D^T_{PQ} sigma_Q at the p x p Gauss-Legendre nodes of P:
  - standard p x p GL on Q
  - dynamically upsampled (p_up >= p, doubling until self-convergence
    < eps_corr) on Q with sigma_Q sampled analytically
  - HCubature reference (rtol = 1e-12, atol = 1e-14)

The upsampled order p_up is chosen exactly as in the production code path
(`check_quad_order3d`):  start at p_up = p, then keep doubling until two
consecutive GL approximations agree to `eps_corr` in the inf-norm over all
target points.  The same p_up is used for every target on P, matching the
per-pair convention in the production solver.

Save fig4_data.jls for fig4_plot.jl.
=#

using LinearAlgebra
using FastGaussQuadrature
using HCubature
using Serialization
using Printf

# ---------------------------------------------------------------------------
# Problem setup
# ---------------------------------------------------------------------------
const p_quad   = 8
const eps_corr = 1e-12             # target tolerance for dynamic p_up
const p_up_max = 512               # safety cap on p_up

const c_near_list = (5, 8)          # integer near-field thresholds

# Panel half-side; both panels are [-1/2, 1/2]^2 in their tangent plane
const L_half = 0.5

const n_P = (0.0, 0.0, 1.0)        # target normal

# Smooth source density on Q; well resolved by p = 8.
σ_Q(η1, η2) = exp(-2.0 * ((η1 - 0.1)^2 + (η2 + 0.15)^2))

# d/h sweep:  h = 1/p
const h_node = 1.0 / p_quad
const dh_list = collect(10 .^ range(log10(0.2), log10(20.0); length = 25))
const d_list  = dh_list .* h_node

const hc_rtol = 1e-12
const hc_atol = 1e-14
const hc_max  = 10_000_000

# ---------------------------------------------------------------------------
# Kernel:  K(x, y) = ((x - y) . n_P) / |x - y|^3  / (4 pi)
# ---------------------------------------------------------------------------
@inline function kernel(x::NTuple{3,Float64}, y::NTuple{3,Float64})
    r = (x[1] - y[1], x[2] - y[2], x[3] - y[3])
    r2 = r[1]*r[1] + r[2]*r[2] + r[3]*r[3]
    invr = 1.0 / sqrt(r2)
    return (r[1]*n_P[1] + r[2]*n_P[2] + r[3]*n_P[3]) * invr^3 / (4π)
end

# Target nodes on P: p x p Gauss-Legendre tensor product
const ns_p, _ws_p = gausslegendre(p_quad)
function target_nodes()
    pts = Vector{NTuple{3,Float64}}(undef, p_quad * p_quad)
    k = 0
    for j in 1:p_quad, i in 1:p_quad
        k += 1
        pts[k] = (ns_p[i] * L_half, ns_p[j] * L_half, 0.0)
    end
    return pts
end

# n x n GL approximation of  int_Q K(x, y) sigma(y) dS_y , sigma analytic.
function gl_quad_apply(x::NTuple{3,Float64}, d::Float64, n::Int)
    ns, ws = gausslegendre(n)
    s = 0.0
    @inbounds for j in 1:n, i in 1:n
        η1 = ns[i] * L_half
        η2 = ns[j] * L_half
        s += ws[i] * ws[j] * σ_Q(η1, η2) * kernel(x, (η1, η2, d))
    end
    return s * L_half * L_half
end

# Vectorized version: returns I(x_t; n) for every target.
function gl_quad_apply_all(targets, d::Float64, n::Int)
    ns, ws = gausslegendre(n)
    out = zeros(length(targets))
    @inbounds for j in 1:n, i in 1:n
        η1 = ns[i] * L_half
        η2 = ns[j] * L_half
        w  = ws[i] * ws[j] * σ_Q(η1, η2)
        y  = (η1, η2, d)
        for (t, x) in pairs(targets)
            out[t] += w * kernel(x, y)
        end
    end
    return out .* (L_half * L_half)
end

# Adaptive upsampling: doubling p_up until ||I(2p) - I(p)||_inf < eps_corr,
# matching the production check_quad_order3d criterion (one p_up per pair).
function adaptive_p_up(targets, d::Float64, eps::Float64, p0::Int, pmax::Int)
    prev = gl_quad_apply_all(targets, d, p0)
    p_try = p0
    while p_try < pmax
        curr = gl_quad_apply_all(targets, d, 2 * p_try)
        if maximum(abs.(curr .- prev)) < eps
            return p_try, prev
        end
        prev = curr
        p_try *= 2
    end
    return pmax, prev
end

# HCubature reference for one target.
function hcubature_ref(x::NTuple{3,Float64}, d::Float64)
    f = η -> σ_Q(η[1], η[2]) * kernel(x, (η[1], η[2], d))
    val, _ = hcubature(f,
                       (-L_half, -L_half), (L_half, L_half);
                       rtol = hc_rtol, atol = hc_atol, maxevals = hc_max)
    return val
end

# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------
const targets = target_nodes()
const Nt = length(targets)

I_std = zeros(Nt, length(d_list))
I_up  = zeros(Nt, length(d_list))
I_ref = zeros(Nt, length(d_list))
p_up_used = zeros(Int, length(d_list))

@info "Running fig4 sweep" p_quad eps_corr Nd=length(d_list)
for (kd, d) in enumerate(d_list)
    t0 = time()
    I_std[:, kd] .= gl_quad_apply_all(targets, d, p_quad)
    p_up_used[kd], I_up_kd = adaptive_p_up(targets, d, eps_corr, p_quad, p_up_max)
    I_up[:, kd] .= I_up_kd
    @inbounds for it in 1:Nt
        I_ref[it, kd] = hcubature_ref(targets[it], d)
    end
    err_std = norm(I_std[:, kd] - I_ref[:, kd]) / norm(I_ref[:, kd])
    err_up  = norm(I_up[:,  kd] - I_ref[:, kd]) / norm(I_ref[:, kd])
    @info @sprintf("  d/h = %7.3f   d = %.4f   p_up = %3d   E_std = %.3e   E_up = %.3e   (%.2fs)",
                   dh_list[kd], d, p_up_used[kd], err_std, err_up, time() - t0)
end

# ---------------------------------------------------------------------------
# Errors and corrected curves
# ---------------------------------------------------------------------------
E_std = [norm(I_std[:, k] - I_ref[:, k]) / norm(I_ref[:, k]) for k in 1:length(d_list)]
E_up  = [norm(I_up[:,  k] - I_ref[:, k]) / norm(I_ref[:, k]) for k in 1:length(d_list)]

E_corr = Dict{Int,Vector{Float64}}()
for c in c_near_list
    Ec = similar(E_std)
    for k in 1:length(d_list)
        Ec[k] = (dh_list[k] <= c) ? E_up[k] : E_std[k]
    end
    E_corr[c] = Ec
end

# Panel (b) data: integer c sweep.
const c_scan = collect(1:12)
max_err_corr_c   = zeros(length(c_scan))
max_resid_far_c  = zeros(length(c_scan))
for (i, c) in enumerate(c_scan)
    near = dh_list .<= c
    far  = dh_list .>  c
    max_err_corr_c[i]  = any(near) ? maximum(E_up[near])  : NaN
    max_resid_far_c[i] = any(far)  ? maximum(E_std[far])  : NaN
end

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out = (
    p_quad         = p_quad,
    eps_corr       = eps_corr,
    p_up_max       = p_up_max,
    p_up_used      = p_up_used,
    h_node         = h_node,
    dh_list        = dh_list,
    d_list         = d_list,
    c_near_list    = collect(c_near_list),
    sigma_label    = "exp(-2 ((y1-0.1)^2 + (y2+0.15)^2))",
    upsample_mode  = "dynamic doubling, analytic sigma on upsampled GL grid",
    hc_rtol        = hc_rtol,
    hc_atol        = hc_atol,
    targets        = targets,
    I_std          = I_std,
    I_up           = I_up,
    I_ref          = I_ref,
    E_std          = E_std,
    E_up           = E_up,
    E_corr         = E_corr,
    c_scan          = c_scan,
    max_err_corr_c  = max_err_corr_c,
    max_resid_far_c = max_resid_far_c,
    geometry = (
        P = "[-1/2, 1/2]^2 x {0}",
        Q = "[-1/2, 1/2]^2 x {d}",
        n_P = n_P,
    ),
)

datapath = joinpath(@__DIR__, "fig4_data.jls")
open(datapath, "w") do io
    serialize(io, out)
end
@info "Saved data" datapath bytes=stat(datapath).size

println("\n=== d/h    p_up   E_std       E_up        E_corr(c=5)   E_corr(c=8) ===")
for k in 1:length(d_list)
    @printf("  %7.3f  %4d   %.3e   %.3e   %.3e     %.3e\n",
            dh_list[k], p_up_used[k],
            E_std[k], E_up[k], E_corr[5][k], E_corr[8][k])
end

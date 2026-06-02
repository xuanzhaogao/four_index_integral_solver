#=
fig4_corr.jl

Companion to fig4_data.jl that adds the *corrected* (upsampled) rule for
fig 4 right panel.  At each d in the sweep, with the same density and
kernel as fig4_data.jl, we evaluate:

  - I_std : standard p × p Gauss-Legendre quadrature (p = p_target).
  - I_up  : upsampled rule with p_up chosen by inverting the
            af-Klinteberg / Trefethen Bernstein bound,
                p_up = ceil(-log(eps) / (2 log rho_min)),
            with rho_min = d + sqrt(1 + d^2)  (worst-case Bernstein
            parameter for the parallel-panel diagnostic with L_half = 1
            and target at the centre of P).
  - I_ref : HCubature reference at rtol = 1e-12.

Saves fig4_corr.jls used by fig4_plot.jl.
=#

using LinearAlgebra
using FastGaussQuadrature
using HCubature
using Serialization
using Printf
using BoundaryIntegral
const BI = BoundaryIntegral

# ---------------------------------------------------------------------------
# Problem setup -- must match fig4_data.jl exactly.
# ---------------------------------------------------------------------------
const L_half = 1.0
σ_Q(η1, η2) = exp(-η1^2 - η2^2)

@inline function kernel(x::NTuple{3,Float64}, y::NTuple{3,Float64})
    r2 = (x[1]-y[1])^2 + (x[2]-y[2])^2 + (x[3]-y[3])^2
    return 1.0 / (4π * sqrt(r2))
end

const p_target = 8
const ε        = 1e-12
const d_list   = collect(10 .^ range(-1.0, 1.0; length = 25))

rho_panel(d) = d + sqrt(1.0 + d^2)

function p_up_at(d::Float64)
    ρ = rho_panel(d)
    p_up = ceil(Int, -log(ε) / (2 * log(ρ)))
    return max(p_up, p_target)
end

# ---------------------------------------------------------------------------
# Density on the p_target grid + barycentric pointwise evaluator.
# ---------------------------------------------------------------------------
const ns_p, ws_p = gausslegendre(p_target)
const λ_p        = BI.gl_barycentric_weights(ns_p, ws_p)

const sigma_p = let M = Matrix{Float64}(undef, p_target, p_target)
    for j in 1:p_target, i in 1:p_target
        M[i, j] = σ_Q(ns_p[i] * L_half, ns_p[j] * L_half)
    end
    M
end

const targets = let T = Vector{NTuple{3,Float64}}(undef, p_target * p_target)
    k = 0
    for j in 1:p_target, i in 1:p_target
        k += 1
        T[k] = (ns_p[i] * L_half, ns_p[j] * L_half, 0.0)
    end
    T
end
const Nt = length(targets)

@inline function sigma_pointwise(u::Float64, v::Float64)
    rx = BI.barycentric_row(ns_p, λ_p, u)
    ry = BI.barycentric_row(ns_p, λ_p, v)
    s = 0.0
    @inbounds for j in 1:p_target, i in 1:p_target
        s += sigma_p[i, j] * rx[i] * ry[j]
    end
    return s
end

# ---------------------------------------------------------------------------
# Per-d evaluators
# ---------------------------------------------------------------------------
function I_std_at_d(d::Float64)
    out = zeros(Nt)
    @inbounds for j in 1:p_target, i in 1:p_target
        η1 = ns_p[i] * L_half
        η2 = ns_p[j] * L_half
        w  = ws_p[i] * ws_p[j] * sigma_p[i, j]
        y  = (η1, η2, d)
        for t in 1:Nt
            out[t] += w * kernel(targets[t], y)
        end
    end
    return out .* (L_half * L_half)
end

function I_up_at_d(d::Float64, p_up::Int)
    ns_u, ws_u = gausslegendre(p_up)
    sigma_u = Matrix{Float64}(undef, p_up, p_up)
    for j in 1:p_up, i in 1:p_up
        sigma_u[i, j] = sigma_pointwise(ns_u[i], ns_u[j])
    end
    out = zeros(Nt)
    @inbounds for j in 1:p_up, i in 1:p_up
        η1 = ns_u[i] * L_half
        η2 = ns_u[j] * L_half
        w  = ws_u[i] * ws_u[j] * sigma_u[i, j]
        y  = (η1, η2, d)
        for t in 1:Nt
            out[t] += w * kernel(targets[t], y)
        end
    end
    return out .* (L_half * L_half)
end

function I_ref_at_d(d::Float64)
    out = zeros(Nt)
    for t in 1:Nt
        x = targets[t]
        f = η -> sigma_pointwise(η[1] / L_half, η[2] / L_half) *
                 kernel(x, (η[1], η[2], d))
        val, _ = hcubature(f,
                           (-L_half, -L_half), (L_half, L_half);
                           rtol = 1e-12, atol = 1e-14, maxevals = 10_000_000)
        out[t] = val
    end
    return out
end

# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------
E_std     = zeros(length(d_list))
E_corr    = zeros(length(d_list))
ρ_used    = zeros(length(d_list))
p_up_used = zeros(Int, length(d_list))

for (k, d) in enumerate(d_list)
    t0  = time()
    p_up = p_up_at(d)
    ρ    = rho_panel(d)
    p_up_used[k] = p_up
    ρ_used[k]    = ρ

    I_s = I_std_at_d(d)
    I_c = I_up_at_d(d, p_up)
    I_r = I_ref_at_d(d)
    n_r = norm(I_r)
    E_std[k]  = norm(I_s - I_r) / n_r
    E_corr[k] = norm(I_c - I_r) / n_r

    @info @sprintf("  d=%.4f  rho=%.3f  p_up=%d  E_std=%.3e  E_corr=%.3e  (%.2fs)",
                   d, ρ, p_up, E_std[k], E_corr[k], time() - t0)
end

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out = (
    p_target = p_target,
    d_list   = d_list,
    ρ        = ρ_used,
    p_up     = p_up_used,
    E_std    = E_std,
    E_corr   = E_corr,
    ε        = ε,
    L_half   = L_half,
)

datapath = joinpath(@__DIR__, "fig4_corr.jls")
open(datapath, "w") do io
    serialize(io, out)
end
@info "Saved" datapath bytes=stat(datapath).size

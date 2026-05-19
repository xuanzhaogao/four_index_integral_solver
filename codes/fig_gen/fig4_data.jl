#=
Data generation for Figure 4.

Two parallel unit panels:
    P = [-1/2, 1/2]^2 x {0},  n_P = e_z   (target)
    Q = [-1/2, 1/2]^2 x {d},              (source, density sigma_Q)

For each p ∈ {4, 6, 8}, sweep d/h ∈ [0.2, 20] and record the relative L2
error of standard p×p Gauss-Legendre quadrature over Q against the
polynomial-σ_p interpolant ground truth (HCubature, rtol=1e-12). The
production solver only sees σ at the p_quad GL grid, so the "true" integrand
on Q is the tensor barycentric Lagrange interpolant P_σ, and the panel-pair
integral is the integral of K · P_σ over Q.

Save fig4_data.jls for fig4_plot.jl.
=#

using LinearAlgebra
using FastGaussQuadrature
using HCubature
using Serialization
using Printf
using BoundaryIntegral
const BI = BoundaryIntegral

# ---------------------------------------------------------------------------
# Problem setup
# ---------------------------------------------------------------------------
const L_half  = 0.5
const n_P     = (0.0, 0.0, 1.0)
σ_Q(η1, η2)   = exp(-η1^2 - η2^2)

@inline function kernel(x::NTuple{3,Float64}, y::NTuple{3,Float64})
    r  = (x[1] - y[1], x[2] - y[2], x[3] - y[3])
    r2 = r[1]*r[1] + r[2]*r[2] + r[3]*r[3]
    invr = 1.0 / sqrt(r2)
    return (r[1]*n_P[1] + r[2]*n_P[2] + r[3]*n_P[3]) * invr^3 / (4π)
end

const p_values = [4, 6, 8]
const dh_list  = collect(10 .^ range(log10(0.2), log10(20.0); length = 25))
const hc_rtol  = 1e-12
const hc_atol  = 1e-14
const hc_max   = 10_000_000

# ---------------------------------------------------------------------------
# Per-p sweep
# ---------------------------------------------------------------------------
function sweep_p(p::Int)
    h_node = 1.0 / p
    d_list = dh_list .* h_node

    ns_p, ws_p = gausslegendre(p)
    λ_p        = BI.gl_barycentric_weights(ns_p, ws_p)

    sigma_p = Matrix{Float64}(undef, p, p)
    for j in 1:p, i in 1:p
        sigma_p[i, j] = σ_Q(ns_p[i] * L_half, ns_p[j] * L_half)
    end

    targets = Vector{NTuple{3,Float64}}(undef, p * p)
    k = 0
    for j in 1:p, i in 1:p
        k += 1
        targets[k] = (ns_p[i] * L_half, ns_p[j] * L_half, 0.0)
    end
    Nt = length(targets)

    @inline function sigma_pointwise(u::Float64, v::Float64)
        rx = BI.barycentric_row(ns_p, λ_p, u)
        ry = BI.barycentric_row(ns_p, λ_p, v)
        s = 0.0
        @inbounds for j in 1:p, i in 1:p
            s += sigma_p[i, j] * rx[i] * ry[j]
        end
        return s
    end

    function I_std_at_d(d::Float64)
        out = zeros(Nt)
        @inbounds for j in 1:p, i in 1:p
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

    function I_ref_at_d(d::Float64)
        out = zeros(Nt)
        for t in 1:Nt
            x = targets[t]
            f = η -> sigma_pointwise(η[1] / L_half, η[2] / L_half) *
                    kernel(x, (η[1], η[2], d))
            val, _ = hcubature(f,
                               (-L_half, -L_half), (L_half, L_half);
                               rtol = hc_rtol, atol = hc_atol, maxevals = hc_max)
            out[t] = val
        end
        return out
    end

    E_std = zeros(length(d_list))
    for (k, d) in enumerate(d_list)
        t0 = time()
        I_s = I_std_at_d(d)
        I_r = I_ref_at_d(d)
        E_std[k] = norm(I_s - I_r) / norm(I_r)
        @info @sprintf("  p=%d  d/h=%7.3f  d=%.4f  E_std=%.3e  (%.2fs)",
                       p, dh_list[k], d, E_std[k], time() - t0)
    end
    return (p = p, h_node = h_node, d_list = d_list, E_std = E_std)
end

# ---------------------------------------------------------------------------
# Run sweep
# ---------------------------------------------------------------------------
results = NamedTuple[]
for p in p_values
    @info "================= p = $p ================="
    push!(results, sweep_p(p))
end

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out = (
    p_values    = p_values,
    dh_list     = dh_list,
    sweeps      = results,
    L_half      = L_half,
    n_P         = n_P,
    sigma_label = "exp(-y1^2 - y2^2)",
    hc_rtol     = hc_rtol,
    hc_atol     = hc_atol,
)

datapath = joinpath(@__DIR__, "fig4_data.jls")
open(datapath, "w") do io
    serialize(io, out)
end
@info "Saved data" datapath bytes=stat(datapath).size

println("\n=== d/h   ", join((@sprintf("E_std(p=%d)", p) for p in p_values), "    "), " ===")
for k in 1:length(dh_list)
    line = @sprintf("  %7.3f", dh_list[k])
    for r in results
        line *= @sprintf("   %.3e", r.E_std[k])
    end
    println(line)
end

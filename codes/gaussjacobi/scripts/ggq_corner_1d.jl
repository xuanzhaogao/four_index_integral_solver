# Generalized-Gauss quadrature driven by corner-singularity series.
#
# Step 1: determine the singularity-power series {γ_k} for a 2D dielectric
#         wedge from the roots of the angular-ODE transmission det
#         (van Bladel / multicorner2dpowers.jl).
# Step 2: build the GGQ rule on [0,1] for weight x^{α_1}, α_k = γ_k - 1,
#         with Müntz basis exponents γ_k - γ_1 = α_k - α_1.
# Step 3: validate on a BIE-type 1D integral
#           I = ∫_0^L K(s, t) σ(s) ds,  σ(s)=s^{α_1} u(s)
#         where K(s,t) = 1/(2π) · t/((s-x0)^2+t^2) is the normal-derivative
#         (D^T) of the 2D Laplace kernel for a target a small distance t off
#         the wedge edge.

include(joinpath(@__DIR__, "..", "src", "gram.jl"))

using LinearAlgebra
using Roots
using FastGaussQuadrature
using Printf
using Random

setprecision(BigFloat, 256)

# ---------- Step 1: series from the angular ODE -----------------------------

"""
    theta_ODE_det(angles, eps, g)

Determinant for the multi-junction θ-ODE periodic eigenvalue problem.
Roots in g are the admissible singularity powers γ at the corner.
`angles` is length nm-1 (last angle = 2π - sum), `eps` length nm.
"""
function theta_ODE_det(a::AbstractVector, e::AbstractVector, g)
    M = Matrix(1.0I, 2, 2)
    aa = [a; 2pi - sum(a)]
    @assert aa[end] >= -1e-12 "angles sum > 2π"
    nm = length(e)
    @assert length(aa) == nm
    nexte = circshift(e, -1)
    for (j, ej) in enumerate(e)
        if !isinf(ej)
            c = cos(g * aa[j])
            s = sin(g * aa[j])
            soverg = g == 0.0 ? aa[j] : s / g
            M = [1.0 0.0; 0.0 ej/nexte[j]] * [c soverg; -g*s c] * M
        end
    end
    return isinf(e[end]) ? M[1, 2] : det(M - I)
end

"""
    corner_power_series(angles, eps; n_powers, gmax, ngrid)

Return the first `n_powers` positive roots of g -> theta_ODE_det(angles,eps,g)
on (0, gmax], found by sign-change bracketing on a fine grid.
"""
function corner_power_series(angles::AbstractVector, eps::AbstractVector;
                             n_powers::Int = 10, gmax::Real = 30.0,
                             ngrid::Int = 50_000)
    f = g -> theta_ODE_det(angles, eps, g)
    gs = range(1e-6, gmax, length=ngrid)
    fs = [f(g) for g in gs]
    roots = Float64[]
    for i in 1:length(gs)-1
        if !isnan(fs[i]) && !isnan(fs[i+1]) && fs[i] * fs[i+1] < 0
            try
                r = find_zero(f, (gs[i], gs[i+1]), Bisection())
                if isempty(roots) || r - roots[end] > 1e-6
                    push!(roots, r)
                end
            catch
            end
        end
        length(roots) >= n_powers && break
    end
    return roots
end

# ---------- Step 2: GGQ from the series --------------------------------------

"""
    build_ggq(γ_target, m; T=BigFloat, nsteps=30)

Return (x, w) — the m-point GGQ on [0,1] with weight x^{a0}, a0 = γ_target[1]-1,
and Müntz basis exponents [γ_target[k] - γ_target[1] for k in 1..2m].
The first basis exponent is forced to 0 so the standard moment problem applies.
"""
function build_ggq(γ_target::AbstractVector, m::Int;
                   T = BigFloat, nsteps::Int = 30, verbose::Bool = false)
    @assert length(γ_target) >= 2m "Need at least 2m powers; got $(length(γ_target)) for m=$m"
    a0 = γ_target[1] - 1                         # weight power: x^{γ_1-1}
    γshift = γ_target[1:2m] .- γ_target[1]       # basis exponents starting at 0
    x, w = generalized_gauss_continuation(T(a0), T.(γshift), m;
                                          nsteps = nsteps, verbose = verbose)
    return Float64.(x), Float64.(w), Float64(a0), Float64.(γshift)
end

# ---------- Step 3: validation on a BIE-style integral -----------------------

# Density and reference integral on [0, L] with target at (x0, t).
# σ(s) = s^{α_1} u(s) where u is a finite combination of the GGQ-basis exponents
# (the "smooth" part u corresponds to the unknown solved for in a BIE).
function make_density(γshift::AbstractVector{Float64}, α1::Float64; seed::Int = 7)
    Nb = length(γshift)
    rng = MersenneTwister(seed)
    coefs = [randn(rng) / (1 + (k-1)) for k in 1:Nb]
    u(s) = sum(coefs[k] * s^γshift[k] for k in 1:Nb)
    σ(s) = s^α1 * u(s)
    return u, σ, coefs
end

# Reference: high-order Gauss-Jacobi for ∫_0^L s^α1 u(s) f(s) ds
function ref_integral(σ_or_u_f::Function, α1::Float64, L::Float64; Nref::Int = 400)
    # ∫_0^L s^α1 g(s) ds where g(s) = u(s) f(s).
    tk, Wk = gaussjacobi(Nref, 0.0, α1)            # Jacobi on [-1,1] with weight (1-x)^0 (1+x)^α1 ?
    # gaussjacobi(N, α, β) integrates ∫_{-1}^1 (1-x)^α (1+x)^β g(x) dx.
    # We want ∫_0^L s^α1 g(s) ds. Change var s = L*(1+x)/2, ds = L/2 dx:
    # = (L/2)^{α1+1} ∫_{-1}^1 (1+x)^α1 g(L(1+x)/2) dx
    s_pts = 0.5 .* L .* (1 .+ tk)
    return (L/2)^(α1+1) * sum(Wk[k] * σ_or_u_f(s_pts[k]) for k in eachindex(tk))
end

# Apply GGQ rule on [0, L]: ∫_0^L s^α u(s)f(s) ds ≈ L^{α+1} Σ_j w_j u(L x_j) f(L x_j).
function ggq_apply(x01::Vector{Float64}, w01::Vector{Float64},
                   α1::Float64, L::Float64, integrand::Function)
    s = L .* x01
    return L^(α1+1) * sum(w01[j] * integrand(s[j]) for j in 1:length(x01))
end

# Standard Gauss-Jacobi on [0,L] for weight s^α1.
function gj_apply(m::Int, α1::Float64, L::Float64, integrand::Function)
    tk, Wk = gaussjacobi(m, 0.0, α1)
    s = 0.5 .* L .* (1 .+ tk)
    return (L/2)^(α1+1) * sum(Wk[k] * integrand(s[k]) for k in eachindex(tk))
end

# Plain Gauss-Legendre on [0,L] folding s^α1 into the integrand.
function gl_apply(m::Int, α1::Float64, L::Float64, integrand::Function)
    tk, Wk = gausslegendre(m)
    s = 0.5 .* L .* (1 .+ tk)
    return 0.5 * L * sum(Wk[k] * s[k]^α1 * integrand(s[k]) for k in eachindex(tk))
end

# ---------- Driver -----------------------------------------------------------

# Right-angle dielectric corner: interior wedge angle α=π/2, eps_in=1, eps_out=200,
# matching the test problems in scripts/density_comparison.jl and benchmark_convergence.jl.
function main()
    eps_in  = 1.0
    eps_out = 200.0
    α_int   = pi/2     # interior wedge angle (eps_in occupies α_int)

    # multicorner2dpowers convention: angles vector has length nm-1.
    angles = [α_int]
    eps    = [eps_in, eps_out]
    n_powers = 16

    println("="^70)
    println("Step 1: corner singularity series for α=π/2, ε_in=1, ε_out=200")
    println("="^70)
    γs = corner_power_series(angles, eps; n_powers = n_powers, gmax = 25.0)
    @printf("Found %d powers γ_k:\n", length(γs))
    for (k, g) in enumerate(γs)
        @printf("  γ_%-2d = %.10f   (α_k = γ-1 = %.10f)\n", k, g, g-1)
    end

    α1 = γs[1] - 1
    γshift = γs .- γs[1]   # density-power shifts; α_k - α_1 = γ_k - γ_1
    println("\nWeight power on [0,1] : a0 = γ_1 - 1 = $(@sprintf("%.10f", α1))")
    println("Basis exponent shifts : γ_k - γ_1 = $(round.(γshift, digits=4))")

    # --------- Step 2: build GGQ for several m, sanity-check moments -------
    println("\n" * "="^70)
    println("Step 2: build GGQ rules from Müntz basis")
    println("="^70)
    rules = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
    for m in 2:8
        try
            x, w, a0, γs_basis = build_ggq(γs, m; nsteps = 40)
            rules[m] = (x, w)
            # consistency check: sum of weights = ∫_0^1 x^{a0} dx = 1/(a0+1)
            err_norm = abs(sum(w) - 1/(a0+1))
            @printf("  m=%d  nodes ∈ [%.3e, %.3e]  Σw err = %.2e  min(w)=%.2e\n",
                    m, minimum(x), maximum(x), err_norm, minimum(w))
        catch e
            @printf("  m=%d FAILED: %s\n", m, sprint(showerror, e))
        end
    end

    # --------- Step 3: BIE-style integral test ----------------------------
    println("\n" * "="^70)
    println("Step 3: BIE integral test")
    println("    I = ∫_0^L K(s,t) σ(s) ds,   σ(s) = s^{α_1} u(s)")
    println("    K(s,t) = (1/2π) · t / ((s-x0)^2 + t^2)   (D^T-style)")
    println("="^70)

    L  = 0.05                    # corner-panel length (matches l_corner≈0.05)
    x0 = -1e-3                   # target s-coordinate (just below corner, on adj edge)
    t  = 1e-3                    # target distance off the edge
    K(s) = (1/(2π)) * t / ((s - x0)^2 + t^2)

    Nb = 2 * maximum(keys(rules))   # use as many basis modes as the largest m supports
    γshift_full = γs[1:Nb] .- γs[1]
    rng = MersenneTwister(7)
    coefs = [randn(rng) / (1 + (k-1)) for k in 1:Nb]
    u(s) = sum(coefs[k] * s^γshift_full[k] for k in 1:Nb)
    integrand(s) = u(s) * K(s)
    σf(s) = s^α1 * u(s) * K(s)         # for ref_integral helper

    I_ref = ref_integral(s -> u(s)*K(s), α1, L; Nref = 400)
    @printf("Reference (400-pt Gauss-Jacobi)  I = %.16e\n\n", I_ref)

    @printf("%-4s  %-14s  %-14s  %-14s\n", "m", "GGQ err", "GJ err", "GL err")
    println("-"^54)
    for m in sort(collect(keys(rules)))
        x01, w01 = rules[m]
        I_ggq = ggq_apply(x01, w01, α1, L, integrand)
        I_gj  = gj_apply(m, α1, L, integrand)
        I_gl  = gl_apply(m, α1, L, integrand)
        @printf("%-4d  %.4e      %.4e      %.4e\n",
                m, abs(I_ggq - I_ref), abs(I_gj - I_ref), abs(I_gl - I_ref))
    end

    # Also report the "smooth K only" case (no density singularity, u≡1) as a
    # control — pure Gauss-Jacobi should already be excellent here.
    println("\nControl: u≡1 (only the s^{α_1} weight matters; no Müntz tail)")
    I_ref2 = ref_integral(K, α1, L; Nref = 400)
    @printf("%-4s  %-14s  %-14s  %-14s\n", "m", "GGQ err", "GJ err", "GL err")
    for m in sort(collect(keys(rules)))
        x01, w01 = rules[m]
        I_ggq = ggq_apply(x01, w01, α1, L, K)
        I_gj  = gj_apply(m, α1, L, K)
        I_gl  = gl_apply(m, α1, L, K)
        @printf("%-4d  %.4e      %.4e      %.4e\n",
                m, abs(I_ggq - I_ref2), abs(I_gj - I_ref2), abs(I_gl - I_ref2))
    end

    return rules, γs
end

if abspath(PROGRAM_FILE) == @__FILE__
    rules, γs = main()
end

#=
Data generation for Figure 5 (TKM volume-solver validation, Coulomb kernel).

We evaluate
    u(x) = (1 / 4 pi) int rho(y) / |x - y| dy
for a normalized 3D Gaussian rho on the box B = [-1, 1]^3 by sampling rho on a
uniform Cartesian grid and passing the samples to two competing methods:
  - TKM   (truncated kernel method)        — `ltkm3dc` from TKM3D
  - FMM   (direct discrete sum at eps tol) — `lfmm3d`  from FMM3D

The analytic reference is
    u_ref(x) = (1 / 4 pi) erf(|x| / (sqrt(2) s)) / |x|,
with limit  (1 / 4 pi) sqrt(2 / pi) / s  at x = 0.

Panels:
  (a) relative L2 error vs grid resolution n for several tolerances eps
  (b) relative L2 error vs eta = min_alpha P_alpha / (l_alpha + L)
      with n fixed; uses a small standalone TKM evaluator so we can drive eta
      below 1 (production TKM enforces eta >= 1).

All errors are measured on a fixed set of random off-grid targets in B, so
the FMM has no self-interaction and the same targets are used at every n.
=#

using LinearAlgebra
using SpecialFunctions
using Serialization
using Printf
using Random
using FINUFFT
using TKM3D
using FMM3D

# ---------------------------------------------------------------------------
# Problem setup
# ---------------------------------------------------------------------------
const s_gauss  = 0.10
const y0       = (0.0, 0.0, 0.0)

const box_half = 1.0
const l_box    = 2 * box_half
const L_diag   = sqrt(3) * l_box

const eps_list = (1e-3, 1e-6, 1e-9, 1e-12)
const n_list   = collect(8:8:64)

const n_for_eta = 64
const eta_list  = collect(10 .^ range(log10(0.5), log10(2.5); length = 25))

# Fixed off-grid target set (same across all n and methods).
const N_targets = 200
Random.seed!(42)
let
    global targets = Matrix{Float64}(undef, 3, N_targets)
    for k in 1:N_targets
        # uniform in [-0.9, 0.9]^3 (well inside B, avoids cell-center alignment)
        targets[1, k] = 1.8 * (rand() - 0.5)
        targets[2, k] = 1.8 * (rand() - 0.5)
        targets[3, k] = 1.8 * (rand() - 0.5)
    end
end

# ---------------------------------------------------------------------------
# Density and analytic reference
# ---------------------------------------------------------------------------
@inline function rho(y::NTuple{3,Float64})
    r2 = (y[1]-y0[1])^2 + (y[2]-y0[2])^2 + (y[3]-y0[3])^2
    return exp(-r2 / (2 * s_gauss^2)) / (2π * s_gauss^2)^(3/2)
end

@inline function u_ref(x::NTuple{3,Float64})
    r = sqrt((x[1]-y0[1])^2 + (x[2]-y0[2])^2 + (x[3]-y0[3])^2)
    if r < 1e-14
        return sqrt(2 / π) / s_gauss / (4π)
    end
    return erf(r / (sqrt(2) * s_gauss)) / r / (4π)
end

function tail_mass(R::Float64, s::Float64)
    a = R / (sqrt(2) * s)
    return erfc(a) + sqrt(2/π) * (R/s) * exp(-R^2 / (2 * s^2))
end

const u_r_t = [u_ref((targets[1,k], targets[2,k], targets[3,k])) for k in 1:N_targets]
const u_r_norm = norm(u_r_t)

# ---------------------------------------------------------------------------
# Build cell-centered grid sources
# ---------------------------------------------------------------------------
function grid_sources(n::Int)
    h  = l_box / n
    xs = collect(-box_half + h/2 .+ h .* (0:n-1))
    sources = Matrix{Float64}(undef, 3, n^3)
    charges = Vector{Float64}(undef, n^3)
    k = 0
    @inbounds for kz in 1:n, ky in 1:n, kx in 1:n
        k += 1
        sources[1, k] = xs[kx]
        sources[2, k] = xs[ky]
        sources[3, k] = xs[kz]
        charges[k]    = rho((xs[kx], xs[ky], xs[kz])) * h^3
    end
    return sources, charges
end

# ---------------------------------------------------------------------------
# Standalone continuous-TKM evaluator with a user-controlled Δk.
# eta = P / (l + L) where l = 2 (box side) and L = sqrt(3) * l.
# eta = 1 is the tightest allowed real-space period; eta < 1 -> aliasing.
# ---------------------------------------------------------------------------
function tkm_eval_eta(sources::Matrix{Float64}, charges::Vector{Float64},
                     targets::Matrix{Float64}, eta::Float64;
                     eps::Float64 = 1e-12)
    l_x = l_y = l_z = l_box
    L   = L_diag

    Δk = 2π / (l_x + L) / eta

    kmax_use = Float64(estimate_kcut3dc(sources;
                                        charges = charges,
                                        tol = eps, eps = eps).kcut)

    kx = TKM3D.centered_mode_axis(Δk, kmax_use)
    ky = kx
    kz = kx

    srcx = Δk .* view(sources, 1, :)
    srcy = Δk .* view(sources, 2, :)
    srcz = Δk .* view(sources, 3, :)
    trgx = Δk .* view(targets, 1, :)
    trgy = Δk .* view(targets, 2, :)
    trgz = Δk .* view(targets, 3, :)

    coeff = nufft3d1(srcx, srcy, srcz, complex.(charges), -1, eps,
                     length(kx), length(ky), length(kz))
    if ndims(coeff) == 4 && size(coeff, 4) == 1
        coeff = dropdims(coeff; dims = 4)
    end

    @inbounds for iz in eachindex(kz), iy in eachindex(ky), ix in eachindex(kx)
        k = sqrt(kx[ix]^2 + ky[iy]^2 + kz[iz]^2)
        if k <= kmax_use
            coeff[ix, iy, iz] *= TKM3D.truncated_laplace3d_hat(k, L)
        else
            coeff[ix, iy, iz] = zero(eltype(coeff))
        end
    end

    prefactor = (Δk * Δk * Δk) / (2π)^3
    pot = nufft3d2(trgx, trgy, trgz, 1, eps, coeff)
    return prefactor .* real.(pot)
end

# ---------------------------------------------------------------------------
# Panel (a): convergence vs n  —  TKM and FMM, several eps
# ---------------------------------------------------------------------------
err_tkm_n = Dict{Float64,Vector{Float64}}()
err_fmm_n = Dict{Float64,Vector{Float64}}()
for eps in eps_list
    err_tkm_n[eps] = Float64[]
    err_fmm_n[eps] = Float64[]
end

@info "Panel (a) — convergence vs n"
for n in n_list
    src, q = grid_sources(n)
    for eps in eps_list
        # TKM
        vals_tkm = ltkm3dc(eps, src; charges = q, targets = targets, pgt = 1)
        u_tkm    = vals_tkm.pottarg
        e_tkm    = norm(u_tkm .- u_r_t) / u_r_norm
        push!(err_tkm_n[eps], e_tkm)

        # FMM (kernel convention: 1/|x-y|, so divide by 4π)
        vals_fmm = lfmm3d(eps, src; charges = q, targets = targets, pgt = 1)
        u_fmm    = vals_fmm.pottarg ./ (4π)
        e_fmm    = norm(u_fmm .- u_r_t) / u_r_norm
        push!(err_fmm_n[eps], e_fmm)

        @info @sprintf("  n = %3d   eps = %.0e   E_TKM = %.3e   E_FMM = %.3e",
                       n, eps, e_tkm, e_fmm)
    end
end

# ---------------------------------------------------------------------------
# Panel (b): aliasing / padding test — vary eta at fixed n, several eps
# ---------------------------------------------------------------------------
err_tkm_eta = Dict{Float64,Vector{Float64}}()
for eps in eps_list
    err_tkm_eta[eps] = Float64[]
end

@info "Panel (b) — eta sweep at n = $n_for_eta"
src_c, q_c = grid_sources(n_for_eta)
for eps in eps_list
    for eta in eta_list
        u_t = tkm_eval_eta(src_c, q_c, targets, eta; eps = eps)
        e   = norm(u_t .- u_r_t) / u_r_norm
        push!(err_tkm_eta[eps], e)
    end
    @info @sprintf("  eps = %.0e :  E(eta=0.5) = %.3e  E(eta=1) ~ %.3e  E(eta=2.5) = %.3e",
                   eps,
                   err_tkm_eta[eps][1],
                   err_tkm_eta[eps][argmin(abs.(eta_list .- 1.0))],
                   err_tkm_eta[eps][end])
end

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out = (
    s_gauss     = s_gauss,
    y0          = y0,
    box_half    = box_half,
    l_box       = l_box,
    L_diag      = L_diag,
    eps_list    = collect(eps_list),
    n_list      = n_list,
    err_tkm_n   = err_tkm_n,
    err_fmm_n   = err_fmm_n,
    n_for_eta   = n_for_eta,
    eta_list    = eta_list,
    err_tkm_eta = err_tkm_eta,
    N_targets   = N_targets,
    tail_mass   = tail_mass(box_half, s_gauss),
)

datapath = joinpath(@__DIR__, "fig5_data.jls")
open(datapath, "w") do io
    serialize(io, out)
end
@info "Saved data" datapath bytes=stat(datapath).size

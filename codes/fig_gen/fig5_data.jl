#=
Data generation for Figure 5 (TKM/FMM hybrid incident-potential validation).

We evaluate
    u(x) = (1 / 4 pi) int rho(y) / |x - y| dy
for a normalized 3D Gaussian rho on B = [-1, 1]^3, sampled on a cell-centered
uniform Cartesian grid, and compare two evaluators at the SAME targets:
  - the Section 3 hybrid: TKM inside B_pad, FMM outside (PrecomputedVolumeField)
  - the pure particle sum of Eq. (3.13), at every target

The padding is h_n = c_pad * h with c_pad = 5, so both panels exercise the
near/far switch. Targets span both regions and keep their classification for
every n in the sweep:
  near: 200 points in [-0.9, 0.9]^3   -- inside B_pad even at n = 64 (B_pad ~ 1.156)
  far:  200 points with max|coord| in [2.5, 4] -- outside B_pad even at n = 8 (B_pad ~ 2.25)

The analytic reference is u_ref(x) = (1 / 4 pi) erf(|x| / (sqrt(2) s)) / |x|,
with limit (1 / 4 pi) sqrt(2 / pi) / s at x = 0, valid at every target.

Panels:
  (a) relative L2 error vs grid resolution n, for several tolerances tau
  (b) relative L2 error vs eta = min_alpha L_alpha / (l_alpha + h_n + L), n fixed.
      Uses a standalone TKM+FMM evaluator with a tunable dk, because the
      production path now sits exactly at eta = 1 and cannot be driven below it.
=#

using LinearAlgebra
using SpecialFunctions
using Serialization
using Printf
using Random
using FMM3D
using TKM3D
using BoundaryIntegral
const BI = BoundaryIntegral

# ---------------------------------------------------------------------------
# Problem setup
# ---------------------------------------------------------------------------
const s_gauss = 0.10
const y0      = (0.0, 0.0, 0.0)

const box_half = 1.0
const l_box    = 2 * box_half
const C_PAD    = 5.0

const eps_list = (1e-3, 1e-6, 1e-9, 1e-12)
const n_list   = collect(8:8:64)

const n_for_eta = 64
const eta_list  = collect(10 .^ range(log10(0.5), log10(1.6); length = 20))

# ---------------------------------------------------------------------------
# Fixed target set: near + far, classification independent of n
# ---------------------------------------------------------------------------
const N_near = 200
const N_far  = 200

function build_targets()
    Random.seed!(42)
    t = Matrix{Float64}(undef, 3, N_near + N_far)
    # near: uniform in [-0.9, 0.9]^3, inside B_pad for every n (B_pad >= 1.156)
    for k in 1:N_near
        t[1, k] = 1.8 * (rand() - 0.5)
        t[2, k] = 1.8 * (rand() - 0.5)
        t[3, k] = 1.8 * (rand() - 0.5)
    end
    # far: max|coord| in [2.5, 4], outside B_pad for every n (B_pad <= 2.25)
    k = N_near
    while k < N_near + N_far
        p = 8.0 .* (rand(3) .- 0.5)          # uniform in [-4, 4]^3
        m = maximum(abs, p)
        (m >= 2.5 && m <= 4.0) || continue
        k += 1
        t[1, k] = p[1]; t[2, k] = p[2]; t[3, k] = p[3]
    end
    return t
end

const targets  = build_targets()
const N_targets = size(targets, 2)
const near_rng = 1:N_near
const far_rng  = (N_near + 1):(N_near + N_far)

# ---------------------------------------------------------------------------
# Density, source construction, analytic reference
# ---------------------------------------------------------------------------
@inline function rho(y::NTuple{3,Float64})
    r2 = (y[1]-y0[1])^2 + (y[2]-y0[2])^2 + (y[3]-y0[3])^2
    return exp(-r2 / (2 * s_gauss^2)) / (2π * s_gauss^2)^(3/2)
end

@inline function u_ref(x::NTuple{3,Float64})
    r = sqrt((x[1]-y0[1])^2 + (x[2]-y0[2])^2 + (x[3]-y0[3])^2)
    r < 1e-14 && return sqrt(2 / π) / s_gauss / (4π)
    return erf(r / (sqrt(2) * s_gauss)) / r / (4π)
end

tail_mass(R, s) = erfc(R / (sqrt(2)*s)) + sqrt(2/π) * (R/s) * exp(-R^2 / (2*s^2))

const u_r_t    = [u_ref((targets[1,k], targets[2,k], targets[3,k])) for k in 1:N_targets]
const u_r_norm = norm(u_r_t)

relerr(u) = norm(u .- u_r_t) / u_r_norm
relerr(u, rng) = norm(u[rng] .- u_r_t[rng]) / norm(u_r_t[rng])

# Cell-centered grid on B = [-1,1]^3. The identity-basis grid constructor stores
# A_rho = h * I, so BI.source_box returns exactly [-1,1]^3 and BI.lattice_spacing
# returns exactly h.
function grid_source(n::Int)
    h  = l_box / n
    xs = collect(-box_half + h/2 .+ h .* (0:n-1))
    weights = fill(h^3, n, n, n)
    density = Array{Float64,3}(undef, n, n, n)
    for k in 1:n, j in 1:n, i in 1:n
        density[i,j,k] = rho((xs[i], xs[j], xs[k]))
    end
    return VolumeSource((xs, xs, xs), weights, density)
end

# ---------------------------------------------------------------------------
# Standalone hybrid with a tunable Fourier spacing (panel b).
# Mirrors BI.near_field_geometry, then scales dk by 1/eta so eta < 1 is reachable.
# ---------------------------------------------------------------------------
function hybrid_eval_eta(vs::VolumeSource{Float64,3}, trg::Matrix{Float64},
                        eta::Float64; eps::Float64)
    g   = BI.near_field_geometry(vs; c_pad = C_PAD)
    src, q = BI._volume_source_fmm_sources(vs)
    km  = BI._estimate_tkm3dc_kmax(BI._estimate_source_spacing(vs))
    dk  = ntuple(d -> 2π / (g.l[d] + g.hn + g.L) / eta, 3)

    kx = TKM3D.centered_mode_axis(dk[1], km)
    ky = TKM3D.centered_mode_axis(dk[2], km)
    kz = TKM3D.centered_mode_axis(dk[3], km)

    out  = Vector{Float64}(undef, size(trg, 2))
    inb  = [BI.in_near_region(g, trg, i) for i in 1:size(trg, 2)]
    bidx = findall(inb); oidx = findall(!, inb)

    if !isempty(bidx)
        srcx = dk[1] .* (vec(view(src,1,:)) .- g.center[1])
        srcy = dk[2] .* (vec(view(src,2,:)) .- g.center[2])
        srcz = dk[3] .* (vec(view(src,3,:)) .- g.center[3])
        coeff0 = TKM3D.FINUFFT.nufft3d1(srcx, srcy, srcz, complex.(q), -1, eps,
                                        length(kx), length(ky), length(kz))
        coeff = ndims(coeff0) == 4 ? dropdims(coeff0; dims = 4) : coeff0
        @inbounds for iz in eachindex(kz), iy in eachindex(ky), ix in eachindex(kx)
            k = sqrt(kx[ix]^2 + ky[iy]^2 + kz[iz]^2)
            coeff[ix,iy,iz] = k <= km ?
                coeff[ix,iy,iz] * TKM3D.truncated_laplace3d_hat(k, g.L) :
                zero(eltype(coeff))
        end
        txn = Float64[dk[1] * (trg[1,i] - g.center[1]) for i in bidx]
        tyn = Float64[dk[2] * (trg[2,i] - g.center[2]) for i in bidx]
        tzn = Float64[dk[3] * (trg[3,i] - g.center[3]) for i in bidx]
        vals = TKM3D._finufft_type2_eval_3d(txn, tyn, tzn, 1, eps, coeff)
        pref = dk[1] * dk[2] * dk[3] / (2π)^3
        for (m, i) in enumerate(bidx)
            out[i] = pref * real(vals[m])
        end
    end
    if !isempty(oidx)
        vals = lfmm3d(eps, src; charges = q, targets = trg[:, oidx], pgt = 1)
        for (m, i) in enumerate(oidx)
            out[i] = vals.pottarg[m] / (4π)
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Panel (a): convergence vs n — hybrid and pure particle sum, several tolerances
# ---------------------------------------------------------------------------
err_tkm_n  = Dict{Float64,Vector{Float64}}()
err_fmm_n  = Dict{Float64,Vector{Float64}}()
err_near_n = Dict{Float64,Vector{Float64}}()
err_far_n  = Dict{Float64,Vector{Float64}}()
for eps in eps_list
    err_tkm_n[eps]  = Float64[]
    err_fmm_n[eps]  = Float64[]
    err_near_n[eps] = Float64[]
    err_far_n[eps]  = Float64[]
end
hn_by_n = Float64[]

@info "Panel (a) — convergence vs n, c_pad = $C_PAD, $N_near near + $N_far far targets"
for n in n_list
    vs = grid_source(n)
    g  = BI.near_field_geometry(vs; c_pad = C_PAD)
    push!(hn_by_n, g.hn)

    # the target set must classify as designed at every n
    nn = count(i -> BI.in_near_region(g, targets, i), 1:N_targets)
    nn == N_near || error("n = $n: $nn targets classified near, expected $N_near " *
                          "(B_pad upper corner = $(g.hi[1]))")

    src, q = BI._volume_source_fmm_sources(vs)
    for eps in eps_list
        field = PrecomputedVolumeField(vs; tol = eps, c_pad = C_PAD, compute_grad = false)
        u_h   = volume_field_potential(field, targets)
        push!(err_tkm_n[eps],  relerr(u_h))
        push!(err_near_n[eps], relerr(u_h, near_rng))
        push!(err_far_n[eps],  relerr(u_h, far_rng))

        # pure particle sum at every target (Eq. 3.13); tau does not enter it,
        # so the four curves coincide
        u_p = lfmm3d(eps, src; charges = q, targets = targets, pgt = 1).pottarg ./ (4π)
        push!(err_fmm_n[eps], relerr(u_p))

        @info @sprintf("  n = %3d  tau = %.0e  E_hyb = %.3e (near %.3e, far %.3e)  E_part = %.3e",
                       n, eps, err_tkm_n[eps][end], err_near_n[eps][end],
                       err_far_n[eps][end], err_fmm_n[eps][end])
    end
end

# ---------------------------------------------------------------------------
# Panel (b): periodization threshold — vary eta at fixed n
# ---------------------------------------------------------------------------
err_tkm_eta = Dict{Float64,Vector{Float64}}()
for eps in eps_list
    err_tkm_eta[eps] = Float64[]
end

@info "Panel (b) — eta sweep at n = $n_for_eta"
vs_c = grid_source(n_for_eta)

# cross-check: at eta = 1 the standalone evaluator must reproduce the production path
let g = BI.near_field_geometry(vs_c; c_pad = C_PAD)
    f  = PrecomputedVolumeField(vs_c; tol = 1e-12, c_pad = C_PAD, compute_grad = false)
    u_prod = volume_field_potential(f, targets)
    u_std  = hybrid_eval_eta(vs_c, targets, 1.0; eps = 1e-12)
    d = maximum(abs.(u_prod .- u_std)) / maximum(abs.(u_prod))
    @info @sprintf("eta = 1 cross-check vs PrecomputedVolumeField: max rel diff = %.3e", d)
    d < 1e-10 || error("standalone eta evaluator disagrees with the production path ($d)")
end

for eps in eps_list
    for eta in eta_list
        u = hybrid_eval_eta(vs_c, targets, eta; eps = eps)
        push!(err_tkm_eta[eps], relerr(u))
    end
    @info @sprintf("  tau = %.0e :  E(0.5) = %.3e   E(~1) = %.3e   E(1.6) = %.3e",
                   eps, err_tkm_eta[eps][1],
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
    c_pad       = C_PAD,
    hn_by_n     = hn_by_n,
    eps_list    = collect(eps_list),
    n_list      = n_list,
    err_tkm_n   = err_tkm_n,
    err_fmm_n   = err_fmm_n,
    err_near_n  = err_near_n,
    err_far_n   = err_far_n,
    n_for_eta   = n_for_eta,
    eta_list    = eta_list,
    err_tkm_eta = err_tkm_eta,
    N_targets   = N_targets,
    n_near      = N_near,
    n_far       = N_far,
    tail_mass   = tail_mass(box_half, s_gauss),
)

datapath = joinpath(@__DIR__, "fig5_data.jls")
open(datapath, "w") do io
    serialize(io, out)
end
@info "Saved data" datapath bytes=stat(datapath).size

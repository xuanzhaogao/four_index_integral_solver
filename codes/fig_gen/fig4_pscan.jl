#=
Diagnostic for the empirical convergence rate of standard p×p GL on the
near-field panel pair from fig4_data.jl.

For each fixed d ∈ d_fixed, sweep p ∈ p_scan, compute the same relative L²
error E_std(d, p) = ‖I_p×p - I_ref‖ / ‖I_ref‖ that fig4 reports, then fit

    log10 E_std ≈ a + b · p          (on the dependable, well-resolved p's)

and compare the empirical slope b against the two textbook predictions

    b_theory    = −2 · log10(ρ),     ρ = d + √(1 + d²)
    b_halfrate  = −  log10(ρ)

(canonical Bernstein-ρ^(-2p) vs. the apparent half-rate ρ^(-p) you fit on
the d-sweep figure). Saves figs/fig4_pscan.{pdf,png}.

Uses the SAME density σ_Q and SAME kernel as fig4_data.jl, so the slope
here is directly comparable to the fig4 right panel.
=#

using LinearAlgebra
using FastGaussQuadrature
using HCubature
using CairoMakie
using LaTeXStrings
using Printf
using BoundaryIntegral
const BI = BoundaryIntegral

# ---------------------------------------------------------------------------
# Problem setup — must mirror fig4_data.jl
# ---------------------------------------------------------------------------
const L_half = 1.0
const n_P    = (0.0, 0.0, 1.0)
σ_Q(η1, η2) = exp(-η1^2 - η2^2)

@inline function kernel(x::NTuple{3,Float64}, y::NTuple{3,Float64})
    r2 = (x[1]-y[1])^2 + (x[2]-y[2])^2 + (x[3]-y[3])^2
    return 1.0 / (4π * sqrt(r2))
end

const hc_rtol = 1e-13
const hc_atol = 1e-15
const hc_max  = 50_000_000

# Diagnostic sweep parameters
const p_scan   = [4, 6, 8, 10, 12, 14]
const d_fixed  = [0.3, 1.0, 3.0]

# ---------------------------------------------------------------------------
# Run E_std(d, p) for each (d, p)
# ---------------------------------------------------------------------------
function run_pair(p::Int, d::Float64)
    ns_p, ws_p = gausslegendre(p)
    λ_p        = BI.gl_barycentric_weights(ns_p, ws_p)

    sigma_p = Matrix{Float64}(undef, p, p)
    for j in 1:p, i in 1:p
        sigma_p[i, j] = σ_Q(ns_p[i] * L_half, ns_p[j] * L_half)
    end

    targets = Vector{NTuple{3,Float64}}(undef, p*p)
    k = 0
    for j in 1:p, i in 1:p
        k += 1
        targets[k] = (ns_p[i] * L_half, ns_p[j] * L_half, 0.0)
    end
    Nt = length(targets)

    @inline function sigma_pointwise(u, v)
        rx = BI.barycentric_row(ns_p, λ_p, u)
        ry = BI.barycentric_row(ns_p, λ_p, v)
        s = 0.0
        @inbounds for j in 1:p, i in 1:p
            s += sigma_p[i, j] * rx[i] * ry[j]
        end
        return s
    end

    # p×p GL of K · P_σ over Q
    I_std = zeros(Nt)
    @inbounds for j in 1:p, i in 1:p
        η1 = ns_p[i] * L_half
        η2 = ns_p[j] * L_half
        w  = ws_p[i] * ws_p[j] * sigma_p[i, j]
        y  = (η1, η2, d)
        for t in 1:Nt
            I_std[t] += w * kernel(targets[t], y)
        end
    end
    I_std .*= L_half * L_half

    # HCubature reference on K · P_σ
    I_ref = zeros(Nt)
    for t in 1:Nt
        x = targets[t]
        f = η -> sigma_pointwise(η[1] / L_half, η[2] / L_half) *
                kernel(x, (η[1], η[2], d))
        val, _ = hcubature(f, (-L_half, -L_half), (L_half, L_half);
                           rtol = hc_rtol, atol = hc_atol, maxevals = hc_max)
        I_ref[t] = val
    end

    return norm(I_std - I_ref) / norm(I_ref)
end

E = Dict{Float64, Vector{Float64}}()
for d in d_fixed
    Ed = zeros(length(p_scan))
    for (i, p) in enumerate(p_scan)
        t0 = time()
        Ed[i] = run_pair(p, d)
        @info @sprintf("  d=%.3f  p=%2d  E=%.3e  (%.2fs)",
                       d, p, Ed[i], time() - t0)
    end
    E[d] = Ed
end

# ---------------------------------------------------------------------------
# Slope fit (least squares on log10 E_std vs p, over the well-resolved range
# where E_std is above the HCubature/round-off floor)
# ---------------------------------------------------------------------------
function fit_slope(ps, Es; floor = 1e-13)
    mask = isfinite.(Es) .& (Es .> 10 * floor)
    n = count(mask)
    n < 2 && return (NaN, NaN, mask)
    xs = Float64.(ps[mask])
    ys = log10.(Es[mask])
    x̄ = sum(xs)/n; ȳ = sum(ys)/n
    b = sum((xs .- x̄) .* (ys .- ȳ)) / sum((xs .- x̄).^2)
    a = ȳ - b*x̄
    return (a, b, mask)
end

println("\n=== Empirical vs. theoretical slope (log10 E_std vs p) ===")
println(rpad("d", 8), rpad("ρ", 10),
        rpad("b_empirical", 14), rpad("-2 log10 ρ", 14),
        rpad("-log10 ρ", 14))
slopes = Dict{Float64, Tuple{Float64,Float64,Float64}}()
for d in d_fixed
    ρ        = d + sqrt(1 + d^2)
    a, b, _  = fit_slope(p_scan, E[d])
    th_full  = -2 * log10(ρ)
    th_half  = -log10(ρ)
    slopes[d] = (b, th_full, th_half)
    @printf("  %.3f  %.4f   %+8.4f      %+8.4f      %+8.4f\n",
            d, ρ, b, th_full, th_half)
end

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig = Figure(size = (700, 500), fontsize = 18)
ax  = Axis(fig[1, 1];
           yscale = log10,
           xlabel = L"p",
           ylabel = L"\mathcal{E}_{\mathrm{near}}",
           xticks = (p_scan, string.(p_scan)))

palette = cgrad(:viridis, length(d_fixed) + 1, categorical = true)
markers = [:circle, :rect, :utriangle, :diamond]

for (i, d) in enumerate(d_fixed)
    ρ = d + sqrt(1 + d^2)
    col = palette[i]
    Ed = max.(E[d], 1e-18)

    # data
    scatterlines!(ax, p_scan, Ed;
                  color = col, marker = markers[i],
                  markersize = 11, linewidth = 2,
                  label = L"d = %$d ,\ \rho = %$(round(ρ; digits=2))")

    # theoretical reference lines: rho^{-2p} and rho^{-p}, anchored to first datapoint
    anchor_p, anchor_E = p_scan[1], Ed[1]
    full = anchor_E .* ρ .^ (-2 .* (p_scan .- anchor_p))
    half = anchor_E .* ρ .^ (-1 .* (p_scan .- anchor_p))
    lines!(ax, p_scan, full;
           color = col, linestyle = :dash, linewidth = 1.5)
    lines!(ax, p_scan, half;
           color = col, linestyle = :dot, linewidth = 1.5)
end

# proxy legend entries for the dashed vs. dotted reference lines
lines!(ax, [NaN, NaN], [NaN, NaN]; color = :black, linestyle = :dash,
       linewidth = 1.5, label = L"\propto \rho^{-2p}")
lines!(ax, [NaN, NaN], [NaN, NaN]; color = :black, linestyle = :dot,
       linewidth = 1.5, label = L"\propto \rho^{-p}")

axislegend(ax; position = :lb, labelsize = 13)
ylims!(ax, 1e-16, 1e1)

outpath = joinpath(@__DIR__, "figs/fig4_pscan.pdf")
save(outpath, fig; px_per_unit = 4)
png_out = replace(outpath, ".pdf" => ".png")
save(png_out, fig; px_per_unit = 4)
@info "Saved diagnostic" outpath png_out

fig

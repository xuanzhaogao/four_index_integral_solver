#=
fig4_plot.jl

Two-panel figure for Section "Near-field evaluation of layer potentials".

(a)  Pointwise relative quadrature error in the (x_t, z_t) plane at y_t = 0
     for the standard n x n Gauss-Legendre rule (n = p_target) applied to
     the diagnostic single-layer integrand on the panel [-1,1]^2 x {0}.
     Filled contours = actual log_10 (|I_n - I_ref| / max|I_ref|),
     overlaid line contours = af-Klinteberg prediction of the same quantity
     using Error estimate 6 (Eq. 81, Klinteberg-Sorgentone-Tornberg 2022)
     integrated along the slower of the two reference axes.

(b)  E_near vs panel-panel gap d at p = p_target = 8 for
       - standard p x p Gauss-Legendre rule (degrades for small d),
       - upsampled rule with p_up dictated by inverting the
         Bernstein-radius bound  (machine-precision floor in the
         near regime up to the p_up cap).
     Loaded from fig4_corr.jls.
=#

using Base.Threads
using Serialization
using FastGaussQuadrature
using HCubature
using SpecialFunctions
using CairoMakie
using LaTeXStrings
using Printf

# ---------------------------------------------------------------------------
# Right-panel data
# ---------------------------------------------------------------------------
const data    = open(deserialize, joinpath(@__DIR__, "fig4_corr.jls"), "r")
const d_list  = data.d_list
const E_std   = data.E_std
const E_corr  = data.E_corr
const p_targ  = data.p_target
const eps_tol = data.ε
const L_half  = data.L_half

const rho_star = eps_tol ^ (-1.0 / (2 * p_targ))

# ---------------------------------------------------------------------------
# Density and single-layer kernel (must match fig4_data.jl / fig4_corr.jl).
# ---------------------------------------------------------------------------
σ_Q(η1, η2) = exp(-η1^2 - η2^2)

@inline function sl_kernel(xt, yt, zt, y1, y2)
    r2 = (xt - y1)^2 + (yt - y2)^2 + zt^2
    return 1.0 / (4π * sqrt(r2))
end

# n_bern: order of the standard rule for the LHS panel; must match the
# scatter markers in the RHS panel.
const n_bern  = p_targ
const yt_bern = 0.0

# Plot window.
const xa, xb = -3.5, 3.5
const za, zb = -3.0, 3.0
const xt_s   = collect(range(xa, xb; length = 180))
const zt_s   = collect(range(za, zb; length = 150))

# ---------------------------------------------------------------------------
# Actual quadrature error on the (x_t, z_t) grid.
# Standard rule:  p_targ-point Gauss-Legendre (matches RHS panel).
# Reference:      adaptive cubature via HCubature, with absolute and relative
#                 tolerance well below the deepest filled contour level
#                 (10^{-15}); the integrand has only an integrable boundary
#                 singularity as z_t -> 0 with |x_t| < 1, which the adaptive
#                 subdivision resolves correctly.
# ---------------------------------------------------------------------------
function integrate_gl(xs_s, ws_s, xt, yt, zt)
    res = 0.0
    @inbounds for j in eachindex(xs_s)
        η2 = xs_s[j] * L_half
        wy = ws_s[j]
        for i in eachindex(xs_s)
            η1 = xs_s[i] * L_half
            res += σ_Q(η1, η2) * ws_s[i] * wy *
                   sl_kernel(xt, yt, zt, η1, η2)
        end
    end
    return res * L_half * L_half
end

function integrate_ref(xt, yt, zt)
    f = η -> σ_Q(η[1], η[2]) * sl_kernel(xt, yt, zt, η[1], η[2])
    val, _ = hcubature(f, (-L_half, -L_half), (L_half, L_half);
                       rtol = 1e-14, atol = 1e-16, maxevals = 50_000_000)
    return val
end

xs_n, ws_n = gausslegendre(n_bern)

res_ref = Matrix{Float64}(undef, length(xt_s), length(zt_s))
res_n   = Matrix{Float64}(undef, length(xt_s), length(zt_s))

# Flatten (i, j) so each thread can grab the next un-claimed cell.
const total_pts = length(xt_s) * length(zt_s)
@info "Computing HCubature reference on grid" total_pts threads=Threads.nthreads()
t_start = time()
@threads for k in 1:total_pts
    j, i = divrem(k - 1, length(xt_s)) .+ (1, 1)
    xt = xt_s[i]; zt = zt_s[j]
    res_ref[i, j] = integrate_ref(xt, yt_bern, zt)
    res_n[i, j]   = integrate_gl(xs_n, ws_n, xt, yt_bern, zt)
end
@info @sprintf("Reference sweep done in %.1fs", time() - t_start)

const denom = maximum(abs.(res_ref))
log_true    = log10.(abs.(res_n .- res_ref) ./ denom .+ eps())

# ---------------------------------------------------------------------------
# af-Klinteberg prediction: Error estimate 6 (Eq. 81) integrated along
# the slower of the two reference axes.
#
# For target x = (xt, yt, zt) and source panel [-L,L]^2 at z = 0:
#   - s1-axis pole, s2 fixed:  t0(s2) = (xt + i sqrt((yt - L s2)^2 + zt^2)) / L
#   - s2-axis pole, s1 fixed:  t0(s1) = (yt + i sqrt((xt - L s1)^2 + zt^2)) / L
#
# Single-layer kernel ⇒ p_k = 1/2.  Geometry factor G = 1 (flat panel).
# Smooth factor f = σ_Q / (4π) analytically continued in the polar variable.
# est(t0, n, p) = (4π / Γ(p)) · (2n+1)^{p-1} / (2 |t0^2 - 1|^p) · ρ(t0)^{-(2n+1)}.
# ---------------------------------------------------------------------------
@inline function bernstein_rho(z::Complex)
    s = sqrt(z * z - one(real(z)))
    return max(abs(z + s), abs(z - s))
end

const p_k     = 0.5
const n_quad  = 24
const xs_quad, ws_quad = gausslegendre(n_quad)
const est_pref = (4π / gamma(p_k)) * (2 * n_bern + 1)^(p_k - 1) / 2

@inline function est_at(t0::Complex)
    ρ      = bernstein_rho(t0)
    t2m1   = abs(t0 * t0 - 1)
    return est_pref / t2m1^p_k * ρ^(-(2 * n_bern + 1))
end

# 1D integrand along s2 (for the s1-axis pole) at target (xt, yt, zt).
@inline function I1_at(xt, yt, zt)
    res = 0.0
    @inbounds for k in eachindex(xs_quad)
        s2 = xs_quad[k]
        η2 = s2 * L_half
        t0 = (xt + im * sqrt((yt - η2)^2 + zt^2)) / L_half
        # σ_Q continued in the first arg: σ_Q(L * t0, η2) = exp(-(L t0)^2 - η2^2)
        fσ = abs(exp(-(L_half * t0)^2 - η2^2)) / (4π)
        res += ws_quad[k] * fσ * est_at(t0) * L_half
    end
    return res
end

@inline function I2_at(xt, yt, zt)
    res = 0.0
    @inbounds for k in eachindex(xs_quad)
        s1 = xs_quad[k]
        η1 = s1 * L_half
        t0 = (yt + im * sqrt((xt - η1)^2 + zt^2)) / L_half
        fσ = abs(exp(-η1^2 - (L_half * t0)^2)) / (4π)
        res += ws_quad[k] * fσ * est_at(t0) * L_half
    end
    return res
end

err_pred = Matrix{Float64}(undef, length(xt_s), length(zt_s))
for j in eachindex(zt_s), i in eachindex(xt_s)
    err_pred[i, j] = I1_at(xt_s[i], yt_bern, zt_s[j]) +
                     I2_at(xt_s[i], yt_bern, zt_s[j])
end
log_pred = log10.(err_pred ./ denom .+ eps())

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
begin
    fig = Figure(size = (1050, 420), fontsize = 18)

    # ----- Panel (a): actual error + Klinteberg prediction overlay -----
    ax_a = Axis(fig[1, 1];
                aspect  = DataAspect(),
                xlabel  = L"x_t",
                ylabel  = L"z_t")

    levels_log = -10.0:10/6:0.0

    hm = contourf!(ax_a, xt_s, zt_s, log_true;
                   levels   = levels_log,
                   colormap = :viridis,
                   rasterize = 4)

    contour!(ax_a, xt_s, zt_s, log_pred;
             levels    = levels_log,
             color     = :black,
             linewidth = 1.3)

    lines!(ax_a, [-L_half, L_half], [0.0, 0.0];
           color = :red, linewidth = 3)

    xlims!(ax_a, -2.3, 2.3)
    ylims!(ax_a, -1.8, 1.8)

    Colorbar(fig[1, 2], hm;
             label = L"\log_{10}\,\mathcal{E}",
             width = 12)

    # ----- Panel (b): standard vs upsampled error at p = 8 -----
    ax_b = Axis(fig[1, 3];
                xscale = log10,
                yscale = log10,
                xlabel = L"d",
                ylabel = L"\mathcal{E}_{\mathrm{near}}")

    floor_y = 1e-17
    clip(y) = max(y, floor_y)

    scatter!(ax_b, d_list, clip.(E_std);
             color = (:indigo, 0.85), marker = :circle,
             markersize = 11, label = "standard")
    scatter!(ax_b, d_list, clip.(E_corr);
             color = (:seagreen, 0.95), marker = :utriangle,
             markersize = 12, label = "upsampled")

    hlines!(ax_b, [eps_tol];
            color = :gray, linestyle = :dash, linewidth = 1.4)
    # text!(ax_b, L"\varepsilon = 10^{-12}";
    #       position = (10^(0.55), eps_tol * 2.0),
    #       fontsize = 14, color = :gray)

    axislegend(ax_b; position = :rt, framevisible = true)
    ylims!(ax_b, 1e-16, 1e0)
    xlims!(ax_b, 10^(-1.1), 10^(1.1))

    colgap!(fig.layout, 1, 6)
    colgap!(fig.layout, 2, 32)

    outpath = joinpath(@__DIR__, "figs/fig4_near_correction.pdf")
    save(outpath, fig; px_per_unit = 4)
    png_out = replace(outpath, ".pdf" => ".png")
    save(png_out, fig; px_per_unit = 4)
    @info "Saved figure" outpath png_out

    fig
end


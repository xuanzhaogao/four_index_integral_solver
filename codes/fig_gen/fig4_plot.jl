#=
Figure 4 plot script.

Left  panel: Bernstein-ellipses picture in the (x_t, z_t) plane at y_t = 0.
             Same density σ_Q and panel Q = [-1, 1]² × {0} as the right
             panel; same Laplace double-layer kernel as fig4_data.jl.
             Filled contours = log10|I_n - I_ref| for n = max(p_values).
Right panel: standard p×p GL near-field error E_std vs the absolute gap d
             for p ∈ {4, 6, 8}, loaded from fig4_data.jls.
=#

using Serialization
using FastGaussQuadrature
using CairoMakie
using LaTeXStrings
using Printf

# ---------------------------------------------------------------------------
# Right-panel data (precomputed by fig4_data.jl)
# ---------------------------------------------------------------------------
const datapath = joinpath(@__DIR__, "fig4_data.jls")
const data     = open(deserialize, datapath, "r")

const d_list    = data.d_list
const p_values  = sort(collect(data.p_values))
const sweeps    = Dict(s.p => s for s in data.sweeps)
const L_half    = data.L_half

# ---------------------------------------------------------------------------
# Shared density / kernel (must match fig4_data.jl)
# ---------------------------------------------------------------------------
σ_Q(η1, η2) = exp(-η1^2 - η2^2)

@inline function dt_kernel(xt, yt, zt, y1, y2)
    # Laplace single-layer kernel 1/(4π r) — matches fig4_data.jl.
    r2 = (xt - y1)^2 + (yt - y2)^2 + zt^2
    return 1.0 / (4π * sqrt(r2))
end

# Source integral over y ∈ [-L_half, L_half]² × {0}, evaluated at the
# target (xt, yt, zt). The n-point GL quadrature samples η = s·L_half with
# s ∈ [-1, 1] and weights ws.
function integrate_dt(xs_s, ws_s, xt, yt, zt)
    res = 0.0
    @inbounds for j in eachindex(xs_s)
        η2 = xs_s[j] * L_half
        wy = ws_s[j]
        for i in eachindex(xs_s)
            η1 = xs_s[i] * L_half
            res += σ_Q(η1, η2) * ws_s[i] * wy *
                   dt_kernel(xt, yt, zt, η1, η2)
        end
    end
    return res * L_half * L_half
end

function bernstein_rho_from_pole(z::Complex)
    s = sqrt(z*z - one(real(z)))
    return max(abs(z + s), abs(z - s))
end

# ---------------------------------------------------------------------------
# Bernstein contour computation (LHS)
# ---------------------------------------------------------------------------
const n_bern  = maximum(p_values)
const yt_bern = 0.0
const n_ref   = 64

# Larger plot window than the panel itself (panel sits in [-L_half, L_half]).
const xa, xb = -3.5, 3.5
const za, zb = -3.0, 3.0
const xt_s   = collect(range(xa, xb; length = 240))
const zt_s   = collect(range(za, zb; length = 200))

xs_ref, ws_ref = gausslegendre(n_ref)
xs_n,  ws_n    = gausslegendre(n_bern)

res_ref  = Matrix{Float64}(undef, length(xt_s), length(zt_s))
res_n    = Matrix{Float64}(undef, length(xt_s), length(zt_s))
err_pred = Matrix{Float64}(undef, length(xt_s), length(zt_s))

for j in eachindex(zt_s), i in eachindex(xt_s)
    xt = xt_s[i]; zt = zt_s[j]
    res_ref[i, j] = integrate_dt(xs_ref, ws_ref, xt, yt_bern, zt)
    res_n[i, j]   = integrate_dt(xs_n,  ws_n,  xt, yt_bern, zt)

    # Bernstein analysis in reference s ∈ [-1, 1] (s = η / L_half).
    # Pole of K(·, ·, zt) in complex η_1, with η_2 real and bounded by
    # [-L_half, L_half], is at η_1 = xt ± i·sqrt((yt - η_2)² + zt²).
    # In s_1 = η_1 / L_half: s_1^pole = (xt + i·sqrt(...)) / L_half.
    # ρ_1 = min over s_2 ∈ [-1, 1] of |s_1^pole + sqrt((s_1^pole)² - 1)|.
    rho_x_min = Inf
    rho_y_min = Inf
    for s2 in -1.0:0.1:1.0
        η2 = s2 * L_half
        s_pole_x = (xt + im * sqrt((yt_bern - η2)^2 + zt^2)) / L_half
        s_pole_y = (yt_bern + im * sqrt((xt - η2)^2 + zt^2)) / L_half
        rho_x_min = min(rho_x_min, bernstein_rho_from_pole(s_pole_x))
        rho_y_min = min(rho_y_min, bernstein_rho_from_pole(s_pole_y))
    end
    err_pred[i, j] = max(rho_x_min^(-2 * n_bern), rho_y_min^(-2 * n_bern))
end

# Use relative error so the colormap is comparable across kernels.
denom    = maximum(abs.(res_ref))
log_true = log10.(abs.(res_n .- res_ref) ./ denom .+ eps())
log_pred = log10.(err_pred .+ eps())

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
begin
    fig = Figure(size = (1000, 400), fontsize = 18)

    # Panel (a): Bernstein ellipses -------------------------------------------------
    levels_log = -15.0:2.0:1.0

    ax_a = Axis(fig[1, 1];
                aspect = DataAspect(),
                xlabel = L"x_t",
                ylabel = L"z_t",
                # title  = L"\log_{10}\,|I_{%$n_bern} - I_{\mathrm{ref}}| / \max|I_{\mathrm{ref}}|,\ y_t = 0"
                )

    hm = contourf!(ax_a, xt_s, zt_s, log_true;
                   levels = levels_log, colormap = :viridis, rasterize = 4)
    # contour!(ax_a, xt_s, zt_s, log_pred;
    #          levels = levels_log,
    #          color = :black, linewidth = 1.4)
    # Mark the integration interval [-L_half, L_half] at z = 0
    lines!(ax_a, [-L_half, L_half], [0.0, 0.0];
           color = :red, linewidth = 3)

    xlims!(ax_a, -3.5, 3.5)
    ylims!(ax_a, -3, 3)

    Colorbar(fig[1, 2], hm; label = L"\log_{10}\,\mathcal{E}", width = 12)

    # Panel (b): E_std vs d/h for p = 4, 6, 8 ---------------------------------------
    ax_b = Axis(fig[1, 3];
                xscale = log10, yscale = log10,
                xlabel = L"d",
                ylabel = L"\mathcal{E}_{\mathrm{near}}")

    floor_y = 1e-16
    clip(y) = max(y, floor_y)

    p_palette = cgrad(:viridis, length(p_values) + 1, categorical = true)
    markers   = [:circle, :rect, :utriangle]
    # Canonical Bernstein prediction: ρ^(-2p)/(ρ²-1), confirmed by the
    # fixed-d / sweep-p diagnostic in fig4_pscan.jl. The per-p constants
    # below are calibrated at d ≈ 1 (well-resolved, above the HCubature
    # floor); the d-dependence of the actual prefactor is more complex,
    # so the dashed lines drift from the data at the extremes.
    factors = [1.1, 2.0, 2.6]
    for (i, p) in enumerate(p_values)
        sw  = sweeps[p]
        col = p_palette[i]
        scatter!(ax_b, d_list, clip.(sw.E_std);
                      color = col, marker = markers[i],
                      markersize = 11, label = L"p = %$p")

        f_temp = x -> (x + sqrt(1 + x^2))^(- 2 * p) / ((x + sqrt(1 + x^2))^2 - 1) / factors[i]
        lines!(ax_b, d_list, f_temp.(d_list);
               color = col, linewidth = 1.5, linestyle = :dash)
    end
    axislegend(ax_b; position = :rt)
    ylims!(ax_b, 1e-16, 1e0)
    xlims!(ax_b, 10^(-1.1), 10^(1.1))

    text!(ax_b, L"O\left( \frac{\rho^{- 2p}}{\rho^2 - 1} \right)";
          position = (10^(-1.0), 1e-12),
          fontsize = 18, color = :black)

    colgap!(fig.layout, 1, 6)
    colgap!(fig.layout, 2, 28)

    outpath = joinpath(@__DIR__, "figs/fig4_near_correction.pdf")
    save(outpath, fig; px_per_unit = 4)
    png_out = replace(outpath, ".pdf" => ".png")
    save(png_out, fig; px_per_unit = 4)
    @info "Saved figure" outpath png_out

    fig
end

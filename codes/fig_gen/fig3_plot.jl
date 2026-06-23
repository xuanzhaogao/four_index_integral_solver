#=
Figure 3 plot script. Loads pre-computed Figure 3 artifacts and generates
figs/fig3_edge_singularity.pdf.

Layout:
  (a) 1D slices of |σ| along x on the +z face for three dielectric contrasts.
      Log y-axis makes the edge singularity at x = ±1 visible.
  (b) Charge accuracy for three dielectric contrasts. For each refinement
      level we plot |∫_Γ σ dS - (1 - 1/ε)| for an interior unit source.
=#

using Serialization
using CairoMakie
using LaTeXStrings
using Printf
using FastGaussQuadrature
using BoundaryIntegral
const BI = BoundaryIntegral

include(joinpath(@__DIR__, "fig_style.jl"))

const surfacepath = joinpath(@__DIR__, "fig3_surface_eps_sweep.jls")
const surf = open(deserialize, surfacepath, "r")

const epspath = joinpath(@__DIR__, "fig3_data_eps_sweep.jls")
const epsdata = open(deserialize, epspath, "r")

const theorypath = joinpath(@__DIR__, "fig3_theory_density.jls")
const theory = open(deserialize, theorypath, "r")
const theory_by_contrast = Dict((r.contrast_num, r.contrast_den) => r for r in theory.theory_results)

# Choose panel (b)'s x-axis without regenerating data: :N or :l_min.
const panel_b_xaxis = :l_min
@assert panel_b_xaxis in (:N, :l_min)

# ---- Resampler: σ at an arbitrary (x, y) on the +z face ---------------------
@inline function panel_frame(panel)
    a, b, c, d = panel.corners
    cc  = ((a[1]+b[1]+c[1]+d[1])/4, (a[2]+b[2]+c[2]+d[2])/4, (a[3]+b[3]+c[3]+d[3])/4)
    bma = (b[1]-a[1], b[2]-a[2], b[3]-a[3])
    dma = (d[1]-a[1], d[2]-a[2], d[3]-a[3])
    return cc, bma, dma
end

function top_face_panels(interface, Lz)
    return [(i, p) for (i, p) in enumerate(interface.panels)
            if all(abs(c[3] - Lz/2) < 1e-12 for c in p.corners)]
end

function sigma_on_top_face(interface, sigma_p::Vector{Float64}, x, y, p_quad;
                           top_panels, ns_p, λ_p)
    for (idx, panel) in top_panels
        cc, bma, dma = panel_frame(panel)
        bnorm2 = bma[1]^2 + bma[2]^2 + bma[3]^2
        dnorm2 = dma[1]^2 + dma[2]^2 + dma[3]^2
        dx = x - cc[1]; dy = y - cc[2]
        u = 2 * (dx * bma[1] + dy * bma[2]) / bnorm2
        v = 2 * (dx * dma[1] + dy * dma[2]) / dnorm2
        if abs(u) <= 1 + 1e-10 && abs(v) <= 1 + 1e-10
            rx = BI.barycentric_row(ns_p, λ_p, u)
            ry = BI.barycentric_row(ns_p, λ_p, v)
            base = (idx - 1) * p_quad^2
            s = 0.0
            @inbounds for i in 1:p_quad, j in 1:p_quad
                s += rx[i] * ry[j] * sigma_p[base + (i - 1) * p_quad + j]
            end
            return s
        end
    end
    return NaN
end

# Build log-spaced d = 1 - x profile for each contrast using saved meshes + σ.
const profile_y = 0.6
const profile_d  = 10 .^ range(-4, 0; length = 400)
let
    profile_sigma = Matrix{Float64}(undef, length(profile_d), length(surf.surface_results))
    for (ir, r) in enumerate(surf.surface_results)
        top = top_face_panels(r.interface, surf.Lz)
        ns_p = gausslegendre(surf.p_quad)[1]
        λ_p  = BI.gl_barycentric_weights(gausslegendre(surf.p_quad)...)
        for (id, dv) in enumerate(profile_d)
            profile_sigma[id, ir] = sigma_on_top_face(r.interface, r.sigma,
                                                      1.0 - dv, profile_y, surf.p_quad;
                                                      top_panels = top,
                                                      ns_p = ns_p, λ_p = λ_p)
        end
    end
    global const PROFILE_SIGMA = profile_sigma
end

begin
    fig = Figure(size = (FIG_W, FIG_H))

    # ----- Panel (a): |σ| vs distance to the x = 1 edge, log-log ----------
    # Log-spaced sampling using the deepest edge-local mesh, so the curve
    # resolves the divergence across several decades in d.
    ax_a = Axis(fig[1, 1];
                xscale = log10, yscale = log10,
                xlabel = L"d", ylabel = L"|\sigma|")

    colors_a = sweep_colors(length(surf.surface_results))
    theory_d = 10 .^ range(-3.5, -0.5; length = 80)
    theory_anchor_d = 1e-2
    anchor_idx = argmin(abs.(log.(profile_d ./ theory_anchor_d)))
    for (ir, r) in enumerate(surf.surface_results)
        line = abs.(PROFILE_SIGMA[:, ir])
        keep = isfinite.(line) .& (line .> 0)
        lines!(ax_a, profile_d[keep], line[keep];
               color = colors_a[ir], linewidth = LW_DATA,
               label = L"\epsilon = %$(round(r.eps_d; sigdigits = 3))")

        theory_r = theory_by_contrast[(r.contrast_num, r.contrast_den)]
        y_anchor = line[anchor_idx]
        theory_line = y_anchor .* (theory_d ./ profile_d[anchor_idx]) .^ theory_r.density_power
        lines!(ax_a, theory_d, theory_line;
               color = colors_a[ir], linewidth = LW_GUIDE, linestyle = :dash)
    end

    vlines!(ax_a, [1.01 / 2^l for l in 7:7]; color = :gray, linewidth = LW_GUIDE, linestyle = :dash)
    text!(ax_a, L"d = l_{\text{min}}", position = (1.2 * 1e-3, 10^(-1.45)), fontsize = FS_ANNOT)

    text!(ax_a, L"O(d^{\beta})", position = (4e-2, 10^(-1.2)),
           color = :black, fontsize = FS_ANNOT)
    axislegend(ax_a; position = :rt)
    ylims!(ax_a, 10^(-1.75), 10^(-0.5))

    # ----- Panel (b): charge-neutrality accuracy --------------------------
    ax_b = Axis(fig[1, 2];
                xscale = log10, yscale = log10,
                xlabel = panel_b_xaxis === :N ? "DOF" : L"\ell_{\min}",
                ylabel = L"\mathcal{E}_{\sigma}",
                xreversed = true
                )
                # ylabel = L"\left|\int_\Gamma \sigma\,dS - (1 - 1/\epsilon)\right|")

    eps_colors = sweep_colors(length(epsdata.epsilon_results))
    for (i, r) in enumerate(epsdata.epsilon_results)
        xs = panel_b_xaxis === :N ? [l.N for l in r.levels] : [l.l_min for l in r.levels]
        errs = [l.charge_error_abs for l in r.levels]
        keep = isfinite.(errs) .& (errs .> 0)
        scatter!(ax_b, xs[keep], errs[keep];
                 color = eps_colors[i],
                 marker = :circle, markersize = MS)

        theory_r = theory_by_contrast[(r.contrast_num, r.contrast_den)]
        slope = 1 + theory_r.density_power
        fit_xs = xs[keep]
        fit_errs = errs[keep]
        logC = sum(log.(fit_errs) .- slope .* log.(fit_xs)) / length(fit_xs)
        guide_xs = 10 .^ range(log10(minimum(fit_xs)) - 0.5, log10(maximum(fit_xs)) + 0.5; length = 80)
        guide_errs = exp(logC) .* guide_xs .^ slope
        lines!(ax_b, guide_xs, guide_errs;
               color = eps_colors[i], linewidth = LW_GUIDE, linestyle = :dash)
    end
    text!(ax_b, L"O(\ell_{\min}^{\,\beta + 1})", position = (0.015, 10^(-1.5)),
           color = :black, fontsize = FS_ANNOT)

    ylims!(ax_b, 10^(-5), 10^(-0))
    xlims!(ax_b, 10^(-0.1), 10^(-2.6))

    colgap!(fig.layout, 1, 30)

    outpath = joinpath(@__DIR__, "figs", "fig3_edge_singularity.pdf")
    save(outpath, fig; px_per_unit = PX_PER_UNIT)
    @info "Saved figure" outpath

    fig
end

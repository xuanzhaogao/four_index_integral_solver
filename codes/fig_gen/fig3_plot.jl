#=
Figure 3 plot script. Loads pre-computed data from fig3_data.jls and
generates fig3_edge_singularity.png.

Layout:
  (a) 1D slices of |σ| along x at several fixed y values on the +z face. Log
      y-axis makes the edge singularity at x = ±1 visible; slices closer to
      y = 1 ride higher overall because they are already near a second edge.
  (b) Self-convergence of the far-field single-layer potential. For each
      refinement level we evaluate u_ℓ at a set of targets on a sphere of
      radius R_far outside the cube; the deepest edge-local level serves as
      reference, and we plot ||u_ℓ − u_ref||_2 / ||u_ref||_2 vs N.
=#

using Serialization
using CairoMakie
using LaTeXStrings
using Printf
using FastGaussQuadrature
using BoundaryIntegral
const BI = BoundaryIntegral

const datapath = joinpath(@__DIR__, "fig3_data.jls")
const data = open(deserialize, datapath, "r")

const surfacepath = joinpath(@__DIR__, "fig3_surface_data.jls")
const surf = open(deserialize, surfacepath, "r")

const pscanpath = joinpath(@__DIR__, "fig3_data_pscan.jls")
const pscan = open(deserialize, pscanpath, "r")

# Optional deeper reference (p, k) = (6, 13) from fig3_data_pref.jl.
const prefpath = joinpath(@__DIR__, "fig3_data_pref.jls")
const pref = isfile(prefpath) ? open(deserialize, prefpath, "r") : nothing

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

# Build log-spaced d = 1 - x profile per y slice using the saved mesh + σ.
const profile_ys = [0.0, 0.3, 0.6, 0.9]
const profile_d  = 10 .^ range(-4, 0; length = 400)
let
    top = top_face_panels(surf.interface, surf.Lz)
    ns_p = gausslegendre(surf.p_quad)[1]
    λ_p  = BI.gl_barycentric_weights(gausslegendre(surf.p_quad)...)
    profile_sigma = Matrix{Float64}(undef, length(profile_d), length(profile_ys))
    for (jy, yv) in enumerate(profile_ys)
        for (id, dv) in enumerate(profile_d)
            profile_sigma[id, jy] = sigma_on_top_face(surf.interface, surf.sigma,
                                                     1.0 - dv, yv, surf.p_quad;
                                                     top_panels = top,
                                                     ns_p = ns_p, λ_p = λ_p)
        end
    end
    global const PROFILE_SIGMA = profile_sigma
end

begin
    fig = Figure(size = (1000, 450), fontsize = 20)

    # ----- Panel (a): |σ| vs distance to the x = 1 edge, log-log ----------
    # Log-spaced sampling using the deepest edge-local mesh, so the curve
    # resolves the divergence across several decades in d.
    ax_a = Axis(fig[1, 1];
                xscale = log10, yscale = log10,
                xlabel = L"d", ylabel = L"|\sigma|")

    colors_a = [:black, :royalblue, :seagreen, :crimson]
    for (jy, y_val) in enumerate(profile_ys)
        line = abs.(PROFILE_SIGMA[:, jy])
        keep = isfinite.(line) .& (line .> 0)
        lines!(ax_a, profile_d[keep], line[keep];
               color = colors_a[jy], linewidth = 2,
               label = L"y = %$(round(y_val; digits = 2))")
    end
    axislegend(ax_a; position = :lb)

    # ----- Panel (b): far-field potential self-convergence ----------------
    ax_b = Axis(fig[1, 2];
                xscale = log10, yscale = log10,
                xlabel = "DOF",
                ylabel = L"\mathcal{E}_u")

    # Reference u_far: prefer fig3_data_pref.jls (deeper solve, e.g. p=6,k=13).
    # Fall back to (pscan.ref_p, pscan.ref_k) from the pscan artifact itself.
    local u_ref, ref_norm, ref_label
    if pref !== nothing
        u_ref = pref.u_far
        ref_label = (p = pref.p_ref, k = pref.k_ref)
    else
        rref = first(r for r in pscan.edge_results if r.p == pscan.ref_p)
        lref = first(l for l in rref.levels if l.k == pscan.ref_k)
        u_ref = lref.u_far
        ref_label = (p = pscan.ref_p, k = pscan.ref_k)
    end
    ref_norm = sqrt(sum(abs2, u_ref))
    @info "panel (b) reference" ref_label

    p_colors = [:black, :royalblue, :seagreen, :darkorange, :crimson]
    for (i, r) in enumerate(pscan.edge_results)
        Ns   = [l.N for l in r.levels]
        errs = [sqrt(sum(abs2, l.u_far .- u_ref)) / ref_norm for l in r.levels]
        keep = isfinite.(errs) .& (errs .> 0)
        scatterlines!(ax_b, Ns[keep], errs[keep];
                      color = p_colors[i],
                      marker = :circle, markersize = 10, linewidth = 2,
                      label = L"p = %$(r.p)")
    end
    axislegend(ax_b; position = :lb)

    colgap!(fig.layout, 1, 30)

    outpath = joinpath(@__DIR__, "figs", "fig3_edge_singularity.pdf")
    save(outpath, fig)
    @info "Saved figure" outpath

    fig
end

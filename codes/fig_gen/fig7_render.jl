#=
Render the controlled stacked-box dielectric model used in the end-to-end
convergence experiment (see slab_model.md and prompts_fig7.md).

Geometry (two axis-aligned dielectric bodies sharing the z = 10 face):
    Omega_2 = [-5,5] x [-5,5] x [ 0, 10]      (bulk substrate, ε_2 = 10)
    Omega_1 = [-5,5] x [-5,5] x [10, 11]      (thin layer,     ε_1 = 4)
    Gamma_12 = [-5,5] x [-5,5] x {10}          (material-material interface)
Exterior permittivity ε_0 = 1.

Densities (smooth normalized 3D Gaussians of width s = 0.12) above the layer:
    x_s = (-0.35, 0, 11.3),   x_t = ( 0.35, 0, 11.3)

For rendering we treat the stack as two visually distinct DielectricInterface
objects (one per box). The shared z = 10 face is drawn by both interfaces;
this is purely a visualization choice. A solver-side construction would
introduce the shared face as a single material-material interface.

Output: fig7_system.png
=#

using FastGaussQuadrature
using SpecialFunctions
using CairoMakie
using LaTeXStrings
using BoundaryIntegral
const BI = BoundaryIntegral
const MColors = CairoMakie.Makie.Colors
const MakieExtMod = Base.get_extension(BoundaryIntegral, :MakieExt)

# ---- Geometry / model parameters -------------------------------------------
const Lxy_half = 5.0
const box2_z = (0.0, 10.0)    # substrate
const box1_z = (10.0, 11.0)   # thin layer

const eps_1, eps_2, eps_0 = 4.0, 10.0, 1.0

const s_gauss = 0.06           # narrowed from 0.12 to suppress geometric leak
const z_src   = 11.3           # 0.3 above the top dielectric face at z = 11
const x_s = (-0.35, 0.0, z_src)
const x_t = ( 0.35, 0.0, z_src)

# Fraction of Gaussian mass that leaks into z < 11 (the upper dielectric).
# 1D z-marginal: leak = Φ(-(z_src - 11)/s) = 0.5 * erfc((z_src - 11)/(s√2)).
const leak_frac = 0.5 * erfc((z_src - 11.0) / (s_gauss * sqrt(2.0)))
@info "Geometric leak fraction (mass below z = 11)" s_gauss leak_frac

const p_quad = 4
const l_ec   = 1.0          # coarse — readability over fidelity
const n_src  = 12
const gauss_tol = 1e-6

# ---- Generic axis-aligned dielectric box constructor -----------------------
function build_box(xrange, yrange, zrange, p_quad, l_ec, eps_in, eps_out;
                   alpha = sqrt(2.0))
    ns, ws = gausslegendre(p_quad)
    xa, xb = xrange; ya, yb = yrange; za, zb = zrange

    v = NTuple{3, Float64}[
        (xb, yb, zb),  # 1  top
        (xa, yb, zb),  # 2
        (xa, ya, zb),  # 3
        (xb, ya, zb),  # 4
        (xb, yb, za),  # 5  bottom
        (xa, yb, za),  # 6
        (xa, ya, za),  # 7
        (xb, ya, za),  # 8
    ]
    faces = [
        (v[1], v[2], v[3], v[4], ( 0.0,  0.0,  1.0)),  # +z
        (v[5], v[8], v[7], v[6], ( 0.0,  0.0, -1.0)),  # -z
        (v[8], v[5], v[1], v[4], ( 1.0,  0.0,  0.0)),  # +x
        (v[7], v[3], v[2], v[6], (-1.0,  0.0,  0.0)),  # -x
        (v[6], v[2], v[1], v[5], ( 0.0,  1.0,  0.0)),  # +y
        (v[7], v[8], v[4], v[3], ( 0.0, -1.0,  0.0)),  # -y
    ]

    panels = BI.FlatPanel{Float64, 3}[]
    for (a, b, c, d, n_hat) in faces
        append!(panels, BI.rect_panel3d_adaptive_panels(
            a, b, c, d, ns, ws, n_hat,
            (true, true, true, true),
            (true, true, true, true),
            alpha, l_ec))
    end
    return BI.DielectricInterface(panels,
                                  fill(eps_in,  length(panels)),
                                  fill(eps_out, length(panels)))
end

iface_2 = build_box((-Lxy_half, Lxy_half), (-Lxy_half, Lxy_half), box2_z,
                    p_quad, l_ec, eps_2, eps_0)
iface_1 = build_box((-Lxy_half, Lxy_half), (-Lxy_half, Lxy_half), box1_z,
                    p_quad, l_ec, eps_1, eps_0)
@info "Stacked boxes" panels_substrate = length(iface_2.panels) panels_layer = length(iface_1.panels)

vs_src = BI.GaussianVolumeSource(x_s, s_gauss, n_src, gauss_tol)
vs_tar = BI.GaussianVolumeSource(x_t, s_gauss, n_src, gauss_tol)

# ---- Render ----------------------------------------------------------------
fig = Figure(size = (900, 700))
ax  = Axis3(fig[1, 1], aspect = :data,
            xlabel = "x", ylabel = "y", zlabel = "z",
            azimuth = 0.28π, elevation = 0.18π)

# Distinct base colors hint at the ε contrast (warmer/saturated → higher ε).
MakieExtMod.viz_3d!(ax, iface_2;
    show_normals = false, show_points = false,
    highlight_edges = true, fill_highlight = true, fill_alpha = 0.18,
    base_color = :steelblue, edge_color = :darkorange)
MakieExtMod.viz_3d!(ax, iface_1;
    show_normals = false, show_points = false,
    highlight_edges = true, fill_highlight = true, fill_alpha = 0.22,
    base_color = :slategray2, edge_color = :darkorange)

# Gaussian volume clouds — distinct colors with density-modulated transparency
# (CairoMakie does not accept array `alpha`, so encode α in per-point RGBA).
function scatter_gaussian!(ax, vs; rgb, label, max_points = 4000,
                           markersize = 6, min_density = 1e-3)
    pos = vs.positions
    dens = vs.density
    dmax = maximum(abs, dens)
    n = length(dens)
    stride = max(1, ceil(Int, n / max_points))
    xs = Float32[]; ys = Float32[]; zs = Float32[]
    cs = MColors.RGBA{Float32}[]
    @inbounds for i in 1:stride:n
        d = abs(dens[i]) / dmax
        d < min_density && continue
        push!(xs, pos[1, i]); push!(ys, pos[2, i]); push!(zs, pos[3, i])
        a = 0.15f0 + 0.75f0 * Float32(d)
        push!(cs, MColors.RGBA{Float32}(rgb.r, rgb.g, rgb.b, a))
    end
    scatter!(ax, xs, ys, zs;
             color = cs, markersize = markersize,
             transparency = true, label = label)
end

src_rgb = MColors.RGB{Float32}(MColors.parse(MColors.Colorant, "royalblue"))
tar_rgb = MColors.RGB{Float32}(MColors.parse(MColors.Colorant, "crimson"))
scatter_gaussian!(ax, vs_src; rgb = src_rgb, label = L"\rho_{\mathrm{src}}")
scatter_gaussian!(ax, vs_tar; rgb = tar_rgb, label = L"\rho_{\mathrm{tar}}")

axislegend(ax; position = :rt, framevisible = false)

outpath = joinpath(@__DIR__, "fig7_system.png")
save(outpath, fig; px_per_unit = 4)
@info "Saved" outpath

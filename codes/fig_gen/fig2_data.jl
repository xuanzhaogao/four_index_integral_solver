#=
Data generation for Figure 2 (RHS-driven adaptive refinement).

Runs the adaptive and uniform sweeps and stores everything needed for the
figure into fig2_data.jls. Run this once; afterwards iterate on the figure by
re-running fig2_plot.jl only.
=#

using LinearAlgebra
using FastGaussQuadrature
using Serialization
using BoundaryIntegral
const BI = BoundaryIntegral

# ---------------------------------------------------------------------------
# Problem setup (must match fig2_plot.jl's expectations)
# ---------------------------------------------------------------------------
const Lx, Ly, Lz = 1.0, 1.0, 1.0
const eps_in, eps_out = 10.0, 1.0
const eps_src = eps_out

const p_quad   = 8
const q_check  = 3 * p_quad
const l_ec     = 10.0   # disable edge/corner pass: isolate RHS-driven refinement

const σ_src    = 0.03
const n_src    = 8
const d_src    = 0.06
const x0_src   = (0.15, -0.10, Lz / 2 + d_src)

const vs = BI.GaussianVolumeSource(x0_src, σ_src, n_src, 1e-10)

function rhs_fn(p::NTuple{3,Float64}, n::NTuple{3,Float64})
    s = 0.0
    @inbounds for i in eachindex(vs.density)
        src = (vs.positions[1, i], vs.positions[2, i], vs.positions[3, i])
        s += vs.density[i] * vs.weights[i] * BI.laplace3d_grad(src, p, n)
    end
    return -s / eps_src
end

# ---------------------------------------------------------------------------
# Per-panel helpers
# ---------------------------------------------------------------------------
@inline function panel_frame(panel)
    a, b, c, d = panel.corners
    cc  = ((a[1]+b[1]+c[1]+d[1])/4, (a[2]+b[2]+c[2]+d[2])/4, (a[3]+b[3]+c[3]+d[3])/4)
    bma = (b[1]-a[1], b[2]-a[2], b[3]-a[3])
    dma = (d[1]-a[1], d[2]-a[2], d[3]-a[3])
    Lab = sqrt(bma[1]^2 + bma[2]^2 + bma[3]^2)
    Lda = sqrt(dma[1]^2 + dma[2]^2 + dma[3]^2)
    return cc, bma, dma, Lab * Lda / 4
end

@inline function panel_point(cc, bma, dma, u, v)
    (cc[1] + bma[1]*u/2 + dma[1]*v/2,
     cc[2] + bma[2]*u/2 + dma[2]*v/2,
     cc[3] + bma[3]*u/2 + dma[3]*v/2)
end

function f_L2_on_interface(interface)
    ns, ws = gausslegendre(q_check)
    s = 0.0
    for panel in interface.panels
        cc, bma, dma, scale = panel_frame(panel)
        for i in 1:q_check, j in 1:q_check
            pt = panel_point(cc, bma, dma, ns[i], ns[j])
            fv = rhs_fn(pt, panel.normal)
            s += ws[i] * ws[j] * fv * fv * scale
        end
    end
    return sqrt(s)
end

function E_f_on_interface(interface, fnorm)
    ns_p, ws_p = gausslegendre(p_quad)
    ns_c, ws_c = gausslegendre(q_check)
    λ          = BI.gl_barycentric_weights(ns_p, ws_p)
    Rcheck     = [BI.barycentric_row(ns_p, λ, ns_c[k]) for k in 1:q_check]

    err2 = 0.0
    for panel in interface.panels
        cc, bma, dma, scale = panel_frame(panel)
        F = Matrix{Float64}(undef, p_quad, p_quad)
        for i in 1:p_quad, j in 1:p_quad
            pt = panel_point(cc, bma, dma, ns_p[i], ns_p[j])
            F[i, j] = rhs_fn(pt, panel.normal)
        end
        for kx in 1:q_check
            rx = Rcheck[kx]
            for ky in 1:q_check
                ry = Rcheck[ky]
                approx = 0.0
                @inbounds for i in 1:p_quad, j in 1:p_quad
                    approx += F[i, j] * rx[i] * ry[j]
                end
                pt = panel_point(cc, bma, dma, ns_c[kx], ns_c[ky])
                exact = rhs_fn(pt, panel.normal)
                err2 += ws_c[kx] * ws_c[ky] * (exact - approx)^2 * scale
            end
        end
    end
    return sqrt(err2) / fnorm
end

function uniform_box_interface(k::Int)
    ns_p, ws_p = gausslegendre(p_quad)
    quads, normals = BI._box3d_face_quads(Lx, Ly, Lz)
    panels = BI.FlatPanel{Float64, 3}[]
    for f in 1:6
        a, b, c, d = quads[f]
        n = normals[f]
        tpl = BI.TempPanel3D(a, b, c, d,
                             true, true, true, true,
                             true, true, true, true, n)
        subs = (k == 0) ? [tpl] : BI.divide_temp_panel3d(tpl, 2^k, 2^k)
        for s in subs
            push!(panels, BI.rect_panel3d_discretize(s.a, s.b, s.c, s.d,
                                                    ns_p, ws_p, s.normal;
                                                    is_edge = true))
        end
    end
    return BI.DielectricInterface(panels,
                                  fill(eps_in,  length(panels)),
                                  fill(eps_out, length(panels)))
end

function panel_depth_on_face(panel, Lface)
    a, b, _, d = panel.corners
    edge = max(norm(b .- a), norm(d .- a))
    return max(0, round(Int, log2(Lface / edge)))
end

# ---------------------------------------------------------------------------
# Run sweeps
# ---------------------------------------------------------------------------
@info "Computing reference ||f||_{L^2(Γ)} on a fine uniform mesh ..."
ref_interface = uniform_box_interface(4)
const fnorm = f_L2_on_interface(ref_interface)
@info "fnorm" fnorm

const total_area = 2 * (Lx*Ly + Ly*Lz + Lz*Lx)
to_atol(ε) = ε * fnorm / sqrt(total_area)

adaptive_tols = [1e-2, 1e-4, 1e-6, 1e-8, 1e-10]
adaptive_data = NamedTuple[]
adaptive_interfaces = Any[]
for ε in adaptive_tols
    atol = to_atol(ε)
    @info "Adaptive sweep" ε atol
    iface = single_dielectric_box3d_rhs_adaptive(Lx, Ly, Lz, p_quad, rhs_fn,
                                                 l_ec, atol, eps_in, eps_out)
    N  = length(iface.panels) * p_quad^2
    Ef = E_f_on_interface(iface, fnorm)
    push!(adaptive_data, (ε = ε, atol = atol, N = N, Ef = Ef,
                          npanels = length(iface.panels)))
    push!(adaptive_interfaces, iface)
    @info "  result" N Ef length(iface.panels)
end

uniform_data = NamedTuple[]
for k in 0:4
    @info "Uniform sweep" k
    iface = uniform_box_interface(k)
    N  = length(iface.panels) * p_quad^2
    Ef = E_f_on_interface(iface, fnorm)
    push!(uniform_data, (k = k, N = N, Ef = Ef,
                         npanels = length(iface.panels)))
    @info "  result" N Ef length(iface.panels)
end

# ---------------------------------------------------------------------------
# Plot-ready artifacts for the tightest adaptive run
# ---------------------------------------------------------------------------
const PLOT_IDX = length(adaptive_data)   # which ε to draw in (a)/(b)
const plot_interface = adaptive_interfaces[PLOT_IDX]

# Lightweight panel records (one per panel) for the 3D panelization plot
panel_records = NamedTuple[]
for panel in plot_interface.panels
    push!(panel_records, (
        corners = (panel.corners[1], panel.corners[2],
                   panel.corners[3], panel.corners[4]),
        normal  = panel.normal,
        depth   = panel_depth_on_face(panel, Lx),
    ))
end

# (b) error grid on the +z face — relative error |f - I_P f| / max_{Γ_top}|f|
top_panels = filter(p -> all(abs(c[3] - Lz/2) < 1e-12 for c in p.corners),
                    plot_interface.panels)

nx_grid = 400
xs_grid = collect(range(-Lx/2, Lx/2; length = nx_grid))
ys_grid = collect(range(-Ly/2, Ly/2; length = nx_grid))
err_grid = fill(NaN, nx_grid, nx_grid)
f_grid   = fill(NaN, nx_grid, nx_grid)

ns_p, ws_p = gausslegendre(p_quad)
λ_p = BI.gl_barycentric_weights(ns_p, ws_p)

for panel in top_panels
    cc, bma, dma, _ = panel_frame(panel)
    F = Matrix{Float64}(undef, p_quad, p_quad)
    for i in 1:p_quad, j in 1:p_quad
        pt = panel_point(cc, bma, dma, ns_p[i], ns_p[j])
        F[i, j] = rhs_fn(pt, panel.normal)
    end

    xs_p = [c[1] for c in panel.corners]
    ys_p = [c[2] for c in panel.corners]
    x_lo, x_hi = extrema(xs_p)
    y_lo, y_hi = extrema(ys_p)
    ix_lo = searchsortedfirst(xs_grid, x_lo - 1e-12)
    ix_hi = searchsortedlast(xs_grid,  x_hi + 1e-12)
    iy_lo = searchsortedfirst(ys_grid, y_lo - 1e-12)
    iy_hi = searchsortedlast(ys_grid,  y_hi + 1e-12)

    bnorm2 = bma[1]^2 + bma[2]^2 + bma[3]^2
    dnorm2 = dma[1]^2 + dma[2]^2 + dma[3]^2

    for ig in ix_lo:ix_hi, jg in iy_lo:iy_hi
        x = xs_grid[ig]; y = ys_grid[jg]
        dx = x - cc[1]; dy = y - cc[2]
        u_local = 2 * (dx * bma[1] + dy * bma[2]) / bnorm2
        v_local = 2 * (dx * dma[1] + dy * dma[2]) / dnorm2
        (abs(u_local) > 1 + 1e-10 || abs(v_local) > 1 + 1e-10) && continue

        rx = BI.barycentric_row(ns_p, λ_p, u_local)
        ry = BI.barycentric_row(ns_p, λ_p, v_local)
        approx = 0.0
        @inbounds for i in 1:p_quad, j in 1:p_quad
            approx += F[i, j] * rx[i] * ry[j]
        end
        exact = rhs_fn((x, y, Lz/2), panel.normal)
        err_grid[ig, jg] = abs(exact - approx)
        f_grid[ig, jg]   = exact
    end
end

f_max_face = maximum(abs.(filter(isfinite, vec(f_grid))))
rel_err_grid = err_grid ./ f_max_face

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out = (
    geometry        = (Lx = Lx, Ly = Ly, Lz = Lz),
    source          = (x0 = x0_src, σ = σ_src),
    p_quad          = p_quad,
    fnorm           = fnorm,
    adaptive        = adaptive_data,
    uniform         = uniform_data,
    plot_ε          = adaptive_data[PLOT_IDX].ε,
    plot_atol       = adaptive_data[PLOT_IDX].atol,
    panel_records   = panel_records,
    xs_grid         = xs_grid,
    ys_grid         = ys_grid,
    rel_err_grid    = rel_err_grid,
    f_max_face      = f_max_face,
)

datapath = joinpath(@__DIR__, "fig2_data.jls")
open(datapath, "w") do io
    serialize(io, out)
end
@info "Saved data" datapath bytes=stat(datapath).size

println("\n=== Adaptive ===")
println("  ε         atol       N      panels    E_f")
for r in adaptive_data
    println("  $(r.ε)   $(round(r.atol; sigdigits=3))   $(r.N)   $(r.npanels)   $(round(r.Ef; sigdigits=4))")
end
println("\n=== Uniform ===")
println("  k    N      panels    E_f")
for r in uniform_data
    println("  $(r.k)   $(r.N)   $(r.npanels)   $(round(r.Ef; sigdigits=4))")
end

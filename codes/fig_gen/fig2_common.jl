#=
fig2_common.jl — shared problem setup and helpers for Figure 2.

Extracted from fig2_data.jl so that the incremental extension script
(fig2_data_uniform_extend.jl) computes E_f with exactly the same code as the
original sweep; there must be only one definition of the RHS and of the
interpolation error in play.

Note on the error definition: `E_f_abs_on_interface` returns the *absolute*
global interpolation error

    E_f = ( Σ_P ‖f - Π_P f‖²_{L²(P)} )^{1/2},

which is Eq. (Ef) of the article. Divide by `fnorm = ‖f‖_{L²(Γ)}` to recover
the relative quantity the original fig2_data.jl stored.
=#

using LinearAlgebra
using FastGaussQuadrature
using BoundaryIntegral
const BI = BoundaryIntegral

# ---------------------------------------------------------------------------
# Problem setup
# ---------------------------------------------------------------------------
const Lx, Ly, Lz = 1.0, 1.0, 1.0
const eps_in, eps_out = 10.0, 1.0
const eps_src = eps_out

const l_ec     = 10.0   # disable edge/corner pass: isolate RHS-driven refinement

const σ_src    = 0.03
const d_src    = 0.06
const x0_src   = (0.15, -0.10, Lz / 2 + d_src)

function rhs_fn(p::NTuple{3,Float64}, n::NTuple{3,Float64})
    g = BI.gaussian_laplace3d_grad(x0_src, p, σ_src)
    return -(g[1]*n[1] + g[2]*n[2] + g[3]*n[3]) / eps_src
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

function f_L2_on_interface(interface, q_check::Int)
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

# Absolute global interpolation error, Eq. (Ef) of the article.
function E_f_abs_on_interface(interface, p_quad::Int, q_check::Int)
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
    return sqrt(err2)
end

function uniform_box_interface(k::Int, p_quad::Int)
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

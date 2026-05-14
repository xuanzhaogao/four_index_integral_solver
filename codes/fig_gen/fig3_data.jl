#=
Data generation for Figure 3 (edge-singularity and edge-local refinement).

Single dielectric cube [-1,1]^3, ε_d = 10, ε_0 = 1, driven by a distant point
source at (0,0,4) so that the incident Neumann data is smooth on Γ.

For each refinement strategy (edge-local vs uniform) we solve the BIE at order
p, then compute the relative BIE residual R_BIE on an oversampled q×q
Gauss-Legendre rule (q = 3p) per leaf panel, by rebuilding the same panel
layout at order q and interpolating σ via tensor-product barycentric.

Also stored: σ on the +z face for the deepest edge-local mesh (panel a) and
edge-normal line profiles σ(r) at increasing edge-refinement depth (panel b).
=#

using LinearAlgebra
using FastGaussQuadrature
using Serialization
using BoundaryIntegral
const BI = BoundaryIntegral

# ---------------------------------------------------------------------------
# Problem setup
# ---------------------------------------------------------------------------
const Lx, Ly, Lz = 2.0, 2.0, 2.0           # cube [-1,1]^3
const eps_d, eps_0 = 10.0, 1.0
const eps_src = eps_0

const p_quad  = 4
const q_check = 3 * p_quad

const ps = BI.PointSource((0.5, 0.6, 100.0), 1.0e4)

const fmm_tol  = 1e-6
const gmres_tol = 1e-6

# Sweeps
const edge_lec_list = [1.01 / 2.0^k for k in 1:10]
const uniform_k_list = [0, 1, 2, 3, 4]

# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------
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
                                  fill(eps_d, length(panels)),
                                  fill(eps_0, length(panels)))
end

function oversample_interface(interface_p, q::Int)
    ns_q, ws_q = gausslegendre(q)
    panels_q = BI.FlatPanel{Float64, 3}[]
    for panel in interface_p.panels
        a, b, c, d = panel.corners
        push!(panels_q,
              BI.rect_panel3d_discretize(a, b, c, d, ns_q, ws_q, panel.normal;
                                         is_edge = panel.is_edge))
    end
    return BI.DielectricInterface(panels_q,
                                  copy(interface_p.eps_in),
                                  copy(interface_p.eps_out))
end

# ---------------------------------------------------------------------------
# Interpolation σ at order p → order q on the same panel layout
# ---------------------------------------------------------------------------
function build_interp_rows(p::Int, q::Int)
    ns_p, ws_p = gausslegendre(p)
    ns_q, _    = gausslegendre(q)
    λ_p = BI.gl_barycentric_weights(ns_p, ws_p)
    return [BI.barycentric_row(ns_p, λ_p, ns_q[k]) for k in 1:q]
end

function interp_sigma(sigma_p::Vector{Float64}, n_panels::Int,
                     p::Int, q::Int, Rq::Vector{<:AbstractVector{Float64}})
    sigma_q = Vector{Float64}(undef, n_panels * q^2)
    off_p, off_q = 0, 0
    for _ in 1:n_panels
        @inbounds for I in 1:q, J in 1:q
            rxI = Rq[I]
            ryJ = Rq[J]
            s = 0.0
            for i in 1:p
                ri = rxI[i]
                base = off_p + (i - 1) * p
                for j in 1:p
                    s += ri * ryJ[j] * sigma_p[base + j]
                end
            end
            sigma_q[off_q + (I - 1) * q + J] = s
        end
        off_p += p^2
        off_q += q^2
    end
    return sigma_q
end

# ---------------------------------------------------------------------------
# Solve + residual
# ---------------------------------------------------------------------------
function solve_and_residual(interface_p)
    Lhs_p = lhs_dielectric_box3d_fmm3d(interface_p, fmm_tol)
    rhs_p = rhs_dielectric_box3d(interface_p, ps, eps_src)
    sigma_p = BI.solve_gmres(Lhs_p, rhs_p, gmres_tol, gmres_tol)
    gmres_res = norm(Lhs_p * sigma_p - rhs_p) / norm(rhs_p)

    Rq = build_interp_rows(p_quad, q_check)
    interface_q = oversample_interface(interface_p, q_check)
    sigma_q = interp_sigma(sigma_p, length(interface_p.panels), p_quad, q_check, Rq)

    Lhs_q = lhs_dielectric_box3d_fmm3d(interface_q, fmm_tol)
    rhs_q = rhs_dielectric_box3d(interface_q, ps, eps_src)
    residual = Lhs_q * sigma_q - rhs_q
    ws_q = BI.all_weights(interface_q)
    r_l2 = sqrt(sum(@. ws_q * residual^2))
    f_l2 = sqrt(sum(@. ws_q * rhs_q^2))
    return sigma_p, r_l2 / f_l2, gmres_res
end

# ---------------------------------------------------------------------------
# Top-face evaluation: barycentric interpolation of σ at arbitrary (x, y, Lz/2)
# ---------------------------------------------------------------------------
@inline function panel_frame(panel)
    a, b, c, d = panel.corners
    cc  = ((a[1]+b[1]+c[1]+d[1])/4, (a[2]+b[2]+c[2]+d[2])/4, (a[3]+b[3]+c[3]+d[3])/4)
    bma = (b[1]-a[1], b[2]-a[2], b[3]-a[3])
    dma = (d[1]-a[1], d[2]-a[2], d[3]-a[3])
    return cc, bma, dma
end

function top_face_panels(interface)
    return [(i, p) for (i, p) in enumerate(interface.panels)
            if all(abs(c[3] - Lz/2) < 1e-12 for c in p.corners)]
end

"""
Evaluate σ at a query point (x, y, Lz/2) by locating the host panel on the
+z face and applying tensor-product barycentric interpolation. Returns NaN if
no host panel is found (point outside the face).
"""
function sigma_on_top_face(interface, sigma_p::Vector{Float64}, x, y;
                           top_panels = top_face_panels(interface),
                           ns_p = gausslegendre(p_quad)[1],
                           λ_p  = BI.gl_barycentric_weights(gausslegendre(p_quad)...))
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

# ---------------------------------------------------------------------------
# Run sweeps
# ---------------------------------------------------------------------------
edge_results = NamedTuple[]
edge_artifacts = Vector{Tuple{Any, Vector{Float64}}}()   # (interface, sigma)
for lec in edge_lec_list
    @info "Edge-local sweep" lec
    iface = single_dielectric_box3d(Lx, Ly, Lz, p_quad, lec, eps_d, eps_0)
    npan  = length(iface.panels)
    N     = npan * p_quad^2
    @info "  panels/N" npan N
    sigma, RBIE, gres = solve_and_residual(iface)
    push!(edge_results, (l_ec = lec, npanels = npan, N = N,
                         R_BIE = RBIE, gmres_res = gres))
    push!(edge_artifacts, (iface, sigma))
    @info "  result" R_BIE = RBIE gmres_res = gres
end

uniform_results = NamedTuple[]
uniform_artifacts = Vector{Tuple{Any, Vector{Float64}}}()
for k in uniform_k_list
    @info "Uniform sweep" k
    iface = uniform_box_interface(k)
    npan  = length(iface.panels)
    N     = npan * p_quad^2
    @info "  panels/N" npan N
    sigma, RBIE, gres = solve_and_residual(iface)
    push!(uniform_results, (k = k, npanels = npan, N = N,
                            R_BIE = RBIE, gmres_res = gres))
    push!(uniform_artifacts, (iface, sigma))
    @info "  result" R_BIE = RBIE gmres_res = gres
end

# ---------------------------------------------------------------------------
# Panel (a): σ on a dense (x, y) grid of the +z face for a 3D surface plot.
# Use a moderate edge-refinement level so the visualization grid is finer than
# the panels (avoids polynomial-extrapolation overshoot at tiny edge panels).
# ---------------------------------------------------------------------------
const PLOT_LEVEL = 3      # index into edge_artifacts (l_ec = 0.25, 312 panels)
iface_a, sigma_a = edge_artifacts[PLOT_LEVEL]
top_a = top_face_panels(iface_a)
ns_p_a = gausslegendre(p_quad)[1]
λ_p_a  = BI.gl_barycentric_weights(gausslegendre(p_quad)...)

nx_grid = 201
xs_a = collect(range(-Lx/2, Lx/2; length = nx_grid))
ys_a = collect(range(-Ly/2, Ly/2; length = nx_grid))
sigma_grid = fill(NaN, nx_grid, nx_grid)
for i in 1:nx_grid, j in 1:nx_grid
    sigma_grid[i, j] = sigma_on_top_face(iface_a, sigma_a, xs_a[i], ys_a[j];
                                         top_panels = top_a,
                                         ns_p = ns_p_a, λ_p = λ_p_a)
end

# ---------------------------------------------------------------------------
# Panel (c): far-field single-layer potential at exterior test points.
# For each refinement level we evaluate
#     u_ℓ(x*) = Σ_j G(y_j, x*) σ_ℓ(y_j) w_j
# at a set of targets on a sphere of radius R_far around the cube. Use the
# deepest level as the reference; report the relative L² error over the targets.
# ---------------------------------------------------------------------------
const R_far = 5.0
const n_theta = 12
const n_phi   = 24
far_targets = Matrix{Float64}(undef, 3, n_theta * n_phi)
let k = 1
    for i in 1:n_theta
        θ = π * (i - 0.5) / n_theta
        for j in 1:n_phi
            φ = 2π * (j - 1) / n_phi
            far_targets[1, k] = R_far * sin(θ) * cos(φ)
            far_targets[2, k] = R_far * sin(θ) * sin(φ)
            far_targets[3, k] = R_far * cos(θ)
            k += 1
        end
    end
end

function far_potential(interface, sigma::Vector{Float64})
    S = BI.laplace3d_pottrg(interface, far_targets)
    return S * sigma
end

# Evaluate at every edge-local and uniform refinement level
edge_potentials = [far_potential(edge_artifacts[k][1], edge_artifacts[k][2])
                   for k in 1:length(edge_artifacts)]
u_ref = edge_potentials[end]

edge_pot_errors = [norm(edge_potentials[k] .- u_ref) / norm(u_ref)
                   for k in 1:length(edge_potentials) - 1]
push!(edge_pot_errors, NaN)   # ref level has no error w.r.t. itself

uniform_potentials = [far_potential(uniform_artifacts[k][1], uniform_artifacts[k][2])
                      for k in 1:length(uniform_artifacts)]
uniform_pot_errors = [norm(uniform_potentials[k] .- u_ref) / norm(u_ref)
                      for k in 1:length(uniform_potentials)]

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out = (
    geometry = (Lx = Lx, Ly = Ly, Lz = Lz),
    materials = (eps_d = eps_d, eps_0 = eps_0),
    source = (point = ps.point, charge = ps.charge),
    p_quad = p_quad,
    q_check = q_check,
    edge_results = edge_results,
    uniform_results = uniform_results,
    edge_lec_list = edge_lec_list,
    surface_xs = xs_a,
    surface_ys = ys_a,
    surface_sigma = sigma_grid,
    far_R = R_far,
    far_n_targets = size(far_targets, 2),
    edge_pot_errors = edge_pot_errors,
    uniform_pot_errors = uniform_pot_errors,
)

datapath = joinpath(@__DIR__, "fig3_data.jls")
open(datapath, "w") do io
    serialize(io, out)
end
@info "Saved data" datapath bytes=stat(datapath).size

println("\n=== Edge-local ===")
println("  l_ec     panels    N        R_BIE        gmres_res")
for r in edge_results
    println("  $(r.l_ec)   $(r.npanels)   $(r.N)   $(round(r.R_BIE; sigdigits=4))   $(round(r.gmres_res; sigdigits=3))")
end
println("\n=== Uniform ===")
println("  k    panels    N        R_BIE        gmres_res")
for r in uniform_results
    println("  $(r.k)   $(r.npanels)   $(r.N)   $(round(r.R_BIE; sigdigits=4))   $(round(r.gmres_res; sigdigits=3))")
end

# 2D BIE on a square dielectric: GGQ on corner panels with HCubature-based
# correction of the panel-to-panel near-interactions at corners.
#
# Density model:        σ(s) = s^{α_1} u(s),   u smooth on [0, L_corner]
# Lagrange in u:        u(s) ≈ Σ_j ℓ_j(s) u_j      (degree m-1 poly through GGQ nodes)
# Therefore             σ(s) ≈ Σ_j (s/s_j)^{α_1} ℓ_j(s) σ_j
# Corrected entry:      A[k,j] = ∫_0^L K(s,t_k) · (s/s_j)^{α_1} · ℓ_j(s) ds
#                       computed by HCubature for any target t_k that is on a
#                       panel sharing the same corner-vertex with the source
#                       panel (the perpendicular-pair near interaction).
# Other entries:        unchanged from the standard quadrature with σ-unknowns,
#                       weights ω_j = L · w_j / s_j^{α_1}.

include(joinpath(@__DIR__, "..", "src", "gram.jl"))
include(joinpath(@__DIR__, "ggq_corner_1d.jl"))

import BoundaryIntegral as BI
using LinearAlgebra, Printf
using HCubature

const FP2 = BI.FlatPanel{Float64, 2}

# ---------- Identify corner panels ------------------------------------------

struct CornerPanelInfo
    panel_idx::Int
    vertex::NTuple{2,Float64}
    endpoint::Symbol
    L::Float64
    a::NTuple{2,Float64}
    b::NTuple{2,Float64}
    dir::NTuple{2,Float64}
    s_nodes::Vector{Float64}
    x_nodes::Vector{Float64}
    w_nodes::Vector{Float64}
    # Müntz Lagrange data.
    γshift::Vector{Float64}   # basis exponents γ_k - γ_1, length n_basis (= m)
    C::Matrix{Float64}        # orthonormal coeff: Q_k(y) = Σ_l C[k,l] y^{γshift[l]}
    Minv::Matrix{Float64}     # inverse of [Q_k(x_j)]^T so that ψ_u_j = Σ_k Minv[k,j] Q_k
end

function find_corner_panels(interface::BI.DielectricInterface{FP2,Float64},
                            vertices, x01, w01,
                            γs::AbstractVector{<:Real}, α1::Float64;
                            tol::Float64 = 1e-10)
    m = length(x01)
    γshift = collect(Float64, γs[1:m] .- γs[1])

    # Orthonormal Müntz basis under weight y^{α1} on [0,1] using exactly the m
    # exponents needed to span the Lagrange interpolant.
    Cbig, _, _ = generalized_power_basis(BigFloat(α1), BigFloat.(γshift))
    Cf = Float64.(Cbig)

    # Q_k(x_j) matrix
    Qmat = Float64.(eval_Q_matrix(BigFloat.(x01), BigFloat.(γshift), Cbig))
    # ψ_u_j(y) = Σ_k Minv[k, j] · Q_k(y);  defined by Σ_k Minv[k,j] Q_k(x_i) = δ_ij
    Minv = inv(Qmat)   # since (Qmat × Minv)[k, i] = Σ_j Q_k(x_j) Minv[j, i] -> wait order

    # Solve for Minv s.t. Σ_k Q_k(x_i) · Minv[k, j] = δ_ij  →  Qmat^T * Minv = I  →  Minv = (Qmat^T)^{-1}
    Minv = inv(Matrix(Qmat'))

    infos = CornerPanelInfo[]
    for (k, p) in enumerate(interface.panels)
        for v in vertices
            d_a = norm(p.corners[1] .- v)
            d_b = norm(p.corners[2] .- v)
            if d_a < tol || d_b < tol
                ep = d_a < tol ? :a : :b
                a = p.corners[1]; b = p.corners[2]
                L = norm(b .- a)
                far = ep == :a ? b : a
                dir = ((far[1] - v[1])/L, (far[2] - v[2])/L)
                s_nodes = L .* x01
                push!(infos, CornerPanelInfo(k, v, ep, L, a, b, dir,
                                             s_nodes, copy(x01), copy(w01),
                                             γshift, Cf, copy(Minv)))
                break
            end
        end
    end
    return infos
end

# ---------- Build interface with GGQ nodes/weights on corner panels ---------

function rebuild_panel_with_ggq(panel::FP2, vertex, endpoint, x01, w01, α1)
    a = panel.corners[1]; b = panel.corners[2]
    L = norm(b .- a)
    far = endpoint == :a ? b : a
    dir = ((far[1] - vertex[1]) / L, (far[2] - vertex[2]) / L)
    s_vals = L .* x01
    points  = [(vertex[1] + s*dir[1], vertex[2] + s*dir[2]) for s in s_vals]
    weights = [L * w01[j] / x01[j]^α1 for j in 1:length(x01)]
    n_quad  = length(x01)
    gl_xs   = collect(2 .* x01 .- 1)
    gl_ws   = copy(w01)
    return BI.FlatPanel(panel.normal, panel.corners, panel.is_edge,
                        n_quad, gl_xs, gl_ws, points, weights)
end

function apply_ggq_to_corners(interface, vertices, x01, w01, α1)
    new_panels = FP2[]
    for p in interface.panels
        replaced = false
        for v in vertices
            if norm(p.corners[1] .- v) < 1e-10
                push!(new_panels, rebuild_panel_with_ggq(p, v, :a, x01, w01, α1))
                replaced = true; break
            elseif norm(p.corners[2] .- v) < 1e-10
                push!(new_panels, rebuild_panel_with_ggq(p, v, :b, x01, w01, α1))
                replaced = true; break
            end
        end
        replaced || push!(new_panels, p)
    end
    return BI.DielectricInterface(new_panels,
                                  copy(interface.eps_in),
                                  copy(interface.eps_out))
end

# ---------- Barycentric Lagrange basis at GGQ nodes -------------------------

function bary_weights(x_nodes::Vector{Float64})
    m = length(x_nodes)
    w = ones(Float64, m)
    for j in 1:m, k in 1:m
        if k != j
            w[j] /= (x_nodes[j] - x_nodes[k])
        end
    end
    return w
end

# Evaluate ℓ_j(x) for x ∈ [0,1] using barycentric form.
function lagrange_eval(j::Int, x::Float64,
                       x_nodes::Vector{Float64}, bw::Vector{Float64};
                       tol::Float64 = 1e-14)
    m = length(x_nodes)
    # exact-node case
    for k in 1:m
        if abs(x - x_nodes[k]) < tol
            return k == j ? 1.0 : 0.0
        end
    end
    num = 0.0; den = 0.0
    for k in 1:m
        diff = x - x_nodes[k]
        c = bw[k] / diff
        den += c
        if k == j
            num = c
        end
    end
    return num / den
end

# ---------- HCubature correction --------------------------------------------

# Apply DT kernel with explicit handling for s -> 0 (s^α factor).
function dt_kernel(src::NTuple{2,Float64}, trg::NTuple{2,Float64},
                   src_normal::NTuple{2,Float64})
    dx = trg[1] - src[1]
    dy = trg[2] - src[2]
    r2 = dx*dx + dy*dy
    return (src_normal[1]*dx + src_normal[2]*dy) / (2π * r2)
end

# Müntz Lagrange basis ψ_u_j(y) = Σ_k Minv[k,j] · Q_k(y),  Q_k(y) = Σ_l C[k,l] y^{γ_l}.
function muntz_lagrange_eval(j::Int, y::Float64, src::CornerPanelInfo)
    m = length(src.x_nodes)
    val = 0.0
    @inbounds for k in 1:m
        # Q_k(y) = Σ_l C[k, l] y^{γshift[l]}
        Qk = 0.0
        for l in 1:m
            γ = src.γshift[l]
            Qk += src.C[k, l] * (γ == 0.0 ? 1.0 : y^γ)
        end
        val += src.Minv[k, j] * Qk
    end
    return val
end

"""
    corrected_entry(target, src_panel, src, j, α1; rtol)

Compute  A[k,j] = ∫_0^{L} K(s, target) · (s/s_j)^{α1} · ψ_u_j(s/L) ds
with ψ_u_j the Müntz Lagrange basis at GGQ nodes.
"""
function corrected_entry(target::NTuple{2,Float64},
                         src_panel::FP2,
                         src::CornerPanelInfo,
                         j::Int, α1::Float64;
                         rtol::Float64 = 1e-9, atol::Float64 = 1e-14,
                         maxevals::Int = 20_000)
    L  = src.L
    nrm = src_panel.normal

    integrand = function (yv)
        y = yv[1]
        y == 0.0 && return 0.0
        s = L * y
        src_pt = (src.vertex[1] + s * src.dir[1],
                  src.vertex[2] + s * src.dir[2])
        K = dt_kernel(src_pt, target, nrm)
        ratio_pow = (y / src.x_nodes[j])^α1
        ψ = muntz_lagrange_eval(j, y, src)
        return K * ratio_pow * ψ * L
    end

    val, _ = hquadrature(integrand, 0.0, 1.0;
                         rtol = rtol, atol = atol, maxevals = maxevals)
    return val
end

# ---------- Apply corrections -----------------------------------------------

"""
For each corner panel `src_info`, replace the LHS rows for targets located on
the *other* panel meeting the same vertex (the perpendicular near-interaction).
Returns the corrected LHS matrix.
"""
function apply_hcub_corrections!(LHS::Matrix{Float64},
                                 interface::BI.DielectricInterface{FP2,Float64},
                                 corner_infos::Vector{CornerPanelInfo},
                                 α1::Float64;
                                 rtol::Float64 = 1e-9,
                                 near_factor::Float64 = 4.0)
    # global offsets for each panel
    offsets = zeros(Int, length(interface.panels))
    off = 0
    for (i, p) in enumerate(interface.panels)
        offsets[i] = off
        off += BI.num_points(p)
    end

    # Group corner panels by their vertex (so we know which pairs share a vertex)
    by_vertex = Dict{NTuple{2,Float64}, Vector{Int}}()
    for (idx, ci) in enumerate(corner_infos)
        v = ci.vertex
        push!(get!(by_vertex, v, Int[]), idx)
    end

    for src in corner_infos
        src_panel = interface.panels[src.panel_idx]
        col_offset = offsets[src.panel_idx]
        m = length(src.x_nodes)
        threshold = near_factor * src.L

        # Correct entries A[k, j_src] for all targets t_k whose distance to
        # the source corner-vertex is < near_factor × L_corner. This covers
        # both the perpendicular partner panels (target right at the vertex)
        # and any nearby panels along either edge whose nodes still lie in
        # the high-mode regime of K(s, ·).
        for (tgt_panel_idx, tgt_panel) in enumerate(interface.panels)
            for k in 1:BI.num_points(tgt_panel)
                target = tgt_panel.points[k]
                d = norm(target .- src.vertex)
                d > threshold && continue
                row_offset = offsets[tgt_panel_idx]
                for j in 1:m
                    if tgt_panel_idx == src.panel_idx && k == j
                        continue
                    end
                    LHS[row_offset + k, col_offset + j] =
                        corrected_entry(target, src_panel, src, j, α1; rtol = rtol)
                end
            end
        end
    end
    return LHS
end

# ---------- Solve helpers ----------------------------------------------------

function solve_box(box, ps, eps_in, exact_flux)
    lhs = BI.lhs_dielectric_box2d(box)
    rhs = BI.rhs_dielectric_box2d(box, ps, eps_in)
    x   = BI.solve_lu(lhs, rhs)
    return x, abs(dot(BI.all_weights(box), x) - exact_flux)
end

function solve_box_corrected(box, corner_infos, α1, ps, eps_in, exact_flux;
                              rtol, near_factor::Float64 = 4.0)
    lhs = BI.lhs_dielectric_box2d(box)
    apply_hcub_corrections!(lhs, box, corner_infos, α1;
                            rtol = rtol, near_factor = near_factor)
    rhs = BI.rhs_dielectric_box2d(box, ps, eps_in)
    x   = BI.solve_lu(lhs, rhs)
    return x, abs(dot(BI.all_weights(box), x) - exact_flux)
end

# ---------- Driver -----------------------------------------------------------

function main_2d_corrected()
    eps_in  = 1.0
    eps_out = 200.0
    ps      = BI.PointSource((0.33, 0.44), 1.0)
    exact_flux = 1.0/eps_out - 1.0/eps_in

    γs = corner_power_series([pi/2], [eps_in, eps_out];
                             n_powers = 20, gmax = 30.0)
    α1 = γs[1] - 1
    @printf("γ_1 = %.6f, α_1 = %.6f\n\n", γs[1], α1)

    Lx = Ly = 1.0
    vertices = [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]

    println("Reference: very-fine GL (n_quad=20, l_corner=0.001)")
    box_ref = BI.single_dielectric_box2d(Lx, Ly, 20, 0.05, 0.001,
                                         eps_in, eps_out, Float64)
    _, err_ref = solve_box(box_ref, ps, eps_in, exact_flux)
    @printf("  N=%d  flux_err=%.3e\n\n", BI.num_points(box_ref), err_ref)

    cases = [
        (8,  0.20, 0.05),
        (12, 0.20, 0.05),
    ]
    # near_factor sweep at p=10 to gauge sensitivity
    nf_sweep = [1.0, 2.0, 4.0, 8.0]

    @printf("%-3s %-6s %-7s %-7s | %-12s %-14s %-14s\n",
            "p", "lp", "lc", "nfac", "GL err", "GGQ-raw err", "GGQ+HCub err")
    println("-"^70)
    for (nq, lp, lc) in cases
        # baseline GL
        box_gl = BI.single_dielectric_box2d(Lx, Ly, nq, lp, lc,
                                             eps_in, eps_out, Float64)
        _, err_gl = solve_box(box_gl, ps, eps_in, exact_flux)

        m = nq
        if 2m > length(γs)
            γs = corner_power_series([pi/2], [eps_in, eps_out];
                                     n_powers = 2m + 4, gmax = 60.0)
        end
        local x01, w01
        try
            x01, w01, _, _ = build_ggq(γs, m; nsteps = 40)
        catch e
            @printf("p=%d  GGQ build failed: %s\n", nq, sprint(showerror, e))
            continue
        end

        # Pure GGQ replacement (no HCubature correction)
        box_ggq = apply_ggq_to_corners(box_gl, vertices, x01, w01, α1)
        _, err_raw = solve_box(box_ggq, ps, eps_in, exact_flux)

        corner_infos = find_corner_panels(box_ggq, vertices, x01, w01, γs, α1)
        # near-factor sweep
        for nf in nf_sweep
            _, err_corr = solve_box_corrected(box_ggq, corner_infos, α1,
                                              ps, eps_in, exact_flux;
                                              rtol = 1e-10, near_factor = nf)
            @printf("%-3d %-6.3f %-7.4f nf=%-4.1f | %-12.4e %-14.4e %-14.4e\n",
                    nq, lp, lc, nf, err_gl, err_raw, err_corr)
            flush(stdout)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_2d_corrected()
end

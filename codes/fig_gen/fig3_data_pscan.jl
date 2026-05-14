#=
Multi-p sweep for fig3 panel (b).

For each polynomial order p ∈ {2,3,4,5,6} and each edge-refinement level
k = 1..10 (l_ec = 1.01/2^k), solve the BIE on the single dielectric cube
and evaluate the far-field single-layer potential on a sphere of radius
R_far outside the cube. An extra deeper level k = 11 is computed at p = 6
to serve as the common reference; relative L² errors of u_ℓ vs u_ref are
plotted in panel (b) as a function of degrees of freedom.

All u_far vectors are serialized, so the reference can be switched in the
plotter without re-running the sweep.

Independent of fig3_data.jl: this script supersedes its panel (b) data.
=#

using LinearAlgebra
using FastGaussQuadrature
using Serialization
using BoundaryIntegral
const BI = BoundaryIntegral

const Lx, Ly, Lz = 2.0, 2.0, 2.0
const eps_d, eps_0 = 10.0, 1.0
const eps_src = eps_0
const ps = BI.PointSource((0.5, 0.6, 100.0), 1.0e4)
const fmm_tol  = 1e-6
const gmres_tol = 1e-6

const p_list = [2, 3, 4, 5, 6]
const k_base = 1:10
const ref_p  = 6
const ref_k  = 11   # extra deep level at p_ref serving as reference

# Per-p depth list: p=6 gets one extra (deeper) level for the reference
k_list_for(p) = (p == ref_p) ? collect(1:ref_k) : collect(k_base)

const R_far = 5.0
const n_theta = 12
const n_phi   = 24
far_targets = Matrix{Float64}(undef, 3, n_theta * n_phi)
let kk = 1
    for i in 1:n_theta
        θ = π * (i - 0.5) / n_theta
        for j in 1:n_phi
            φ = 2π * (j - 1) / n_phi
            far_targets[1, kk] = R_far * sin(θ) * cos(φ)
            far_targets[2, kk] = R_far * sin(θ) * sin(φ)
            far_targets[3, kk] = R_far * cos(θ)
            kk += 1
        end
    end
end

function solve_and_far(p_quad::Int, l_ec::Float64)
    iface = single_dielectric_box3d(Lx, Ly, Lz, p_quad, l_ec, eps_d, eps_0)
    Lhs = lhs_dielectric_box3d_fmm3d(iface, fmm_tol)
    rhs = rhs_dielectric_box3d(iface, ps, eps_src)
    sigma = BI.solve_gmres(Lhs, rhs, gmres_tol, gmres_tol)
    S = BI.laplace3d_pottrg(iface, far_targets)
    u_far = S * sigma
    return length(iface.panels), u_far
end

edge_results = NamedTuple[]
for p in p_list
    @info "===== p sweep =====" p
    levels = NamedTuple[]
    for k in k_list_for(p)
        l_ec = 1.01 / 2.0^k
        @info "  level" p k l_ec
        t0 = time()
        npan, u_far = solve_and_far(p, l_ec)
        dt = time() - t0
        N = npan * p^2
        push!(levels, (k = k, l_ec = l_ec, npanels = npan, N = N, u_far = u_far))
        @info "  done" p k npan N dt
    end
    push!(edge_results, (p = p, levels = levels))
end

out = (
    p_list = p_list,
    ref_p = ref_p,
    ref_k = ref_k,
    far_R = R_far,
    edge_results = edge_results,
)

path = joinpath(@__DIR__, "fig3_data_pscan.jls")
open(io -> serialize(io, out), path, "w")
@info "Saved" path bytes=stat(path).size

println("\n=== Multi-p edge-local sweep (reference: p=$ref_p, k=$ref_k) ===")
for r in edge_results
    println("\n  p = $(r.p)")
    println("    k    l_ec        N")
    for l in r.levels
        println("    $(l.k)    $(round(l.l_ec; sigdigits=3))   $(l.N)")
    end
end

#=
Multi-p sweep for fig3 panel (b)  --  EDGE-CORRECTED solve.

For each polynomial order p ∈ {2,3,4,5,6} and each edge-refinement level
k = 1..8 (l_ec = 1.01/2^k), solve the dielectric-cube BIE with the
near-corrected DT operator (FMM base + sparse near correction + adaptive
quadtree on edge-touching panel pairs, i.e. correct_edges=true), then
evaluate the far-field single-layer potential on a sphere of radius R_far
outside the cube. p=6 gets one extra level k=9 to serve as the common
reference; relative L² errors of u_ℓ vs u_ref are plotted in panel (b) as a
function of degrees of freedom.

Why correct_edges=true: without it the edge-region quadrature error decays
only ~2^-k and dominates the metric (it is NOT the geometric edge
singularity we want to study). The adaptive correction removes that floor so
the convergence reflects σ resolution. The previous (uncorrected) version of
this figure needed k up to 13 to fight that floor; the corrected solve
converges by k≈8, so the deep 18.8M-DOF reference is no longer needed.

All u_far vectors are serialized, so the reference can be switched in the
plotter without re-running the sweep. CHECKPOINTING: the .jls is rewritten
after every completed (p,k) level and a flushed line appended to the
progress file, so a long deep level loses nothing on interruption.
=#

using LinearAlgebra
using FastGaussQuadrature
using Serialization
using Printf
using BoundaryIntegral
const BI = BoundaryIntegral

const Lx, Ly, Lz = 2.0, 2.0, 2.0
const eps_d, eps_0 = 10.0, 1.0
const eps_src = eps_0
const ps = BI.PointSource((0.5, 0.6, 100.0), 1.0e4)
const fmm_tol   = 1e-6
const gmres_tol = 1e-6
const up_tol    = 1e-6   # near-correction upsample atol
const max_order = 64     # cap on near-correction GL order

const p_list = [2, 3, 4, 5, 6]
const k_base = 1:8
const ref_p  = 6
const ref_k  = 9    # extra deep level at p_ref serving as reference

# Per-p depth list: p=ref_p gets one extra (deeper) level for the reference
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

const jls_path  = joinpath(@__DIR__, "fig3_data_pscan.jls")
const prog_path = joinpath(@__DIR__, "fig3_data_pscan_progress.txt")

function solve_and_far(p_quad::Int, l_ec::Float64)
    iface = single_dielectric_box3d(Lx, Ly, Lz, p_quad, l_ec, eps_d, eps_0)
    Lhs = lhs_dielectric_box3d_fmm3d_corrected(iface, fmm_tol, up_tol, max_order;
                                               correct_edges = true)
    rhs = rhs_dielectric_box3d(iface, ps, eps_src)
    sigma = BI.solve_gmres(Lhs, rhs, gmres_tol, gmres_tol)
    gres = norm(Lhs * sigma - rhs) / norm(rhs)
    S = BI.laplace3d_pottrg(iface, far_targets)
    u_far = S * sigma
    return length(iface.panels), u_far, gres
end

# checkpoint: rewrite the .jls with everything collected so far
function checkpoint(edge_results)
    out = (
        p_list = p_list,
        ref_p = ref_p,
        ref_k = ref_k,
        far_R = R_far,
        correct_edges = true,
        max_order = max_order,
        edge_results = edge_results,
    )
    open(io -> serialize(io, out), jls_path, "w")
end

open(prog_path, "w") do io
    println(io, "# fig3 panel(b) multi-p sweep  (correct_edges=true, max_order=$max_order)")
    println(io, "# fmm_tol=$fmm_tol gmres_tol=$gmres_tol up_tol=$up_tol  ref=(p=$ref_p,k=$ref_k)")
    println(io, "# p   k   l_ec        npanels   N          gmres_res     dt[s]")
    flush(io)
end

edge_results = NamedTuple[]
for p in p_list
    @info "===== p sweep =====" p
    levels = NamedTuple[]
    for k in k_list_for(p)
        l_ec = 1.01 / 2.0^k
        @info "  level" p k l_ec
        t0 = time()
        npan, u_far, gres = solve_and_far(p, l_ec)
        dt = time() - t0
        N = npan * p^2
        push!(levels, (k = k, l_ec = l_ec, npanels = npan, N = N, u_far = u_far,
                       gmres_res = gres, dt = dt))
        @info "  done" p k npan N gres dt

        # checkpoint after every completed level (rewrite full results so far)
        checkpoint(vcat(edge_results, [(p = p, levels = levels)]))
        open(prog_path, "a") do io
            @printf(io, "  %d   %d   %-9.4g   %-7d   %-9d   %-11.3g   %.1f\n",
                    p, k, l_ec, npan, N, gres, dt)
            flush(io)
        end
    end
    push!(edge_results, (p = p, levels = levels))
end

checkpoint(edge_results)
@info "Saved" jls_path bytes=stat(jls_path).size

println("\n=== Multi-p edge-corrected sweep (reference: p=$ref_p, k=$ref_k) ===")
for r in edge_results
    println("\n  p = $(r.p)")
    println("    k    l_ec        N          gmres_res")
    for l in r.levels
        println("    $(l.k)    $(round(l.l_ec; sigdigits=3))   $(l.N)    $(round(l.gmres_res; sigdigits=3))")
    end
end

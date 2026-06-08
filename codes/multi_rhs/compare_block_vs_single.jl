# Run a .bie center group and compute V_ijkl two ways on the SAME interface / operator / RHS:
#   (A) block GMRES  — all K layer densities solved together (the multi-RHS path)
#   (B) K separate single-RHS GMRES solves
# Then compare a subset of V (the center-on-site row V[1, :]) block-vs-separate, and report the
# solve speedup of (A) over (B).
#
# Usage: julia --project scripts/compare_block_vs_single.jl [bie_path] [center_id]
# Defaults: scripts/graphene_cell.bie, center 1.

using BoundaryIntegral
using Krylov, LinearAlgebra, Printf
const BI = BoundaryIntegral

bie    = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "graphene_cell.bie")
center = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

si = read_system_input(bie)
sp = si.solve
group   = BI.assemble_rhs_group(si, center; support_rtol = sp.support_rtol)
sources = group_volume_sources(group)
labels  = ["rho_$(center)_$(j)" for j in group.neighbor_ids]
K = length(sources)
@printf("center %d:  K=%d  union-support points=%d\n", center, K, size(group.positions, 2))

interface = build_group_interface(si, group; n_quad = sp.n_quad, rhs_atol = sp.rhs_tol,
                                   l_ec = resolved_l_ec(si), eps_out = si.eps_out, max_depth = sp.max_depth)
@printf("interface: %d panels, %d points\n", length(interface.panels), BI.num_points(interface))

op = batched_lhs_dielectric_box3d_fmm3d_corrected(interface, sp.lhs_tol, sp.lhs_tol, sp.max_order)
F  = rhs_dielectric_box3d_fmm3d(interface, sources, sp.lhs_tol)

# warm up the matvec (so timings exclude first-call compilation)
op * F[:, 1]

# (A) block multi-RHS solve
t0 = time_ns()
Σ_block, bstats = Krylov.block_gmres(op, F; rtol = sp.gmres_rtol, itmax = 500)
t_block = (time_ns() - t0) / 1e9

# (B) K separate single-RHS solves
Σ_single = similar(Σ_block)
single_iters = Int[]
t0 = time_ns()
for k in 1:K
    xk, st = Krylov.gmres(op, F[:, k]; rtol = sp.gmres_rtol, itmax = 500)
    Σ_single[:, k] .= xk
    push!(single_iters, st.niter)
end
t_single = (time_ns() - t0) / 1e9
@printf("\nSOLVE:  block=%.2fs (%d iters)   separate=%.2fs (iters %s)   speedup=%.2fx\n",
        t_block, bstats.niter, t_single, single_iters, t_single / t_block)

# shared V evaluation: incident potential (TKM) per source + layer-potential operator (built once)
targets = sources[1].positions
pottrg  = BI.laplace3d_pottrg_fmm3d_corrected_hcubature(interface, targets, sp.lhs_tol, sp.lhs_tol, 5.0)
u_inc = Vector{Vector{Float64}}(undef, K)
for b in 1:K
    sb = BI.screened_volume_source(interface, sources[b], BI.SharpScreening())
    vals = BI.TKM3D.ltkm3dc(sp.volume_tol, sb.positions; charges = sb.weights .* sb.density,
                            targets = targets, pgt = 1, kmax = BI._estimate_tkm3dc_kmax(sb))
    u_inc[b] = real.(vals.pottarg)
end
tw = [sources[a].weights .* sources[a].density for a in 1:K]
Vmat(Σ) = begin
    V = Matrix{Float64}(undef, K, K)
    for b in 1:K
        φb = u_inc[b] .+ (pottrg * Σ[:, b])
        for a in 1:K
            V[a, b] = dot(tw[a], φb)
        end
    end
    V
end
V_block  = Vmat(Σ_block)
V_single = Vmat(Σ_single)

println("\nSubset V[1, b]  (center on-site vs each group member): block vs separate")
@printf("  %-12s %18s %18s %10s\n", "pair", "V_block", "V_separate", "rel|Δ|")
for b in 1:K
    rel = abs(V_block[1, b] - V_single[1, b]) / max(abs(V_single[1, b]), eps())
    @printf("  %-12s % .10e % .10e %10.2e\n", labels[b], V_block[1, b], V_single[1, b], rel)
end
@printf("\nmax rel |ΔV| over full K×K matrix = %.3e\n",
        maximum(abs.(V_block .- V_single) ./ (abs.(V_single) .+ eps())))
println("done.")

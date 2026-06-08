# Per-RHS runtime benchmark for the multi-RHS (block GMRES) solve.
#
# Fixed-interface K-sweep: build ONE interface/operator/RHS for a large neighbor set of a single
# center (graphene A + both-sublattice lattice images over 2 shells), then block-solve the K
# nearest right-hand sides for K = 1,2,4,8,16,Kmax. Reports the per-RHS solve time t_block(K)/K
# (which should fall as K grows: shared Krylov space + one nd=K batched FMM per matvec) against
# the single-RHS baseline.
#
# Usage: julia --project scripts/bench_per_rhs.jl

using BoundaryIntegral
using Krylov, LinearAlgebra, Printf
const BI = BoundaryIntegral

const DATA = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const A = joinpath(DATA, "graphene_00001.xsf")
const B = joinpath(DATA, "graphene_00002.xsf")

# --- generate a 2-shell graphene .bie (center A home + A,B images over (n1,n2) in [-2,2]^2) ---
biepath = joinpath(@__DIR__, "graphene_2shell.bie")
open(biepath, "w") do io
    println(io, "UNITS bohr\n")
    println(io, "BEGIN_DIELECTRICS\nEPS_OUT 1.0\n  0.0 0.0 7.5    90.0 90.0 3.35    3.5\nEND_DIELECTRICS\n")
    println(io, "BEGIN_ORBITALS")
    id = 0
    id += 1; println(io, "  $id   $A")
    id += 1; println(io, "  $id   $B")
    for n1 in -2:2, n2 in -2:2
        (n1 == 0 && n2 == 0) && continue
        id += 1; println(io, "  $id   $A   LATTICE $n1 $n2 0")
        id += 1; println(io, "  $id   $B   LATTICE $n1 $n2 0")
    end
    println(io, "END_ORBITALS\n")
    println(io, "BEGIN_GROUPING\nCUTOFF 5.0\nEND_GROUPING\n")
    println(io, "BEGIN_SOLVE")
    println(io, "  N_QUAD 6\n  EDGE_REFINE_LEVEL 2\n  RHS_TOL 1e-3\n  LHS_TOL 1e-5")
    println(io, "  GMRES_RTOL 1e-5\n  SUPPORT_RTOL 1e-4\n  VOLUME_TOL 1e-5\nEND_SOLVE")
end

si = read_system_input(biepath)
sp = si.solve
group = BI.assemble_rhs_group(si, 1; support_rtol = sp.support_rtol)

# order the group's densities by orbital distance from center 1 (nearest first), so F[:,1:K]
# is the K nearest neighbors.
c1 = si.orbitals[1].center
perm = sortperm([norm(si.orbitals[j].center .- c1) for j in group.neighbor_ids])
sources = group_volume_sources(group)[perm]
Kmax = length(sources)
@printf("center 1: Kmax=%d  union-support points=%d\n", Kmax, size(group.positions, 2))

interface = build_group_interface(si, group; n_quad = sp.n_quad, rhs_atol = sp.rhs_tol,
                                   l_ec = resolved_l_ec(si), eps_out = si.eps_out, max_depth = sp.max_depth)
@printf("interface: %d panels, %d points\n", length(interface.panels), BI.num_points(interface))

op = batched_lhs_dielectric_box3d_fmm3d_corrected(interface, sp.lhs_tol, sp.lhs_tol, sp.max_order)
@printf("operator: %d rows, %d cols\n", size(op, 1), size(op, 2))
F  = rhs_dielectric_box3d_fmm3d(interface, sources, sp.lhs_tol)     # N x Kmax (nearest-first)
@printf("RHS: %d rows, %d cols\n", size(F, 1), size(F, 2))

op * F[:, 1]   # warm up the matvec / compilation

# --- single-RHS baseline: solve all Kmax separately ---
single_iters = Int[]
t0 = time_ns()
for k in 1:Kmax
    println("single-RHS solve $k / $Kmax")
    _, st = Krylov.gmres(op, F[:, k]; rtol = sp.gmres_rtol, itmax = 500, verbose = 1)
    push!(single_iters, st.niter)
end
t_single = (time_ns() - t0) / 1e9
per_rhs_single = t_single / Kmax

# --- block sweep ---
Ks = sort(unique(filter(K -> 1 <= K <= Kmax, [1, 2, 4, 8, 16, Kmax])))
println("\n  K   t_block(s)   per-RHS(s)   block-iters   speedup-vs-single")
for K in Ks
    println("block solve $K / $Kmax")
    t0 = time_ns()
    _, bst = Krylov.block_gmres(op, F[:, 1:K]; rtol = sp.gmres_rtol, itmax = 500, verbose = 1)
    tK = (time_ns() - t0) / 1e9
    @printf("%4d  %9.2f   %9.3f   %8d        %6.2fx\n", K, tK, tK / K, bst.niter, per_rhs_single / (tK / K))
end
@printf("\nsingle-RHS baseline: per-RHS = %.3f s  (total %.1f s over %d solves, iters %s)\n",
        per_rhs_single, t_single, Kmax, single_iters)
println("done.")

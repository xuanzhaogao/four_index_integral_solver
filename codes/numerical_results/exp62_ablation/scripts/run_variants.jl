# 6.2.1 / 6.2.2 — component ablation on the gap-stressed slab (A = 10).
# Geometry: 10 x 10 x 0.5 slab (eps 10) + 0.6 x 0.6 x 0.2 material box (eps 2)
# at gap 0.05 above the slab center; Gaussian source (s = 0.05) inside the box.
# This stresses the near correction across the GAP (parallel close faces),
# a different mechanism from the edge handling studied in 6.1.
#
# Variants, identical in everything else (same adaptive interface as A unless noted):
#   A   full method (Bernstein upsampling + adaptive edge-pair correction)
#   B   no near correction at all (D^{T,corr} = 0)
#   Bp  upsampling on, edge correction off  (isolates the two correction kinds)
#   C   no RHS adaptivity: greedy quasi-uniform dyadic refinement DOF-matched
#       to A (within the +3p^2 granularity of one split), same edge refinement,
#       full corrections
# Sweep eps in {1e-3, 1e-4, 1e-5, 1e-6} at p = 6, r = 4 (regime set with 6.1).
# Reference: p = 8, eps = 1e-7, r = 6, finer source/TKM grid, edge-corrected.
# Variant B's predicted stagnation rho_min^(-2p) is logged from A's neighbor list.
#
# Run:  cd numerical_results
#       GMRES_VERBOSE=0 JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 \
#         julia --project=. exp62_ablation/scripts/run_variants.jl

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using LinearAlgebra, Printf

const SYS = slab_system()
const DATA = joinpath(@__DIR__, "..", "data")
const CSVPATH = joinpath(DATA, "variants.csv")

const EPS_LIST = [1e-3, 1e-4, 1e-5, 1e-6]
const P, R = 6, 4
const REF = (p = 8, eps = 1e-7, r = 6, margin = 1.4)
const GMRES_VERBOSE = parse(Int, get(ENV, "GMRES_VERBOSE", "0"))

function log_run(res, variant, eps, V, wall; extra...)
    append_csv_row(CSVPATH, run_cols(res; V = V, variant = variant,
                   wall = round(wall; digits = 2), extra...))
    @printf("[6.2] %-3s eps=%.0e N=%d niter=%d V=%.10e wall=%.0fs\n",
            variant, eps, res.N, res.niter, V, wall)
    flush(stdout)
end

println(">>> warm-up"); flush(stdout)
solve_system(slab_internal(); eps = 1e-2, p = 2, r = 1)

# reference (edge-corrected; differs from all tests in p, eps, r, grids)
if !isfile(joinpath(DATA, "raw", "slab_ref.jls"))
    println(">>> reference p=$(REF.p) eps=$(REF.eps) r=$(REF.r)"); flush(stdout)
    t0 = time()
    res = solve_system(SYS; eps = REF.eps, p = REF.p, r = REF.r,
                       src_margin = REF.margin, gmres_verbose = GMRES_VERBOSE)
    V = eval_V(res; t_out = res.times, margin = REF.margin)
    save_ref(joinpath(DATA, "raw", "slab_ref.jls"),
             (; V, N = res.N, niter = res.niter, eps = REF.eps, p = REF.p, r = REF.r))
    log_run(res, "ref", REF.eps, V, time() - t0)
end

for eps in EPS_LIST
    donefile = joinpath(DATA, "raw", @sprintf("variants_eps%.0e.jls", eps))
    isfile(donefile) && (println("eps=$eps exists, skipping"); continue)

    # --- A: full method
    t0 = time()
    resA = solve_system(SYS; eps = eps, p = P, r = R, gmres_verbose = GMRES_VERBOSE)
    VA = eval_V(resA; t_out = resA.times)
    wallA = time() - t0
    rho_min = bernstein_rho_min(resA.interface, resA.nb)
    log_run(resA, "A", eps, VA, wallA; rho_min = round(rho_min; digits = 5),
            stagnation_est = @sprintf("%.3e", rho_min^(-2P)))

    # --- B: same interface, no near correction at all
    t0 = time()
    resB = solve_system(SYS; eps = eps, p = P, r = R,
                        near_correction = false,
                        interface_override = resA.interface,
                        screened_vs_override = resA.screened_vs,
                        tkm_kmax_override = resA.tkm_kmax,
                        gmres_verbose = GMRES_VERBOSE)
    VB = eval_V(resB; t_out = resB.times)
    log_run(resB, "B", eps, VB, time() - t0)

    # --- Bp: upsampling on, edge correction off
    t0 = time()
    resBp = solve_system(SYS; eps = eps, p = P, r = R,
                         correct_edges = false,
                         interface_override = resA.interface,
                         screened_vs_override = resA.screened_vs,
                         tkm_kmax_override = resA.tkm_kmax,
                         gmres_verbose = GMRES_VERBOSE)
    VBp = eval_V(resBp; t_out = resBp.times)
    log_run(resBp, "Bp", eps, VBp, time() - t0)

    # --- C: quasi-uniform refinement DOF-matched to A, same l_ec, full corrections
    t0 = time()
    base = BI.multi_dielectric_box3d(P, Harness.l_ec_of(SYS, R), SYS.boxes, SYS.epses, SYS.eps_out)
    ifaceC = refine_to_dof(base, resA.N)
    resC = solve_system(SYS; eps = eps, p = P, r = R,
                        interface_override = ifaceC,
                        screened_vs_override = resA.screened_vs,
                        tkm_kmax_override = resA.tkm_kmax,
                        gmres_verbose = GMRES_VERBOSE)
    VC = eval_V(resC; t_out = resC.times)
    log_run(resC, "C", eps, VC, time() - t0; dof_match = round(resC.N / resA.N; digits = 3))

    save_ref(donefile, (; eps, VA, VB, VBp, VC, NA = resA.N, NC = resC.N,
             niterA = resA.niter, niterB = resB.niter, niterBp = resBp.niter,
             niterC = resC.niter,
             rho_min, n_near = resA.n_near_pairs, n_adaptive = resA.n_adaptive_pairs,
             p_up_max = resA.p_up_max,
             timesA = resA.times, timesB = resB.times, timesBp = resBp.times,
             timesC = resC.times, corr_bytes = resA.corr_bytes))
end

println("VARIANTS DONE")

# 6.2.1 / 6.2.2 — component ablation on the high-aspect-ratio slab (A = 10).
# Variants, identical in everything else:
#   A  full method
#   B  no near correction (D^{T,corr} = 0), same adaptive interface as A
#   C  no RHS adaptivity: greedy quasi-uniform dyadic refinement DOF-matched
#      to A (within the +3p^2 granularity of one split), same edge refinement
# Sweep eps in {1e-3..1e-11} at p = 8. Reference: p=12, eps=1e-13, r=6, finer grids.
# The Bernstein stagnation estimate rho_min^(-2p) is logged from variant A's
# neighbor list at each eps.

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using LinearAlgebra, Printf

const SYS = slab_system()        # 10 x 10 x 0.5 slab (eps 10) + material box (eps 2), gap 0.05
const DATA = joinpath(@__DIR__, "..", "data")
const CSVPATH = joinpath(DATA, "variants.csv")

const EPS_LIST = [1e-3, 1e-5, 1e-7, 1e-9, 1e-11]
const P, R = 8, 4

function log_run(res, variant, eps, V, wall; extra...)
    append_csv_row(CSVPATH, run_cols(res; V = V, variant = variant,
                   wall = round(wall; digits = 2), extra...))
    @printf("[6.2] %s eps=%.0e N=%d niter=%d V=%.10e wall=%.0fs\n",
            variant, eps, res.N, res.niter, V, wall)
    flush(stdout)
end

println(">>> warm-up")
solve_system(system1(); eps = 1e-2, p = 4, r = 1)

# reference
if !isfile(joinpath(DATA, "raw", "slab_ref.jls"))
    println(">>> reference p=12 eps=1e-13 r=6")
    t0 = time()
    res = solve_system(SYS; eps = 1e-13, p = 12, r = R + 2, src_margin = 1.4)
    V = eval_V(res; t_out = res.times, margin = 1.4)
    save_ref(joinpath(DATA, "raw", "slab_ref.jls"), (; V, N = res.N, niter = res.niter, eps = 1e-13, p = 12, r = R + 2))
    log_run(res, "ref", 1e-13, V, time() - t0)
end

for eps in EPS_LIST
    donefile = joinpath(DATA, "raw", @sprintf("variants_eps%.0e.jls", eps))
    isfile(donefile) && (println("eps=$eps exists, skipping"); continue)

    # --- A: full method
    t0 = time()
    resA = solve_system(SYS; eps = eps, p = P, r = R)
    VA = eval_V(resA; t_out = resA.times)
    wallA = time() - t0
    rho_min = bernstein_rho_min(resA.interface, resA.nb)
    log_run(resA, "A", eps, VA, wallA; rho_min = round(rho_min; digits = 5),
            stagnation_est = @sprintf("%.3e", rho_min^(-2P)))

    # --- B: same interface, no near correction
    t0 = time()
    resB = solve_system(SYS; eps = eps, p = P, r = R,
                        near_correction = false,
                        interface_override = resA.interface,
                        screened_vs_override = resA.screened_vs,
                        tkm_kmax_override = resA.tkm_kmax)
    VB = eval_V(resB; t_out = resB.times)
    log_run(resB, "B", eps, VB, time() - t0)

    # --- C: quasi-uniform refinement DOF-matched to A, same l_ec
    t0 = time()
    base = BI.multi_dielectric_box3d(P, l_ec_of(SYS, R), SYS.boxes, SYS.epses, SYS.eps_out)
    ifaceC = refine_to_dof(base, resA.N)
    resC = solve_system(SYS; eps = eps, p = P, r = R,
                        interface_override = ifaceC,
                        screened_vs_override = resA.screened_vs,
                        tkm_kmax_override = resA.tkm_kmax)
    VC = eval_V(resC; t_out = resC.times)
    log_run(resC, "C", eps, VC, time() - t0; dof_match = round(resC.N / resA.N; digits = 3))

    save_ref(donefile, (; eps, VA, VB, VC, NA = resA.N, NC = resC.N,
             niterA = resA.niter, niterB = resB.niter, niterC = resC.niter,
             rho_min, n_near = resA.n_near_pairs, p_up_max = resA.p_up_max,
             timesA = resA.times, timesB = resB.times, timesC = resC.times,
             corr_bytes = resA.corr_bytes))
end

println("VARIANTS DONE")

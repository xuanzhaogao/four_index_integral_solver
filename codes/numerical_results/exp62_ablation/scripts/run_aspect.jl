# 6.2.3 — aspect-ratio sweep: slab A x A x 0.5, A in {1, 5, 10, 20, 50},
# fixed gap 0.05, eps = 1e-9, p = 8. Variants A (full) and B (no correction),
# with a per-aspect reference (p=12, eps=1e-13, r=6, finer grids).

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using Printf

const DATA = joinpath(@__DIR__, "..", "data")
const CSVPATH = joinpath(DATA, "aspect.csv")
const EPS, P, R = 1e-9, 8, 4
const A_LIST = [1.0, 5.0, 10.0, 20.0, 50.0]

println(">>> warm-up")
solve_system(system1(); eps = 1e-2, p = 4, r = 1)

for A in A_LIST
    donefile = joinpath(DATA, "raw", "aspect_A$(A).jls")
    isfile(donefile) && (println("A=$A exists, skipping"); continue)
    sys = slab_system(A = A)

    t0 = time()
    resA = solve_system(sys; eps = EPS, p = P, r = R)
    VA = eval_V(resA; t_out = resA.times)
    wallA = time() - t0
    rho_min = bernstein_rho_min(resA.interface, resA.nb)
    append_csv_row(CSVPATH, run_cols(resA; V = VA, variant = "A", aspect = A,
                   wall = round(wallA; digits = 2), rho_min = round(rho_min; digits = 5)))
    @printf("[aspect] A=%g full   N=%d niter=%d near=%d p_up_max=%d V=%.10e (%.0fs)\n",
            A, resA.N, resA.niter, resA.n_near_pairs, resA.p_up_max, VA, wallA)
    flush(stdout)

    t0 = time()
    resB = solve_system(sys; eps = EPS, p = P, r = R, near_correction = false,
                        interface_override = resA.interface,
                        screened_vs_override = resA.screened_vs,
                        tkm_kmax_override = resA.tkm_kmax)
    VB = eval_V(resB; t_out = resB.times)
    append_csv_row(CSVPATH, run_cols(resB; V = VB, variant = "B", aspect = A,
                   wall = round(time() - t0; digits = 2)))

    t0 = time()
    resR = solve_system(sys; eps = 1e-13, p = 12, r = R + 2, src_margin = 1.4)
    VR = eval_V(resR; t_out = resR.times, margin = 1.4)
    append_csv_row(CSVPATH, run_cols(resR; V = VR, variant = "ref", aspect = A,
                   wall = round(time() - t0; digits = 2)))
    @printf("[aspect] A=%g ref    N=%d niter=%d V=%.12e\n", A, resR.N, resR.niter, VR)
    flush(stdout)

    save_ref(donefile, (; A, VA, VB, VR, NA = resA.N, NR = resR.N,
             niterA = resA.niter, niterB = resB.niter,
             n_near = resA.n_near_pairs, p_up_max = resA.p_up_max,
             p_up_mean = resA.p_up_mean, rho_min, corr_bytes = resA.corr_bytes))
end

println("ASPECT SWEEP DONE")

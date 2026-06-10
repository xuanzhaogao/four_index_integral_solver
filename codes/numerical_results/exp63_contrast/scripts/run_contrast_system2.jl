# 6.3.3 — System II with (eps1, eps2) jointly scaled by {1, 2, 8} (eps_m = 2 fixed).
# Iteration behavior with multiple coexisting gamma values; one N_iter table.

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
using Printf

const DATA = joinpath(@__DIR__, "..", "data")
const CSVPATH = joinpath(DATA, "contrast_system2.csv")
const EPS, P, R = 1e-9, 8, 4

println(">>> warm-up")
solve_system(system1(); eps = 1e-2, p = 4, r = 1)

for k in (1.0, 2.0, 8.0)
    donefile = joinpath(DATA, "raw", "contrast2_k$(k).jls")
    isfile(donefile) && continue
    sys = system2(eps1 = 4.0 * k, eps2 = 12.0 * k)
    t0 = time()
    res = solve_system(sys; eps = EPS, p = P, r = R)
    V = eval_V(res; t_out = res.times)
    wall = time() - t0
    g1 = (4k - 1) / (4k + 1); g2 = (12k - 1) / (12k + 1); g12 = (4k - 12k) / (4k + 12k)
    append_csv_row(CSVPATH, run_cols(res; V = V, scale = k, eps1 = 4k, eps2 = 12k,
                   gamma1 = round(g1; digits = 4), gamma2 = round(g2; digits = 4),
                   gamma12 = round(g12; digits = 4), wall = round(wall; digits = 2)))
    @printf("[6.3.3] k=%g N=%d niter=%d V=%.10e (%.0fs)\n", k, res.N, res.niter, V, wall)
    flush(stdout)
    save_ref(donefile, (; k, V, N = res.N, niter = res.niter, history = res.gmres_history))
end

println("CONTRAST SYSTEM2 DONE")

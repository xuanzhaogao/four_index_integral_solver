# 6.3.1 / 6.3.2 — dielectric contrast sweep on System I (eps = 1e-9, p = 8).
# eps1 in {2, 5, 10, 20, 80, 1000} + 1e12 (gamma -> 1 grounded-conductor limit).
# Per-contrast reference (p = 12, eps = 1e-13, r = 6, finer grids).
# Also stores GMRES residual histories and V_vacuum (direct TKM, no interface).

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
using Printf

const DATA = joinpath(@__DIR__, "..", "data")
const CSVPATH = joinpath(DATA, "contrast.csv")
const EPS, P, R = 1e-9, 8, 4
const EPS1_LIST = [2.0, 5.0, 10.0, 20.0, 80.0, 1000.0, 1e12]

println(">>> warm-up")
solve_system(system1(); eps = 1e-2, p = 4, r = 1)

# vacuum interaction (no interface), high accuracy
let f = joinpath(DATA, "raw", "v_vacuum.jls")
    if !isfile(f)
        Vvac = vacuum_V(system1(), 1e-13)
        save_ref(f, (; Vvac))
        @printf("V_vacuum = %.12e\n", Vvac)
    end
end

for eps1 in EPS1_LIST
    donefile = joinpath(DATA, "raw", "contrast_eps1_$(eps1).jls")
    isfile(donefile) && (println("eps1=$eps1 exists, skipping"); continue)
    sys = system1(eps1 = eps1)
    gamma = (eps1 - 1.0) / (eps1 + 1.0)

    t0 = time()
    res = solve_system(sys; eps = EPS, p = P, r = R)
    V = eval_V(res; t_out = res.times)
    wall = time() - t0
    append_csv_row(CSVPATH, run_cols(res; V = V, eps1 = eps1, gamma = gamma,
                   variant = "test", wall = round(wall; digits = 2)))
    @printf("[contrast] eps1=%g gamma=%.6f N=%d niter=%d V=%.10e (%.0fs)\n",
            eps1, gamma, res.N, res.niter, V, wall)
    flush(stdout)

    t0 = time()
    resR = solve_system(sys; eps = 1e-13, p = 12, r = R + 2, src_margin = 1.4)
    VR = eval_V(resR; t_out = resR.times, margin = 1.4)
    append_csv_row(CSVPATH, run_cols(resR; V = VR, eps1 = eps1, gamma = gamma,
                   variant = "ref", wall = round(time() - t0; digits = 2)))
    @printf("[contrast] eps1=%g ref N=%d niter=%d V=%.12e\n", eps1, resR.N, resR.niter, VR)
    flush(stdout)

    save_ref(donefile, (; eps1, gamma, V, VR, N = res.N, niter = res.niter,
             niter_ref = resR.niter, history = res.gmres_history,
             residual = res.residual))
end

println("CONTRAST DONE")

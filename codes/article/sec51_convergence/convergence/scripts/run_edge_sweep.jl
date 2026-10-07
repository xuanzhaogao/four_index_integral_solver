# 6.1.2 — edge-refinement isolation on System I.
# p = 8, all tolerances at the eps = 1e-11 setting; sweep edge depth r = 0..8 alone.
# Reference: p = 12, eps = 1e-13, r = 10, doubled source grid.

include(joinpath(@__DIR__, "..", "..", "..", "common", "Harness.jl"))
using .Harness
using Printf

const SYS = system1()
const DATA = joinpath(@__DIR__, "..", "data")
const CSVPATH = joinpath(DATA, "edge_sweep.csv")

const EPS = 1e-11
const P = 8
const R_LIST = 0:8

function one_run(eps, p, r; margin = 1.25, tag = "run")
    t_all = time()
    res = solve_system(SYS; eps = eps, p = p, r = r, src_margin = margin)
    V = eval_V(res; t_out = res.times, margin = margin)
    wall = time() - t_all
    save_ref(joinpath(DATA, "raw", "edge_$(tag).jls"),
             (; eps, p, r, margin, V, N = res.N, niter = res.niter,
                residual = res.residual, times = res.times, wall))
    append_csv_row(CSVPATH, run_cols(res; V = V, wall = round(wall; digits = 2), tag = tag))
    @printf("[edge] %s r=%d N=%d niter=%d V=%.12e wall=%.1fs\n", tag, r, res.N, res.niter, V, wall)
    flush(stdout)
end

println(">>> warm-up")
solve_system(system1(); eps = 1e-2, p = 4, r = 1)

if !isfile(joinpath(DATA, "raw", "edge_ref.jls"))
    println(">>> reference: p=12 eps=1e-13 r=10")
    one_run(1e-13, 12, 10; margin = 1.4, tag = "ref")
end

for r in R_LIST
    tag = "r$(r)"
    isfile(joinpath(DATA, "raw", "edge_$(tag).jls")) && continue
    one_run(EPS, P, r; tag = tag)
end

println("EDGE SWEEP DONE")

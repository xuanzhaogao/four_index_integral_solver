# 6.1.1 / 6.1.4 — convergence sweep eps x p for System I or System II.
# Usage: julia --project=. exp61_convergence/scripts/run_sweep.jl system1|system2
#
# Reference run first (p=12, eps=1e-13, r = test_r + 2, doubled source/TKM grid),
# then the eps x p sweep. Raw values stored per run (.jls); errors are computed
# at analysis time. CSV rows appended incrementally (crash-safe).

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using LinearAlgebra, Printf

const WHICH = isempty(ARGS) ? "system1" : ARGS[1]
const SYS = WHICH == "system1" ? system1() : system2()
const DATA = joinpath(@__DIR__, "..", "data")
const CSVPATH = joinpath(DATA, "sweep_$(WHICH).csv")

const EPS_LIST = [1e-3, 1e-5, 1e-7, 1e-9, 1e-11]
const P_LIST = [4, 6, 8]
const R_TEST = 4

const REF = (p = 12, eps = 1e-13, r = R_TEST + 2, margin = 1.4)

function one_run(sys, eps, p, r; margin = 1.25, tag = "run")
    t_all = time()
    res = solve_system(sys; eps = eps, p = p, r = r, src_margin = margin)
    tdict = res.times
    V = eval_V(res; t_out = tdict, margin = margin)
    zt = zone_targets(sys)
    phi = Dict{Symbol, Vector{Float64}}()
    for (k, X) in pairs((near = zt.near, supp = zt.supp, far = zt.far))
        phi[k] = eval_phi(res, X; t_out = tdict)
    end
    wall = time() - t_all
    payload = (; sys = sys.name, eps, p, r, margin, V,
               phi_near = phi[:near], phi_supp = phi[:supp], phi_far = phi[:far],
               N = res.N, niter = res.niter, residual = res.residual,
               history = res.gmres_history, times = res.times, wall,
               n_near_pairs = res.n_near_pairs, p_up_max = res.p_up_max,
               p_up_mean = res.p_up_mean, corr_bytes = res.corr_bytes)
    save_ref(joinpath(DATA, "raw", "$(WHICH)_$(tag).jls"), payload)
    append_csv_row(CSVPATH, run_cols(res; V = V, wall = round(wall; digits = 2), tag = tag))
    @printf("[%s] %s eps=%.0e p=%d r=%d N=%d niter=%d V=%.10e wall=%.1fs\n",
            WHICH, tag, eps, p, r, res.N, res.niter, V, wall)
    flush(stdout)
    return payload
end

println(">>> warm-up (compile)")
solve_system(system1(); eps = 1e-2, p = 4, r = 1)

println(">>> reference run: p=$(REF.p) eps=$(REF.eps) r=$(REF.r) margin=$(REF.margin)")
if !isfile(joinpath(DATA, "raw", "$(WHICH)_ref.jls"))
    one_run(SYS, REF.eps, REF.p, REF.r; margin = REF.margin, tag = "ref")
else
    println("    reference exists, skipping")
end

for p in P_LIST, eps in EPS_LIST
    tag = @sprintf("p%d_eps%.0e", p, eps)
    if isfile(joinpath(DATA, "raw", "$(WHICH)_$(tag).jls"))
        println("    $tag exists, skipping"); continue
    end
    one_run(SYS, eps, p, R_TEST; tag = tag)
end

println("SWEEP $(WHICH) DONE")

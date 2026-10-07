# 6.1 Fig.-1 system benchmark: 10x10x1 slab (eps 10) on two 10x10x10 cubes
# (eps 4 | eps 12), Gaussian source (s = 0.05) at the slab center.
# Same protocol as the slab sweep: fixed tied eps = 1e-4,
# sweep p in {2, 4, 6} x r in 1..5; edge-corrected operator throughout.
# Reference: p = 8, eps = 1e-6, r = 6, finer source/TKM grid (margin 1.4).
#
# Per run: V (target Gaussian at (0.2, 0, 0.5)), phi at the three zone sets
# (near plane z = 0.01 spanning cube tops + buried contact faces; source-support
# ball; far sphere), DOF, N_iter, near-pair stats, phase timings.
# Results: data/sweep_fig1.csv + data/raw/fig1_<tag>.jls (crash-safe, resumable).
#
# Run:  cd codes/article
#       GMRES_VERBOSE=0 JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 \
#         julia --project=. sec51_convergence/convergence/scripts/run_fig1_sweep.jl

include(joinpath(@__DIR__, "..", "..", "..", "common", "Harness.jl"))
using .Harness
using Printf

const SYS = system_fig1()
# RERUN_TAG: appends a suffix to this experiment's output directory so a rerun
# never overwrites the data behind the submitted manuscript. Empty = original paths.
const TAG = get(ENV, "RERUN_TAG", "")
const DATA = joinpath(@__DIR__, "..", "data" * TAG)
mkpath(joinpath(DATA, "raw"))
const CSVPATH = joinpath(DATA, "sweep_fig1.csv")

const EPS = 1e-4
const P_LIST = [2, 4, 6]
const R_LIST = 1:5
const REF = (p = 8, eps = 1e-6, r = 6, margin = 1.4)

const GMRES_VERBOSE = parse(Int, get(ENV, "GMRES_VERBOSE", "1"))

function one_run(eps, p, r; margin = 1.25, tag = "run")
    t_all = time()
    res = solve_system(SYS; eps = eps, p = p, r = r, src_margin = margin,
                       gmres_verbose = GMRES_VERBOSE)
    tdict = res.times
    # batched evaluation: V-target grid + all three zone sets in one operator
    vs_t = gaussian_source(SYS.tgt_center, SYS.tgt_sigma, eps; margin = margin)
    zt = zone_targets(SYS)
    nV = size(vs_t.positions, 2)
    X = hcat(Matrix{Float64}(vs_t.positions), zt.near, zt.supp, zt.far)
    phi_all = eval_phi(res, X; t_out = tdict)
    V = sum((vs_t.weights .* vs_t.density) .* phi_all[1:nV])
    nz = size(zt.near, 2)
    phi = Dict(:near => phi_all[(nV + 1):(nV + nz)],
               :supp => phi_all[(nV + nz + 1):(nV + 2nz)],
               :far  => phi_all[(nV + 2nz + 1):(nV + 3nz)])
    wall = time() - t_all
    save_ref(joinpath(DATA, "raw", "fig1_$(tag).jls"),
             (; sys = SYS.name, eps, p, r, margin, V,
                phi_near = phi[:near], phi_supp = phi[:supp], phi_far = phi[:far],
                N = res.N, niter = res.niter, residual = res.residual,
                history = res.gmres_history, times = res.times, wall,
                n_near_pairs = res.n_near_pairs, n_adaptive_pairs = res.n_adaptive_pairs,
                p_up_max = res.p_up_max, p_up_mean = res.p_up_mean,
                corr_bytes = res.corr_bytes))
    append_csv_row(CSVPATH, run_cols(res; V = V, wall = round(wall; digits = 2), tag = tag))
    @printf("[fig1] %s eps=%.0e p=%d r=%d N=%d niter=%d V=%.10e wall=%.1fs\n",
            tag, eps, p, r, res.N, res.niter, V, wall)
    flush(stdout)
end

println(">>> warm-up (compile)"); flush(stdout)
solve_system(slab_internal(); eps = 1e-2, p = 2, r = 1)

if !isfile(joinpath(DATA, "raw", "fig1_ref.jls"))
    println(">>> reference: p=$(REF.p) eps=$(REF.eps) r=$(REF.r) (edge-corrected)"); flush(stdout)
    one_run(REF.eps, REF.p, REF.r; margin = REF.margin, tag = "ref")
end

for p in P_LIST, r in R_LIST
    tag = "p$(p)_r$(r)"
    isfile(joinpath(DATA, "raw", "fig1_$(tag).jls")) && continue
    one_run(EPS, p, r; tag = tag)
end

println("FIG1 SWEEP DONE")

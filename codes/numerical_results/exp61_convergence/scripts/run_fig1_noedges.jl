# 6.1 Fig.-1 system — NO-EDGE-CORRECTION contrast sweep (mirrors the slab
# noedges archive in data/raw/noedges_run1/). Same geometry/protocol as
# run_fig1_sweep.jl, but with correct_edges = false. Output goes to
# data/raw/fig1_noedges/fig1_p<p>_r<r>.jls. Errors are computed at plot time
# against the existing edge-corrected reference data/raw/fig1_ref.jls (the
# converged solution), exactly as plot_slab.jl does for the slab.
#
# Run on a worker (e.g. worker7014); FINUFFT capped via OMP to dodge the
# high-thread FFTW-plan pathology:
#   GMRES_VERBOSE=0 JULIA_NUM_THREADS=96 OMP_NUM_THREADS=16 \
#     julia --project=codes/numerical_results \
#       codes/numerical_results/exp61_convergence/scripts/run_fig1_noedges.jl

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
using Printf

const SYS = system_fig1()
const DATA = joinpath(@__DIR__, "..", "data")
const OUTDIR = joinpath(DATA, "raw", "fig1_noedges")
const CSVPATH = joinpath(DATA, "sweep_fig1_noedges.csv")

const EPS = 1e-4
const P_LIST = [2, 4, 6]
const R_LIST = 1:5

const GMRES_VERBOSE = parse(Int, get(ENV, "GMRES_VERBOSE", "1"))

function one_run(eps, p, r; margin = 1.25, tag = "run")
    t_all = time()
    res = solve_system(SYS; eps = eps, p = p, r = r, src_margin = margin,
                       gmres_verbose = GMRES_VERBOSE, correct_edges = false)
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
    save_ref(joinpath(OUTDIR, "fig1_$(tag).jls"),
             (; sys = SYS.name, eps, p, r, margin, V,
                phi_near = phi[:near], phi_supp = phi[:supp], phi_far = phi[:far],
                N = res.N, niter = res.niter, residual = res.residual,
                history = res.gmres_history, times = res.times, wall,
                n_near_pairs = res.n_near_pairs, p_up_max = res.p_up_max,
                p_up_mean = res.p_up_mean, corr_bytes = res.corr_bytes))
    append_csv_row(CSVPATH, run_cols(res; V = V, wall = round(wall; digits = 2), tag = tag))
    @printf("[fig1-noedge] %s eps=%.0e p=%d r=%d N=%d niter=%d V=%.10e wall=%.1fs\n",
            tag, eps, p, r, res.N, res.niter, V, wall)
    flush(stdout)
end

println(">>> warm-up (compile)"); flush(stdout)
solve_system(slab_internal(); eps = 1e-2, p = 2, r = 1, correct_edges = false)

for p in P_LIST, r in R_LIST
    tag = "p$(p)_r$(r)"
    isfile(joinpath(OUTDIR, "fig1_$(tag).jls")) && continue
    one_run(EPS, p, r; tag = tag)
end

println("FIG1 NOEDGES SWEEP DONE")

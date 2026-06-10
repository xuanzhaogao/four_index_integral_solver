# 6.1 (revised per user 2026-06-10) — slab 10 x 10 x 1 (eps1 = 10), Gaussian
# source (s = 0.05) at the slab center, support fully inside the slab
# (truncation radius 0.24 at eps = 1e-4 < half-thickness 0.5).
#
# Fixed tied tolerance eps = 1e-4; sweep p in {2, 4, 6} x r in 1..5.
# Reference (differs in every discretization parameter): p = 8, eps = 1e-6,
# r = 7 (deepest test + 2), finer source/TKM grid (margin 1.4).
#
# Per run: V (target Gaussian at (0.2, 0, 0)), phi at the three fixed zone sets,
# DOF, N_iter, near-pair stats, phase timings. Raw values -> .jls; errors are
# computed at analysis time against the reference.

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
using Printf

const SYS = slab_internal()
const DATA = joinpath(@__DIR__, "..", "data")
const CSVPATH = joinpath(DATA, "sweep_slab.csv")

const EPS = 1e-4
const P_LIST = [2, 4, 6]
const R_LIST = 1:5
# Reference uses correct_edges=true (adaptive quadtree on edge-touching pairs):
# without it the edge-region quadrature error decays only ~2^-r and the
# reference would be the least-converged run of the suite (verified 2026-06-10).
const REF = (p = 8, eps = 1e-6, r = 6, margin = 1.4, edges = true)

const GMRES_VERBOSE = parse(Int, get(ENV, "GMRES_VERBOSE", "1"))

function one_run(eps, p, r; margin = 1.25, tag = "run", correct_edges = true)
    t_all = time()
    res = solve_system(SYS; eps = eps, p = p, r = r, src_margin = margin,
                       gmres_verbose = GMRES_VERBOSE, correct_edges = correct_edges)
    tdict = res.times
    # batched evaluation: V-target grid + all three zone sets in ONE operator
    # (single post-refinement + FMM + hcubature assembly), then split.
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
    save_ref(joinpath(DATA, "raw", "slab_$(tag).jls"),
             (; sys = SYS.name, eps, p, r, margin, V,
                phi_near = phi[:near], phi_supp = phi[:supp], phi_far = phi[:far],
                N = res.N, niter = res.niter, residual = res.residual,
                history = res.gmres_history, times = res.times, wall,
                n_near_pairs = res.n_near_pairs, p_up_max = res.p_up_max,
                p_up_mean = res.p_up_mean, corr_bytes = res.corr_bytes))
    append_csv_row(CSVPATH, run_cols(res; V = V, wall = round(wall; digits = 2), tag = tag))
    @printf("[slab] %s eps=%.0e p=%d r=%d N=%d niter=%d V=%.10e wall=%.1fs\n",
            tag, eps, p, r, res.N, res.niter, V, wall)
    flush(stdout)
end

println(">>> warm-up (compile)"); flush(stdout)
solve_system(slab_internal(); eps = 1e-2, p = 2, r = 1)

if !isfile(joinpath(DATA, "raw", "slab_ref.jls"))
    println(">>> reference: p=$(REF.p) eps=$(REF.eps) r=$(REF.r) edges=$(REF.edges)"); flush(stdout)
    one_run(REF.eps, REF.p, REF.r; margin = REF.margin, tag = "ref", correct_edges = REF.edges)
end

for p in P_LIST, r in R_LIST
    tag = "p$(p)_r$(r)"
    isfile(joinpath(DATA, "raw", "slab_$(tag).jls")) && continue
    one_run(EPS, p, r; tag = tag)
end

println("SLAB SWEEP DONE")

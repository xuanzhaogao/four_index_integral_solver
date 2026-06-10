# 6.4 Part A — solver scaling on the 6.2 slab problem (eps = 1e-6, p = 8).
# N grown by uniform dyadic refinement (x4 per level) beyond what adaptivity
# requires. Per level: setup (FMM plan, KD-tree pair search + upsampled block
# assembly), per-iteration FMM matvec and near-correction apply (median of
# reps), GMRES total, sparse-correction memory, p_up distribution.
# GMRES is skipped (and flagged) if the projected time exceeds TIME_CAP_S.

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using Krylov, LinearAlgebra, LinearMaps, SparseArrays, Statistics, Printf

const SYS = slab_system()
const DATA = joinpath(@__DIR__, "..", "data")
const CSVPATH = joinpath(DATA, "scaling.csv")
const EPS, P, R = 1e-6, 8, 4
const K_LIST = 0:3
const TIME_CAP_S = 3 * 3600.0
const N_CAP = 2.0e7

median_time(f, x; reps = 5) = begin
    f(x)  # warm
    ts = [(t0 = time(); f(x); time() - t0) for _ in 1:reps]
    median(ts)
end

println(">>> warm-up")
solve_system(system1(); eps = 1e-2, p = 4, r = 1)

println(">>> base adaptive solve")
base = solve_system(SYS; eps = EPS, p = P, r = R)
@printf("base N=%d niter=%d\n", base.N, base.niter)
expected_iters = base.niter

for k in K_LIST
    donefile = joinpath(DATA, "raw", "scaling_k$(k).jls")
    isfile(donefile) && (println("k=$k exists, skipping"); continue)

    iface = k == 0 ? base.interface : uniform_refine(base.interface, k)
    N = BI.num_points(iface)
    if N > N_CAP
        println("k=$k N=$N exceeds N_CAP=$(N_CAP): stopping (reported, not silently skipped)")
        break
    end
    @printf(">>> level k=%d  N=%d\n", k, N)

    t0 = time(); rhs = BI.rhs_dielectric_box3d_hybrid(iface, base.screened_vs, 1.0,
                  fmm_tol_of(EPS); tkm_kmax = base.tkm_kmax); t_rhs = time() - t0
    t0 = time(); D_base = BI.laplace3d_DT_fmm3d(iface, fmm_tol_of(EPS)); t_fmm_setup = time() - t0
    t0 = time(); nb = BI.build_neighbor_list(iface, 64, EPS); t_pairs = time() - t0
    t0 = time(); corr = BI.laplace3d_DT_corrections(iface, nb.upsample, nb.adaptive); t_blocks = time() - t0

    tvec = Harness._diag_coeffs(iface)
    t_matvec_fmm = median_time(x -> D_base * x, randn(N))
    t_matvec_corr = median_time(x -> corr * x, randn(N))
    apply = x -> (D_base * x) .+ (corr * x) .+ (tvec .* x)
    A = LinearMap{Float64}(apply, N, N)

    p_ups = collect(values(nb.upsample))
    projected = (t_matvec_fmm + t_matvec_corr) * expected_iters * 1.6
    gmres_ok = projected < TIME_CAP_S
    niter, t_gmres, resid = -1, NaN, NaN
    if gmres_ok
        t0 = time()
        sigma, stats = Krylov.gmres(A, rhs; atol = 1e-14, rtol = EPS, itmax = 500)
        t_gmres = time() - t0
        niter = stats.niter
        resid = norm(A * sigma - rhs) / norm(rhs)
        global expected_iters = niter
    else
        @printf("    GMRES SKIPPED: projected %.0fs > cap %.0fs\n", projected, TIME_CAP_S)
    end

    row = (k = k, N = N, n_panels = length(iface.panels),
           t_rhs = round(t_rhs; digits = 3), t_fmm_setup = round(t_fmm_setup; digits = 3),
           t_pair_search = round(t_pairs; digits = 3), t_block_assembly = round(t_blocks; digits = 3),
           t_matvec_fmm = round(t_matvec_fmm; digits = 4), t_matvec_corr = round(t_matvec_corr; digits = 4),
           t_gmres = round(t_gmres; digits = 2), niter = niter, gmres_run = gmres_ok,
           residual = resid, n_near_pairs = length(nb.upsample),
           p_up_max = isempty(p_ups) ? 0 : maximum(p_ups),
           p_up_mean = isempty(p_ups) ? 0.0 : round(mean(p_ups); digits = 3),
           corr_nnz = nnz(corr), corr_mb = round(Base.summarysize(corr) / 1e6; digits = 1),
           max_order = 64, eps = EPS, p = P,
           hostname = gethostname(), nthreads = Base.Threads.nthreads())
    append_csv_row(CSVPATH, row)
    save_ref(donefile, merge(row, (; p_up_values = p_ups)))
    @printf("    fmm=%.3fs corr=%.4fs pairs=%.1fs blocks=%.1fs gmres=%.1fs (niter=%d) corr_mb=%.1f\n",
            t_matvec_fmm, t_matvec_corr, t_pairs, t_blocks, t_gmres, niter,
            Base.summarysize(corr) / 1e6)
    flush(stdout)
end

# 6.4.3: p_up cap documentation comes from code inspection (laplace3d_near.jl:
# order clamped to max_order, no fallback for non-edge pairs) + the logged
# p_up_max hitting 64 if the cap binds.
println("SCALING DONE")

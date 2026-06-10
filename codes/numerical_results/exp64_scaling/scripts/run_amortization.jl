# 6.4 Part B — post-refinement amortization on the 6.2 slab problem.
# One solved sigma (eps = 1e-6, p = 8). M in {1, 10, 100, 1000} Gaussian target
# densities (s = 0.05, centers uniform inside Omega_m, seed 20260610), each on
# a 100^3 grid (~1e6 points before ball truncation).
#
#  (i)  no post-refinement: coarse-interface FMM + adaptive (hcubature) near
#       corrections per target. Timed on the first NSUB targets, linearly
#       extrapolated to larger M (extrapolated points are marked).
#  (ii) one-time post-refinement at h0 (+ prolongated sigma), then per target
#       FMM far field + residual hcubature near set.
# Validation: first 10 targets against strategy (i) with hcubature atol 1e-12.
# h0 sweep at M = 100 justifies the h0 choice.

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using LinearAlgebra, SparseArrays, Random, Printf, Statistics

const SYS = slab_system()
const DATA = joinpath(@__DIR__, "..", "data")
const EPS, P, R = 1e-6, 8, 4
const SEED_B = 20260610
const NGRID = 100
const TRUNC = 1e-10
const M_MAX = 1000
const M_MARKS = [1, 10, 100, 1000]
const NSUB = 10
const H0_DEFAULT = 0.025
const H0_LIST = [0.1, 0.05, 0.025, 0.0125, 0.00625]

println(">>> warm-up")
solve_system(system1(); eps = 1e-2, p = 4, r = 1)

println(">>> base solve (fixed sigma)")
res = solve_system(SYS; eps = EPS, p = P, r = R)
@printf("N=%d niter=%d\n", res.N, res.niter)

rng = MersenneTwister(SEED_B)
centers = [(0.5 * rand(rng) - 0.25, 0.5 * rand(rng) - 0.25, 0.35 + 0.1 * rand(rng)) for _ in 1:M_MAX]
target_grid(m) = gaussian_source(centers[m], SYS.src_sigma, 1e-9; n = NGRID)

"u_inc + given scattered values -> V_m; returns (V, t_inc)."
function v_of(res, vs_t, u_sc)
    t0 = time()
    u_inc = Harness.eval_incident(res, Matrix{Float64}(vs_t.positions))
    t_inc = time() - t0
    return dot(vs_t.weights .* vs_t.density, u_inc .+ u_sc), t_inc
end

# ---------------------------------------------------------------- strategy (ii)
"One-time post-refinement setup at threshold h0; returns refined data."
function setup_ii(res, h0)
    rt = sqrt(2 * log(1 / TRUNC)) * SYS.src_sigma
    lo = (-0.25 - rt, -0.25 - rt, 0.35 - rt)
    hi = (0.25 + rt, 0.25 + rt, 0.45 + rt)
    g = range(0, 1; length = 17)
    proxy = Matrix{Float64}(undef, 3, 17^3)
    c = 0
    for a in g, b in g, cc in g
        c += 1
        proxy[:, c] .= (lo[1] + a * (hi[1] - lo[1]), lo[2] + b * (hi[2] - lo[2]), lo[3] + cc * (hi[3] - lo[3]))
    end
    t0 = time()
    refined, parent_ids, from_split = BI._refine_interface_for_targets(res.interface, proxy, h0; range_factor = 5.0)
    prol = BI._refined_interface_prolongation(res.interface, refined, parent_ids, from_split)
    sig_ref = prol * res.sigma
    t_setup = time() - t0
    return (; refined, sig_ref, t_setup, n_refined = BI.num_points(refined))
end

"Per-target scattered eval on the pre-refined interface."
function eval_ii(setup, X)
    t0 = time()
    u = BI.laplace3d_pottrg_fmm3d(setup.refined, X, fmm_tol_of(EPS)) * setup.sig_ref
    t_fmm = time() - t0
    t0 = time()
    tnl = BI.build_target_neighbor_list(setup.refined, X, false; range_factor = 5.0)
    n_hcub = isempty(tnl) ? 0 : sum(length(v) for v in values(tnl))
    u .+= BI.laplace3d_pottrg_corrections_hcubature(setup.refined, X, tnl, EPS) * setup.sig_ref
    t_near = time() - t0
    return u, t_fmm, t_near, n_hcub
end

# ------------------------------------------------- validation truth (10 targets)
println(">>> validation truth: strategy (i) at hcubature atol 1e-12, 10 targets")
truth = Float64[]
truth_file = joinpath(DATA, "raw", "amort_truth.jls")
if isfile(truth_file)
    truth = load_ref(truth_file).truth
else
    for m in 1:NSUB
        vs_t = target_grid(m)
        u_sc, _ = eval_scatter_with_h0(res, Matrix{Float64}(vs_t.positions), Inf; hcub_atol = 1e-12)
        V, _ = v_of(res, vs_t, u_sc)
        push!(truth, V)
        @printf("  truth m=%d V=%.12e\n", m, V); flush(stdout)
    end
    save_ref(truth_file, (; truth, centers = centers[1:NSUB]))
end

# ------------------------------------------------------------- strategy (i) subset
println(">>> strategy (i), $NSUB targets")
ti_rows = NamedTuple[]
for m in 1:NSUB
    vs_t = target_grid(m)
    X = Matrix{Float64}(vs_t.positions)
    t0 = time()
    u_sc, st = eval_scatter_with_h0(res, X, Inf; hcub_atol = EPS)
    t_sc = time() - t0
    V, t_inc = v_of(res, vs_t, u_sc)
    err = abs(V - truth[m]) / abs(truth[m])
    push!(ti_rows, (; m, V, t_sc, t_inc, n_hcub = st.n_hcub, err))
    @printf("  (i) m=%d t_sc=%.1fs t_inc=%.1fs n_hcub=%d err=%.2e\n", m, t_sc, t_inc, st.n_hcub, err)
    flush(stdout)
end
save_ref(joinpath(DATA, "raw", "amort_strategy_i.jls"), (; rows = ti_rows))
t_i_per_target = mean(r.t_sc + r.t_inc for r in ti_rows)
for M in M_MARKS
    append_csv_row(joinpath(DATA, "amortization.csv"), (
        strategy = "i", M = M, h0 = Inf,
        t_setup = 0.0,
        t_total = round(M <= NSUB ? sum(r.t_sc + r.t_inc for r in ti_rows[1:M]) : t_i_per_target * M; digits = 1),
        extrapolated = M > NSUB,
        err_max = round(maximum(r.err for r in ti_rows); sigdigits = 3)))
end

# ------------------------------------------------------- strategy (ii) full pass
println(">>> strategy (ii), h0=$H0_DEFAULT, cumulative pass to M=$M_MAX")
setup = setup_ii(res, H0_DEFAULT)
@printf("  setup: %.1fs  (N %d -> %d)\n", setup.t_setup, res.N, setup.n_refined)
cum = 0.0
errs_ii = Float64[]
csvpath = joinpath(DATA, "amortization.csv")
for m in 1:M_MAX
    vs_t = target_grid(m)
    X = Matrix{Float64}(vs_t.positions)
    u_sc, t_fmm, t_near, n_hcub = eval_ii(setup, X)
    V, t_inc = v_of(res, vs_t, u_sc)
    global cum += t_fmm + t_near + t_inc
    m <= NSUB && push!(errs_ii, abs(V - truth[m]) / abs(truth[m]))
    if m in M_MARKS
        append_csv_row(csvpath, (
            strategy = "ii", M = m, h0 = H0_DEFAULT,
            t_setup = round(setup.t_setup; digits = 1),
            t_total = round(setup.t_setup + cum; digits = 1),
            extrapolated = false,
            err_max = round(maximum(errs_ii); sigdigits = 3)))
        @printf("  (ii) M=%d cumulative=%.1fs (+setup %.1fs)\n", m, cum, setup.t_setup)
        flush(stdout)
    end
end
save_ref(joinpath(DATA, "raw", "amort_strategy_ii.jls"),
         (; h0 = H0_DEFAULT, t_setup = setup.t_setup, n_refined = setup.n_refined, errs_ii))

# ------------------------------------------------------------------- h0 sweep
println(">>> h0 sweep at M=100")
for h0 in H0_LIST
    donefile = joinpath(DATA, "raw", "amort_h0_$(h0).jls")
    isfile(donefile) && continue
    s = setup_ii(res, h0)
    cumh = 0.0
    errs = Float64[]
    for m in 1:100
        vs_t = target_grid(m)
        X = Matrix{Float64}(vs_t.positions)
        u_sc, t_fmm, t_near, _ = eval_ii(s, X)
        V, t_inc = v_of(res, vs_t, u_sc)
        cumh += t_fmm + t_near + t_inc
        m <= NSUB && push!(errs, abs(V - truth[m]) / abs(truth[m]))
    end
    append_csv_row(joinpath(DATA, "amort_h0.csv"), (
        h0 = h0, M = 100, n_refined = s.n_refined,
        t_setup = round(s.t_setup; digits = 1), t_total = round(s.t_setup + cumh; digits = 1),
        err_max = round(maximum(errs); sigdigits = 3)))
    save_ref(donefile, (; h0, t_setup = s.t_setup, t_total = s.t_setup + cumh, errs))
    @printf("  h0=%g: setup=%.1fs total=%.1fs err=%.2e\n", h0, s.t_setup, s.t_setup + cumh, maximum(errs))
    flush(stdout)
end

println("AMORTIZATION DONE")

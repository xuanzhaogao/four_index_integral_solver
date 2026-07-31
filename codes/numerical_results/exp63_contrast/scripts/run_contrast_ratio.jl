# 6.3 — dielectric contrast on the Fig.-1 geometry (10x10x1 slab on two
# 10x10x10 cubes, source inside the slab): per-contrast convergence study.
#
# eps1 = 4 and eps_slab = 10 held fixed; sweep the second substrate
#   eps2 in {2, 6, 20, 60, 200}   (log-spaced 2 -> 200)
#   => buried-face contrast gamma_12 = (eps2-eps1)/(eps2+eps1)
#      in {-0.333, 0.2, 0.667, 0.875, 0.961}   (eps2 = 2 inverts the sign).
# For EACH eps2: p = 6 fixed, edge depth r swept 1..5; per-eps2 reference at
# p = 8, eps = 1e-6, r = 6, margin 1.4 (the 6.1 protocol). Recorded per run:
# V (accuracy vs the eps2's own reference), N_iter, DOF, GMRES history, times.
# eps2 = 4 or 10 are EXCLUDED by construction (they would degenerate the
# cube|cube resp. slab|cube_2 face; remap_eps errors on that).
#
# Mesh protocol — ONE mesh per (p, r) protocol, eps-remapped per eps2
# (Harness.remap_eps): the screened source lives entirely in the slab
# (eps_slab unchanged), so the RHS-adaptive mesh is IDENTICAL for every eps2;
# building it once per r and remapping the per-panel permittivities gives
# exactly the production mesh while isolating the contrast effect at fixed
# discretization. gmres_atol = 0: pure relative stopping, uniform everywhere.
#
# Cost estimate (96 threads): tests ~1 h (25 runs, r = 1..5 x 5 eps2);
# references ~20-35 min each x 5  =>  ~2.5-3 h total.
# Resumable per case and per variant.
#
# Run:    GMRES_VERBOSE=0 JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 \
#           julia --project=. exp63_contrast/scripts/run_contrast_ratio.jl
# Smoke:  CONTRAST_SMOKE=1 julia --project=. exp63_contrast/scripts/run_contrast_ratio.jl
#         (tiny parameters, validates the remap/override path end to end)

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using Printf

const SMOKE = get(ENV, "CONTRAST_SMOKE", "0") == "1"
const GMRES_VERBOSE = parse(Int, get(ENV, "GMRES_VERBOSE", "0"))

# RERUN_TAG: appends a suffix to this experiment's output directory so a rerun
# never overwrites the data behind the submitted manuscript. Empty = original paths.
const TAG = get(ENV, "RERUN_TAG", "")
const DATA = joinpath(@__DIR__, "..", "data" * TAG, SMOKE ? "smoke" : "")
mkpath(joinpath(DATA, "raw"))
const CSVPATH = joinpath(DATA, "contrast_ratio.csv")
mkpath(joinpath(DATA, "raw"))

const E2_LIST = SMOKE ? [2.0, 200.0] : [2.0, 6.0, 20.0, 60.0, 200.0]
const R_LIST = SMOKE ? (1:2) : (1:5)
const TST = SMOKE ? (p = 2, eps = 1e-2, margin = 1.25) :
                    (p = 6, eps = 1e-4, margin = 1.25)
const REF = SMOKE ? (p = 4, eps = 1e-3, r = 3, margin = 1.25) :
                    (p = 8, eps = 1e-6, r = 6, margin = 1.4)
const VAC_EPS = SMOKE ? 1e-6 : 1e-13

const SYS1 = system_fig1()                      # baseline (4, 12, 10)
sys_of(e2) = system_fig1(eps2 = e2)
mapping_of(e2) = Dict(12.0 => e2)
gamma12(e2) = (e2 - 4.0) / (e2 + 4.0)

# baseline mesh exactly as solve_system would build it for protocol (p, eps, margin) at depth r
function build_mesh(p, r, eps, margin)
    vs = gaussian_source(SYS1.src_center, SYS1.src_sigma, eps; margin = margin)
    svs = Harness.screened_source(SYS1, vs)
    kmax = BI._estimate_tkm3dc_kmax(BI._estimate_source_spacing(svs))
    t0 = time()
    iface = BI.multi_dielectric_box3d_rhs_adaptive(
        p, Harness.l_ec_of(SYS1, r), SYS1.boxes, SYS1.epses, svs,
        1.0, eps, SYS1.eps_out; max_depth = 128, tkm_kmax = kmax)
    @printf(">>> mesh p=%d r=%d eps=%.0e: %d panels, %d points (%.0fs)\n",
            p, r, eps, length(iface.panels), BI.num_points(iface), time() - t0)
    flush(stdout)
    return (; iface, kmax)
end

const MESH = Dict{Tuple{Symbol, Int}, Any}()
getmesh(kind, p, r, eps, margin) =
    get!(() -> build_mesh(p, r, eps, margin), MESH, (kind, r))

function run_case(e2, p, r, eps, margin, mesh, variant)
    syse = sys_of(e2)
    vs = gaussian_source(syse.src_center, syse.src_sigma, eps; margin = margin)
    svs = Harness.screened_source(syse, vs)
    ife = remap_eps(mesh.iface, mapping_of(e2))
    t0 = time()
    res = solve_system(syse; eps = eps, p = p, r = r, src_margin = margin,
                       interface_override = ife, screened_vs_override = svs,
                       tkm_kmax_override = mesh.kmax,
                       gmres_atol = 0.0, gmres_verbose = GMRES_VERBOSE)
    V = eval_V(res; t_out = res.times, margin = margin)
    wall = time() - t0
    append_csv_row(CSVPATH, run_cols(res; V = V, variant = variant, eps2 = e2,
                   gamma12 = gamma12(e2), wall = round(wall; digits = 2)))
    @printf("[6.3] %-4s eps2=%-5g r=%d gamma12=%+.4f N=%d niter=%d V=%.12e (%.0fs)\n",
            variant, e2, r, gamma12(e2), res.N, res.niter, V, wall)
    flush(stdout)
    return (; eps2 = e2, r, V, N = res.N, niter = res.niter, residual = res.residual,
            history = copy(res.gmres_history), times = copy(res.times), wall,
            n_near = res.n_near_pairs, n_adaptive = res.n_adaptive_pairs,
            p_up_max = res.p_up_max, corr_bytes = res.corr_bytes)
end

println(">>> warm-up"); flush(stdout)
solve_system(slab_internal(); eps = 1e-2, p = 2, r = 1)

# vacuum interaction (no interfaces) — contrast-independent, computed once
let f = joinpath(DATA, "raw", "v_vacuum.jls")
    if !isfile(f)
        Vvac = vacuum_V(SYS1, VAC_EPS)
        save_ref(f, (; Vvac))
        @printf("V_vacuum = %.12e\n", Vvac); flush(stdout)
    end
end

println(">>> phase 1: tests (p=$(TST.p), eps=$(TST.eps), r=$(first(R_LIST))..$(last(R_LIST)))"); flush(stdout)
for r in R_LIST, e2 in E2_LIST
    donefile = joinpath(DATA, "raw", "ratio_test_e2_$(e2)_r$(r).jls")
    isfile(donefile) && (println("test eps2=$e2 r=$r exists, skipping"); continue)
    out = run_case(e2, TST.p, r, TST.eps, TST.margin,
                   getmesh(:test, TST.p, r, TST.eps, TST.margin), "test")
    save_ref(donefile, out)
    GC.gc()
end

println(">>> phase 2: per-eps2 references (p=$(REF.p), eps=$(REF.eps), r=$(REF.r))"); flush(stdout)
for e2 in E2_LIST
    donefile = joinpath(DATA, "raw", "ratio_ref_e2_$(e2).jls")
    isfile(donefile) && (println("ref eps2=$e2 exists, skipping"); continue)
    out = run_case(e2, REF.p, REF.r, REF.eps, REF.margin,
                   getmesh(:ref, REF.p, REF.r, REF.eps, REF.margin), "ref")
    save_ref(donefile, out)
    GC.gc()
end

println("RATIO SWEEP DONE")

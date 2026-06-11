# 6.3 — dielectric contrast on the Fig.-1 geometry (10x10x1 slab on two
# 10x10x10 cubes, source inside the slab).
#
# eps1 = 4 and eps_slab = 10 held fixed; sweep the second substrate
#   eps2 in {6, 12, 40, 400, 4000, 1e11}
#   => buried-face contrast gamma_12 = (eps2-eps1)/(eps2+eps1)
#      in {0.2, 0.5, 0.818, 0.980, 0.998, 1 - 8e-11}.
# This changes the cube|cube, slab|cube_2 and cube_2|vacuum contrasts jointly
# while the source-side faces stay fixed — coexisting gammas at every point.
# eps2 = 1e11 is the conductor-limit anchor (cube 2 effectively a floating
# conductor). eps2 = 4 or 10 are EXCLUDED: they would degenerate the
# cube|cube resp. slab|cube_2 face (eps_in == eps_out; remap_eps errors).
#
# Experimental control — ONE mesh per protocol, built at the eps2 = 12
# baseline and eps-remapped per case (Harness.remap_eps). The screened source
# lives entirely in the slab (eps_slab unchanged), so the baseline mesh is the
# EXACT RHS-adaptive mesh for every eps2; only per-panel permittivities
# change. This isolates the conditioning effect of contrast at fixed
# discretization (per-case rebuilds would also drift because the adaptive
# thresholds are absolute). gmres_atol = 0: pure relative stopping, uniform
# across the sweep.
#
# Outputs per eps2: V test + per-case reference (p=8, eps=1e-6, r=6, margin
# 1.4 — the 6.1 protocol), N_iter, full GMRES residual histories; V_vacuum
# once (physics curve V/V_vac vs eps2). The eps2 = 12 test/ref reproduce the
# 6.1 fig1 p6_r4 run and reference (V_ref = 0.055678122198) — analysis
# cross-check.
#
# Cost estimate (96 threads): tests ~3-5 min each, references ~20-35 min
# each, x 6  =>  ~2.5-3 h total. Resumable per case and per variant.
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

const DATA = joinpath(@__DIR__, "..", "data", SMOKE ? "smoke" : "")
const CSVPATH = joinpath(DATA, "contrast_ratio.csv")
mkpath(joinpath(DATA, "raw"))

const E2_LIST = SMOKE ? [6.0, 1.0e11] : [6.0, 12.0, 40.0, 400.0, 4000.0, 1.0e11]
const TST = SMOKE ? (p = 2, eps = 1e-2, r = 1, margin = 1.25) :
                    (p = 6, eps = 1e-4, r = 4, margin = 1.25)
const REF = SMOKE ? (p = 4, eps = 1e-3, r = 2, margin = 1.25) :
                    (p = 8, eps = 1e-6, r = 6, margin = 1.4)
const VAC_EPS = SMOKE ? 1e-6 : 1e-13

const SYS1 = system_fig1()                      # baseline (4, 12, 10)
sys_of(e2) = system_fig1(eps2 = e2)
mapping_of(e2) = Dict(12.0 => e2)
gamma12(e2) = (e2 - 4.0) / (e2 + 4.0)

# baseline mesh exactly as solve_system would build it for this protocol
function build_mesh(prot)
    vs = gaussian_source(SYS1.src_center, SYS1.src_sigma, prot.eps; margin = prot.margin)
    svs = Harness.screened_source(SYS1, vs)
    kmax = BI._estimate_tkm3dc_kmax(BI._estimate_source_spacing(svs))
    t0 = time()
    iface = BI.multi_dielectric_box3d_rhs_adaptive(
        prot.p, Harness.l_ec_of(SYS1, prot.r), SYS1.boxes, SYS1.epses, svs,
        1.0, prot.eps, SYS1.eps_out; max_depth = 128, tkm_kmax = kmax)
    @printf(">>> mesh p=%d r=%d eps=%.0e: %d panels, %d points (%.0fs)\n",
            prot.p, prot.r, prot.eps, length(iface.panels), BI.num_points(iface), time() - t0)
    flush(stdout)
    return (; iface, kmax)
end

const MESH = Dict{Symbol, Any}()
getmesh(key, prot) = get!(() -> build_mesh(prot), MESH, key)

function run_case(e2, prot, mesh, variant)
    syse = sys_of(e2)
    vs = gaussian_source(syse.src_center, syse.src_sigma, prot.eps; margin = prot.margin)
    svs = Harness.screened_source(syse, vs)
    ife = remap_eps(mesh.iface, mapping_of(e2))
    t0 = time()
    res = solve_system(syse; eps = prot.eps, p = prot.p, r = prot.r,
                       src_margin = prot.margin,
                       interface_override = ife, screened_vs_override = svs,
                       tkm_kmax_override = mesh.kmax,
                       gmres_atol = 0.0, gmres_verbose = GMRES_VERBOSE)
    V = eval_V(res; t_out = res.times, margin = prot.margin)
    wall = time() - t0
    append_csv_row(CSVPATH, run_cols(res; V = V, variant = variant, eps2 = e2,
                   gamma12 = gamma12(e2), wall = round(wall; digits = 2)))
    @printf("[6.3] %-4s eps2=%-8s gamma12=%.6f N=%d niter=%d V=%.12e (%.0fs)\n",
            variant, string(e2), gamma12(e2), res.N, res.niter, V, wall)
    flush(stdout)
    return (; eps2 = e2, V, N = res.N, niter = res.niter, residual = res.residual,
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

println(">>> phase 1: tests (p=$(TST.p), eps=$(TST.eps), r=$(TST.r))"); flush(stdout)
for e2 in E2_LIST
    donefile = joinpath(DATA, "raw", "ratio_test_e2_$(e2).jls")
    isfile(donefile) && (println("test eps2=$e2 exists, skipping"); continue)
    out = run_case(e2, TST, getmesh(:test, TST), "test")
    save_ref(donefile, out)
    GC.gc()
end

println(">>> phase 2: per-case references (p=$(REF.p), eps=$(REF.eps), r=$(REF.r))"); flush(stdout)
for e2 in E2_LIST
    donefile = joinpath(DATA, "raw", "ratio_ref_e2_$(e2).jls")
    isfile(donefile) && (println("ref eps2=$e2 exists, skipping"); continue)
    out = run_case(e2, REF, getmesh(:ref, REF), "ref")
    save_ref(donefile, out)
    GC.gc()
end

println("RATIO SWEEP DONE")

# 6.3.1 / 6.3.2 — dielectric contrast on the Fig.-1 geometry (10x10x1 slab on
# two 10x10x10 cubes, source inside the slab).
#
# Joint permittivity scaling: (eps1, eps2, eps_slab) = m * (4, 12, 10), vacuum
# fixed at 1, m in {0.5, 1, 2, 8, 32, 100, 1e4, 1e11}. External-face contrasts
# approach the conductor limit (gamma -> 1) as m grows while the INTERNAL
# contrasts (cube|cube, slab|cube) stay fixed — multiple coexisting gammas at
# every m. m = 1e11 is the conductor-limit anchor (1 - gamma_ext ~ 2e-12).
#
# Experimental control — ONE mesh per protocol, built at m = 1 and
# eps-remapped per m (Harness.remap_eps). This is required, not just cheaper:
# the RHS-adaptive thresholds are absolute and the screened source is exactly
# svs(m=1)/m (Gaussian support fully inside the slab), so per-m mesh rebuilds
# would coarsen with contrast and conflate discretization with conditioning.
# (At m = 0.5 the RHS is 2x larger than the mesh was built for — well inside
# one dyadic refinement step of granularity.)
#
# gmres_atol = 0 everywhere (pure relative stopping): at m = 1e11 the RHS norm
# is ~1e-11 and the default 1e-14 absolute floor would contaminate N_iter.
#
# Outputs per m: V test + per-contrast reference (p=8, eps=1e-6, r=6, margin
# 1.4 — the 6.1 protocol), N_iter, full GMRES residual histories; V_vacuum
# once (screening curve V/V_vac: ~1/(10 m) bulk decay, then saturation at the
# compound-conductor capacitance value). The m = 1 test/ref reproduce the 6.1
# fig1 p6_r4 run and reference (V_ref = 0.055678122198) — analysis cross-check.
#
# Cost estimate (96 threads): tests ~30-45 min total; references dominate at
# ~20-35 min each x 8  =>  ~4-5 h total. Resumable per m and per variant.
#
# Run:    GMRES_VERBOSE=0 JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 \
#           julia --project=. exp63_contrast/scripts/run_contrast.jl
# Smoke:  CONTRAST_SMOKE=1 julia --project=. exp63_contrast/scripts/run_contrast.jl
#         (tiny parameters, validates the remap/override path end to end)

include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
import BoundaryIntegral as BI
using Printf

const SMOKE = get(ENV, "CONTRAST_SMOKE", "0") == "1"
const GMRES_VERBOSE = parse(Int, get(ENV, "GMRES_VERBOSE", "0"))

const DATA = joinpath(@__DIR__, "..", "data", SMOKE ? "smoke" : "")
const CSVPATH = joinpath(DATA, "contrast.csv")
mkpath(joinpath(DATA, "raw"))

const M_LIST = SMOKE ? [1.0, 100.0, 1.0e11] :
                       [0.5, 1.0, 2.0, 8.0, 32.0, 100.0, 1.0e4, 1.0e11]
const TST = SMOKE ? (p = 2, eps = 1e-2, r = 1, margin = 1.25) :
                    (p = 6, eps = 1e-4, r = 4, margin = 1.25)
const REF = SMOKE ? (p = 4, eps = 1e-3, r = 2, margin = 1.25) :
                    (p = 8, eps = 1e-6, r = 6, margin = 1.4)
const VAC_EPS = SMOKE ? 1e-6 : 1e-13

const SYS1 = system_fig1()                      # m = 1 baseline (4, 12, 10)
sys_of(m) = system_fig1(eps1 = 4.0m, eps2 = 12.0m, eps_slab = 10.0m)
mapping_of(m) = Dict(4.0 => 4.0m, 12.0 => 12.0m, 10.0 => 10.0m)
gamma_slab(m) = (10.0m - 1.0) / (10.0m + 1.0)   # slab|vacuum external contrast

# m = 1 mesh exactly as solve_system would build it for this protocol
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

function run_case(m, prot, mesh, variant)
    sysm = sys_of(m)
    vs = gaussian_source(sysm.src_center, sysm.src_sigma, prot.eps; margin = prot.margin)
    svs = Harness.screened_source(sysm, vs)
    ifm = remap_eps(mesh.iface, mapping_of(m))
    t0 = time()
    res = solve_system(sysm; eps = prot.eps, p = prot.p, r = prot.r,
                       src_margin = prot.margin,
                       interface_override = ifm, screened_vs_override = svs,
                       tkm_kmax_override = mesh.kmax,
                       gmres_atol = 0.0, gmres_verbose = GMRES_VERBOSE)
    V = eval_V(res; t_out = res.times, margin = prot.margin)
    wall = time() - t0
    append_csv_row(CSVPATH, run_cols(res; V = V, variant = variant, m = m,
                   gamma_slab = gamma_slab(m), wall = round(wall; digits = 2)))
    @printf("[6.3] %-4s m=%-7s gamma=%.6f N=%d niter=%d V=%.12e (%.0fs)\n",
            variant, string(m), gamma_slab(m), res.N, res.niter, V, wall)
    flush(stdout)
    return (; m, V, N = res.N, niter = res.niter, residual = res.residual,
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
for m in M_LIST
    donefile = joinpath(DATA, "raw", "contrast_test_m$(m).jls")
    isfile(donefile) && (println("test m=$m exists, skipping"); continue)
    out = run_case(m, TST, getmesh(:test, TST), "test")
    save_ref(donefile, out)
    GC.gc()
end

println(">>> phase 2: per-contrast references (p=$(REF.p), eps=$(REF.eps), r=$(REF.r))"); flush(stdout)
for m in M_LIST
    donefile = joinpath(DATA, "raw", "contrast_ref_m$(m).jls")
    isfile(donefile) && (println("ref m=$m exists, skipping"); continue)
    out = run_case(m, REF, getmesh(:ref, REF), "ref")
    save_ref(donefile, out)
    GC.gc()
end

println("CONTRAST DONE")

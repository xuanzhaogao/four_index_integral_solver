# 6.5 — walltime + peak-RAM benchmark of the production pipeline on the REAL
# system: graphene monolayer slab (L = 90, Lz = 3.35, eps_in = 3.5 calibrated)
# with the real Wannier pz orbital (k_323201_nb_144_c_15), single RHS.
#
# Source = |phi_1|^2 (squared pz orbital, centered, z = Lz/2), one solve, one
# evaluation (onsite channel: target density = source density). Mirrors
# ScreenedOrbitalSolve.solve_screened_mode stage by stage (same BI calls, same
# calibrated parameters as screened_monolayer_hund_calibrated.jl), inlined so
# each stage gets a walltime AND a peak-RSS snapshot.
#
# Phase grouping reported:
#   load       XSF read + VolumeSource truncation (one-time, I/O + cold)
#   precompute screened source + TKM kmax + RHS-adaptive interface + LHS
#              assembly (near-correction upsampling) — reusable across RHS
#   solve      RHS assembly (FMM) + GMRES
#   eval       u_int (TKM volume potential) + u_scatter (FMM + hcubature
#              corrected target eval) + reductions
#
# RAM: Sys.maxrss() is the process high-water mark (monotone); the per-phase
# increment attributes NEW peak memory to that phase. Base.gc_live_bytes()
# is also recorded. A warm-up solve (tiny Gaussian, same code path) runs
# first so timings are warm; the post-warm-up baseline RSS is recorded.
# Wrap the launch in `/usr/bin/time -v` for an OS-level peak-RSS cross-check.
#
# Env knobs: CORRECT_EDGES (default 1; 6.1 showed edge correction is
# load-bearing — set 0 to benchmark the legacy uncorrected-edges config),
# ORBBENCH_SMOKE=1 (coarse tolerances, same geometry + real orbital).
#
# Run:   JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 /usr/bin/time -v \
#          julia --project=. exp65_orbital_bench/scripts/run_single_rhs.jl
# Smoke: ORBBENCH_SMOKE=1 julia --project=. exp65_orbital_bench/scripts/run_single_rhs.jl

using BoundaryIntegral
import BoundaryIntegral as BI
using Krylov
using LinearAlgebra
using Printf
using Serialization

const GRAPHENE = normpath(joinpath(@__DIR__, "..", "..", "..", "graphene"))
include(joinpath(GRAPHENE, "bilayer_slab", "src", "ScreenedOrbitalSolve.jl"))
using .ScreenedOrbitalSolve
include(joinpath(GRAPHENE, "monolayer", "src", "MonolayerOrbitalLoader.jl"))
using .MonolayerOrbitalLoader

const SMOKE = get(ENV, "ORBBENCH_SMOKE", "0") == "1"
const CORRECT_EDGES = get(ENV, "CORRECT_EDGES", "1") == "1"

const DATA = joinpath(@__DIR__, "..", "data", SMOKE ? "smoke" : "")
mkpath(joinpath(DATA, "raw"))

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")

# calibrated production parameters (screened_monolayer_hund_calibrated.jl)
const LZ, Z_CENTER, EPS_IN, EPS_OUT, L = 3.35, 3.35 / 2, 3.5, 1.0, 90.0
const MODE = BI.SharpScreening()
const P = SMOKE ?
    (source_tol = 1e-2, n_quad = 6, edge_level = 2, rhs_tol = 1e-2,
     lhs_tol = 1e-3, gmres_atol = 1e-3, gmres_rtol = 1e-3,
     max_order = 64, max_depth = 12) :
    (source_tol = 1e-3, n_quad = 6, edge_level = 4, rhs_tol = 1e-3,
     lhs_tol = 1e-5, gmres_atol = 1e-5, gmres_rtol = 1e-5,
     max_order = 64, max_depth = 12)
const VOLUME_TOL = P.rhs_tol

gb(x) = x / 2^30
rss_gb() = gb(Sys.maxrss())
live_gb() = gb(Base.gc_live_bytes())

const STAGES = NamedTuple[]
function stage!(label, t)
    push!(STAGES, (; label, t, maxrss_gb = rss_gb(), live_gb = live_gb()))
    @printf("  %-38s %9.2f s   maxrss %6.2f GB   live %6.2f GB\n",
            label, t, rss_gb(), live_gb())
    flush(stdout)
end

# one full pass through every pipeline stage (used for tiny warm-up too)
function pipeline(vs, Lx, Ly, Lz, prm; record::Bool)
    note = record ? stage! : (l, t) -> nothing
    l_ec = Lz / 2.0^prm.edge_level * 1.01

    t = @elapsed screened_vs = BI.screened_volume_source(
        Lx, Ly, Lz, vs, EPS_IN, EPS_OUT, MODE; tol = prm.rhs_tol)
    note("screened_volume_source", t)
    t = @elapsed kmax = BI._estimate_tkm3dc_kmax(screened_vs)
    note("tkm kmax estimate", t)
    t = @elapsed interface = BI.single_dielectric_box3d_rhs_adaptive(
        Lx, Ly, Lz, prm.n_quad, screened_vs, 1.0, l_ec, prm.rhs_tol,
        EPS_IN, EPS_OUT, Float64; max_depth = prm.max_depth, tkm_kmax = kmax)
    note("interface build (RHS-adaptive)", t)
    t = @elapsed lhs = BI.Lhs_dielectric_box3d_fmm3d_corrected(
        interface, prm.lhs_tol, prm.lhs_tol, prm.max_order;
        correct_edges = CORRECT_EDGES)
    note("LHS assembly (near corrections)", t)

    t = @elapsed rhs = BI.Rhs_dielectric_box3d_fmm3d(interface, screened_vs, 1.0, prm.rhs_tol)
    note("RHS assembly (FMM)", t)
    local sigma, stats
    t = @elapsed begin
        sigma, stats = Krylov.gmres(lhs, rhs; atol = prm.gmres_atol,
                                    rtol = prm.gmres_rtol, verbose = 0)
    end
    note("GMRES solve", t)
    residual = norm(lhs * sigma - rhs) / max(norm(rhs), eps(Float64))

    # onsite evaluation: target density = source density
    tmat = Matrix{Float64}(vs.positions)
    tw = vs.weights .* vs.density
    t = @elapsed u_int = ScreenedOrbitalSolve.evaluate_volume_potential(
        screened_vs, tmat; tol = VOLUME_TOL, kmax = kmax)
    note("eval: volume potential (TKM)", t)
    local u_scatter
    t = @elapsed begin
        op = BI.laplace3d_pottrg_fmm3d_corrected_hcubature(
            interface, tmat, prm.lhs_tol, prm.lhs_tol, 5.0)
        u_scatter = op * sigma
    end
    note("eval: scattered potential (FMM+hcub)", t)

    u_int_raw = dot(u_int, tw)
    u_scatter_raw = dot(u_scatter, tw)
    return (; n_points = BI.num_points(interface), n_src = length(screened_vs.density),
            n_targets = size(tmat, 2), niter = stats.niter, residual,
            u_int_raw, u_scatter_raw, u_total_raw = u_int_raw + u_scatter_raw,
            tkm_kmax = kmax)
end

@printf("threads = %d   correct_edges = %s   smoke = %s\n",
        Threads.nthreads(), CORRECT_EDGES, SMOKE)

println(">>> warm-up (tiny Gaussian through every stage)"); flush(stdout)
let g = BI.GaussianVolumeSource((0.0, 0.0, 0.5), 0.3, 8, 1e-3)
    tw = @elapsed pipeline(g, 5.0, 5.0, 1.0,
        (; P..., edge_level = 2); record = false)
    @printf("  warm-up: %.1f s   baseline maxrss %.2f GB\n", tw, rss_gb())
end
const RSS_BASELINE = rss_gb()

println(">>> load: real pz orbital |phi_1|^2  ($(basename(XSF_1)))"); flush(stdout)
t_load = @elapsed begin
    dg = load_squared_xsf(XSF_1)
    sh = MonolayerOrbitalLoader._centering_shift(dg; tol = P.source_tol)
    dg = shift_datagrid(dg, (sh[1], sh[2], sh[3] + Z_CENTER))
    vs1 = BI.VolumeSource(dg, tol = P.source_tol)
    Nphi1 = sum(vs1.weights .* vs1.density)
end
stage!("load (XSF read + truncation)", t_load)
@printf("  Nphi1 = %.6f   source points = %d\n", Nphi1, length(vs1.density))

println(">>> benchmark run  L=$L  Lz=$LZ  eps_in=$EPS_IN  (warm)"); flush(stdout)
t_total = @elapsed res = pipeline(vs1, L, L, LZ, P; record = true)

u_total_ev = to_eV(res.u_total_raw, Nphi1, Nphi1)

# phase grouping
tof(lbl) = first(s.t for s in STAGES if s.label == lbl)
t_precompute = tof("screened_volume_source") + tof("tkm kmax estimate") +
               tof("interface build (RHS-adaptive)") + tof("LHS assembly (near corrections)")
t_solve = tof("RHS assembly (FMM)") + tof("GMRES solve")
t_eval = tof("eval: volume potential (TKM)") + tof("eval: scattered potential (FMM+hcub)")

println("=" ^ 72)
@printf("  interface points %d   src points %d   targets %d   niter %d   residual %.3e\n",
        res.n_points, res.n_src, res.n_targets, res.niter, res.residual)
@printf("  u_onsite = %.6e raw = %.6f eV\n", res.u_total_raw, u_total_ev)
println("-" ^ 72)
@printf("  %-12s %9.2f s\n", "load", t_load)
@printf("  %-12s %9.2f s\n", "precompute", t_precompute)
@printf("  %-12s %9.2f s\n", "solve", t_solve)
@printf("  %-12s %9.2f s\n", "eval", t_eval)
@printf("  %-12s %9.2f s   (sum of stages; total wall %.2f s)\n",
        "TOTAL", t_load + t_precompute + t_solve + t_eval, t_load + t_total)
@printf("  peak RSS %.2f GB   (baseline after warm-up %.2f GB)\n", rss_gb(), RSS_BASELINE)
println("=" ^ 72)

out = (; smoke = SMOKE, correct_edges = CORRECT_EDGES,
       L, Lz = LZ, eps_in = EPS_IN, eps_out = EPS_OUT, z_center = Z_CENTER,
       params = P, volume_tol = VOLUME_TOL,
       Nphi1, n_points = res.n_points, n_src = res.n_src, n_targets = res.n_targets,
       niter = res.niter, residual = res.residual, tkm_kmax = res.tkm_kmax,
       u_int_raw = res.u_int_raw, u_scatter_raw = res.u_scatter_raw,
       u_total_raw = res.u_total_raw, u_total_ev,
       t_load, t_precompute, t_solve, t_eval,
       stages = copy(STAGES), rss_baseline_gb = RSS_BASELINE, rss_peak_gb = rss_gb(),
       nthreads = Threads.nthreads(), hostname = gethostname())
serialize(joinpath(DATA, "raw", "single_rhs_edges$(Int(CORRECT_EDGES)).jls"), out)

let csv = joinpath(DATA, "single_rhs.csv"),
    row = (; out.hostname, out.nthreads, out.smoke, out.correct_edges,
           out.n_points, out.n_src, out.niter, out.residual,
           out.t_load, out.t_precompute, out.t_solve, out.t_eval,
           out.rss_baseline_gb, out.rss_peak_gb, out.u_total_ev)
    newfile = !isfile(csv)
    open(csv, "a") do io
        newfile && println(io, join(string.(keys(row)), ","))
        println(io, join([sprint(print, v) for v in values(row)], ","))
    end
end
println("saved -> $(joinpath(DATA, "raw", "single_rhs_edges$(Int(CORRECT_EDGES)).jls"))")
println("ORBITAL BENCH DONE")

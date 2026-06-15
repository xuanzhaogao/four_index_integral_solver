# 6.5k — FULL production pipeline, per-stage runtime, with the NEW code:
# PrecomputedVolumeField (cache_fft) for the volume-field stages (interface
# build, RHS, u_int), and the UNCHANGED surface stages (LHS assembly, GMRES,
# scattered-potential eval — these use FMM3D, not FINUFFT, and never hit the
# 96-thread pathology). Mirrors run_single_rhs.jl's config and stages exactly
# so it is a clean old-vs-new comparison. Real monolayer pz orbital,
# calibrated parameters. All warm; on a clean idle node.
#
#   JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 \
#     julia --project=/mnt/home/xgao1/codes/BoundaryIntegral.jl \
#       exp65_orbital_bench/scripts/bench_full_field.jl
# env: CACHE_FFT (default 1), CORRECT_EDGES (default 1)

using BoundaryIntegral
import BoundaryIntegral as BI
using Krylov, LinearAlgebra, Printf, Serialization

const GRAPHENE = normpath(joinpath(@__DIR__, "..", "..", "..", "graphene"))
include(joinpath(GRAPHENE, "monolayer", "src", "MonolayerScreenedSolve.jl"))
using .MonolayerScreenedSolve
include(joinpath(GRAPHENE, "monolayer", "src", "MonolayerOrbitalLoader.jl"))
using .MonolayerOrbitalLoader

const CACHE_FFT = get(ENV, "CACHE_FFT", "1") == "1"
const CORRECT_EDGES = get(ENV, "CORRECT_EDGES", "1") == "1"
const DATA = joinpath(@__DIR__, "..", "data"); mkpath(joinpath(DATA, "raw"))

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")
const LZ, Z_CENTER, EPS_IN, EPS_OUT, L = 3.35, 3.35 / 2, 3.5, 1.0, 90.0
const P = (source_tol = 1e-3, n_quad = 6, edge_level = 4, rhs_tol = 1e-3,
           lhs_tol = 1e-5, gmres_atol = 1e-5, gmres_rtol = 1e-5,
           max_order = 64, max_depth = 12)
const E2_4PIEPS0 = 14.3996

rss() = round(Sys.maxrss() / 2^30; digits = 2)
const STAGES = Pair{String,Float64}[]
function note(lbl, t)
    push!(STAGES, lbl => t)
    @printf("  %-40s %8.2f s   maxrss %6.2f GB\n", lbl, t, rss()); flush(stdout)
end

# full pipeline; record=false for warm-up
function pipeline(svs, src_density_vs, Lx, Ly, Lz, prm; record::Bool)
    n(l, t) = record ? note(l, t) : nothing
    l_ec = Lz / 2.0^prm.edge_level * 1.01

    t = @elapsed field = PrecomputedVolumeField(svs; tol = prm.rhs_tol * 0.1,
                                                cache_fft = CACHE_FFT)
    n("field construction (cache_fft=$CACHE_FFT)", t)
    t = @elapsed interface = BI.single_dielectric_box3d_rhs_adaptive(
        Lx, Ly, Lz, prm.n_quad, field, 1.0, l_ec, prm.rhs_tol, EPS_IN, EPS_OUT, Float64;
        max_depth = prm.max_depth)
    n("interface build (field overload)", t)
    t = @elapsed lhs = BI.Lhs_dielectric_box3d_fmm3d_corrected(
        interface, prm.lhs_tol, prm.lhs_tol, prm.max_order; correct_edges = CORRECT_EDGES)
    n("LHS assembly (near corrections)", t)
    t = @elapsed rhs = rhs_dielectric_box3d_field(interface, field, 1.0)
    n("RHS assembly (field)", t)
    local sigma, stats
    t = @elapsed begin
        sigma, stats = Krylov.gmres(lhs, rhs; atol = prm.gmres_atol, rtol = prm.gmres_rtol, verbose = 0)
    end
    n("GMRES solve", t)
    residual = norm(lhs * sigma - rhs) / max(norm(rhs), eps(Float64))

    tmat = Matrix{Float64}(src_density_vs.positions)
    tw = src_density_vs.weights .* src_density_vs.density
    t = @elapsed u_int = volume_field_potential(field, tmat)
    n("eval: volume potential (field)", t)
    local u_scatter
    t = @elapsed begin
        op = BI.laplace3d_pottrg_fmm3d_corrected_hcubature(interface, tmat, prm.lhs_tol, prm.lhs_tol, 5.0)
        u_scatter = op * sigma
    end
    n("eval: scattered potential (FMM+hcub)", t)

    return (; n_points = BI.num_points(interface), n_src = length(svs.density),
            niter = stats.niter, residual,
            u_int_raw = dot(u_int, tw), u_scatter_raw = dot(u_scatter, tw))
end

@printf("threads = %d   cache_fft = %s   correct_edges = %s\n", Threads.nthreads(), CACHE_FFT, CORRECT_EDGES)

println(">>> warm-up (tiny Gaussian through every stage)"); flush(stdout)
let g = BI.GaussianVolumeSource((0.0, 0.0, 0.5), 0.3, 8, 1e-3)
    tw = @elapsed pipeline(g, g, 5.0, 5.0, 1.0, (; P..., edge_level = 2); record = false)
    @printf("  warm-up: %.1f s   baseline maxrss %.2f GB\n", tw, rss())
end

println(">>> load real pz orbital + screened source"); flush(stdout)
t_load = @elapsed begin
    dg = load_squared_xsf(XSF_1); sh = MonolayerOrbitalLoader._centering_shift(dg; tol = P.source_tol)
    dg = shift_datagrid(dg, (sh[1], sh[2], sh[3] + Z_CENTER))
    global vs1 = BI.VolumeSource(dg, tol = P.source_tol)
    global svs = BI.screened_volume_source(L, L, LZ, vs1, EPS_IN, EPS_OUT, BI.SharpScreening(); tol = P.rhs_tol)
end
note("load (XSF read + truncation + screening)", t_load)
Nphi1 = sum(vs1.weights .* vs1.density)
@printf("  Nphi1 = %.6f   source points = %d\n", Nphi1, length(vs1.density)); flush(stdout)

println(">>> full solve  L=$L  Lz=$LZ  eps_in=$EPS_IN  (warm)"); flush(stdout)
t_total = @elapsed res = pipeline(svs, vs1, L, L, LZ, P; record = true)

tof(l) = first(s.second for s in STAGES if s.first == l)
t_precompute = tof("field construction (cache_fft=$CACHE_FFT)") + tof("interface build (field overload)") + tof("LHS assembly (near corrections)")
t_solve = tof("RHS assembly (field)") + tof("GMRES solve")
t_eval = tof("eval: volume potential (field)") + tof("eval: scattered potential (FMM+hcub)")
u_onsite_ev = (res.u_int_raw + res.u_scatter_raw) * 4π * E2_4PIEPS0 / (Nphi1 * Nphi1)

println("=" ^ 72)
@printf("  interface points %d   src points %d   niter %d   residual %.3e\n",
        res.n_points, res.n_src, res.niter, res.residual)
@printf("  u_onsite = %.6f eV\n", u_onsite_ev)
println("-" ^ 72)
@printf("  %-12s %9.2f s\n", "load", t_load)
@printf("  %-12s %9.2f s\n", "precompute", t_precompute)
@printf("  %-12s %9.2f s\n", "solve", t_solve)
@printf("  %-12s %9.2f s\n", "eval", t_eval)
@printf("  %-12s %9.2f s   (end-to-end wall %.2f s)\n", "TOTAL", t_load + t_precompute + t_solve + t_eval, t_load + t_total)
@printf("  peak RSS %.2f GB\n", rss())
println("=" ^ 72)

serialize(joinpath(DATA, "raw", "bench_full_field.jls"),
          (; cache_fft = CACHE_FFT, correct_edges = CORRECT_EDGES, stages = copy(STAGES),
             t_load, t_precompute, t_solve, t_eval, n_points = res.n_points,
             niter = res.niter, u_onsite_ev, rss_gb = Sys.maxrss()/2^30, nthreads = Threads.nthreads()))
println("FULL FIELD BENCH DONE")

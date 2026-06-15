# 6.5c — production-scale benchmark of the MERGED-QUALITY PrecomputedVolumeField
# implementation (branch feature/precomputed-volume-field), on the real
# monolayer pz orbital. Run with the WORKTREE as the project so the new API is
# loaded without touching the live checkout or the shared numerical_results env:
#
#   JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 \
#     julia --project=/mnt/home/xgao1/codes/BoundaryIntegral.jl-wt-pvf \
#       exp65_orbital_bench/scripts/bench_field_path.jl
#
# Measures (96 threads, production parameters identical to run_single_rhs.jl):
#   1. field construction (type-1 + scaling + gradient coeffs, once)
#   2. adaptive mesh build via the field overload  [headline: was 201 s]
#   3. RHS assembly: rhs_dielectric_box3d_field vs Rhs_dielectric_box3d_fmm3d
#   4. volume potential at the 356k source points: volume_field_potential vs
#      ltkm3dc  [the u_int evaluation: was 44 s]
#   5. (RUN_VS=1, default) the production VolumeSource-path build in the same
#      session for an apples-to-apples A/B, plus mesh-equality check.
# BI deps only (no Printf): plain println formatting.

using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
using Serialization

const GRAPHENE = normpath(joinpath(@__DIR__, "..", "..", "..", "graphene"))
include(joinpath(GRAPHENE, "monolayer", "src", "MonolayerOrbitalLoader.jl"))
using .MonolayerOrbitalLoader

const RUN_VS = get(ENV, "RUN_VS", "1") == "1"
const DATA = joinpath(@__DIR__, "..", "data")
mkpath(joinpath(DATA, "raw"))

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")
const LZ, Z_CENTER, EPS_IN, EPS_OUT, L = 3.35, 3.35 / 2, 3.5, 1.0, 90.0
const N_QUAD, EDGE_LEVEL, RHS_TOL = 6, 4, 1e-3
const FMM_TOL = RHS_TOL * 0.1
const L_EC = LZ / 2.0^EDGE_LEVEL * 1.01

r2(x) = round(x; digits = 2)
rss() = round(Sys.maxrss() / 2^30; digits = 2)
line(lbl, t) = (println("  ", rpad(lbl, 46), lpad(string(r2(t)), 8), " s   maxrss ", rss(), " GB"); flush(stdout))

println("threads = ", Threads.nthreads(), "   run_vs = ", RUN_VS); flush(stdout)

# warm-up: tiny Gaussian through every code path, OLD and NEW (JIT).
# Both the production ltkm3dc u_int call and Rhs_dielectric_box3d_fmm3d are
# warmed here so the headline comparison is steady-state (the earlier 43.5 s
# u_int number was a COLD ltkm3dc first call — this removes that artifact).
println(">>> warm-up (old + new paths)"); flush(stdout)
let g = BI.GaussianVolumeSource((0.0, 0.0, 0.0), 0.3, 8, 1e-3)
    gpos = Matrix{Float64}(g.positions)
    gq = g.weights .* g.density
    gk = BI._estimate_tkm3dc_kmax(BI._estimate_source_spacing(g))
    f = PrecomputedVolumeField(g; tol = 1e-3)
    iface = BI.single_dielectric_box3d_rhs_adaptive(5.0, 5.0, 1.0, 2, f, 1.0, 0.6, 1e-2, EPS_IN, EPS_OUT, Float64; max_depth = 4)
    rhs_dielectric_box3d_field(iface, f, 1.0)
    volume_field_potential(f, gpos[:, 1:50])
    # OLD paths:
    BI.TKM3D.ltkm3dc(1e-2, gpos; charges = gq, targets = gpos[:, 1:50], pgt = 1, kmax = gk)
    BI.Rhs_dielectric_box3d_fmm3d(iface, g, 1.0, 1e-2)
    RUN_VS && BI.single_dielectric_box3d_rhs_adaptive(5.0, 5.0, 1.0, 2, g, 1.0, 0.6, 1e-2, EPS_IN, EPS_OUT, Float64; max_depth = 4)
end
println("  warm-up done   baseline maxrss ", rss(), " GB"); flush(stdout)

println(">>> load orbital + screened source"); flush(stdout)
t_load = @elapsed begin
    dg = load_squared_xsf(XSF_1)
    sh = MonolayerOrbitalLoader._centering_shift(dg; tol = 1e-3)
    dg = shift_datagrid(dg, (sh[1], sh[2], sh[3] + Z_CENTER))
    global vs1 = BI.VolumeSource(dg, tol = 1e-3)
    global svs = BI.screened_volume_source(L, L, LZ, vs1, EPS_IN, EPS_OUT, BI.SharpScreening(); tol = RHS_TOL)
end
const KMAX = BI._estimate_tkm3dc_kmax(BI._estimate_source_spacing(svs))
line("load (XSF + truncation + screening)", t_load)
println("  src points ", length(svs.density), "   kmax ", round(KMAX; digits = 2)); flush(stdout)

println(">>> NEW PATH (PrecomputedVolumeField)"); flush(stdout)
t_field = @elapsed field = PrecomputedVolumeField(svs; tol = FMM_TOL)
line("field construction (coeff + grad, once)", t_field)
println("  modes ", field.nmodes, " = ", prod(field.nmodes)); flush(stdout)

t_build = @elapsed iface = BI.single_dielectric_box3d_rhs_adaptive(
    L, L, LZ, N_QUAD, field, 1.0, L_EC, RHS_TOL, EPS_IN, EPS_OUT, Float64; max_depth = 12)
line("adaptive build (field overload)", t_build)
println("  interface points ", BI.num_points(iface), "   (production run: 960768)"); flush(stdout)

t_rhs_f = @elapsed rhs_f = rhs_dielectric_box3d_field(iface, field, 1.0)
line("RHS assembly (field)", t_rhs_f)

t_rhs_p = @elapsed rhs_p = BI.Rhs_dielectric_box3d_fmm3d(iface, svs, 1.0, RHS_TOL)
line("RHS assembly (production FMM)", t_rhs_p)
println("  RHS max rel dev: ", maximum(abs, rhs_f .- rhs_p) / maximum(abs, rhs_p)); flush(stdout)

tmat = Matrix{Float64}(vs1.positions)
tw = vs1.weights .* vs1.density
t_pot_f = @elapsed pot_f = volume_field_potential(field, tmat)
line("volume potential at 356k targets (field)", t_pot_f)

t_pot_p = @elapsed begin
    vals = BI.TKM3D.ltkm3dc(RHS_TOL, svs.positions; charges = svs.weights .* svs.density,
                            targets = tmat, pgt = 1, kmax = KMAX)
    vals.ier == 0 || error("ltkm3dc failed")
    global pot_p = real.(vals.pottarg)
end
line("volume potential (production ltkm3dc)", t_pot_p)
u_f = dot(tw, pot_f); u_p = dot(tw, pot_p)
println("  u_int (field) = ", u_f, "   u_int (ltkm3dc) = ", u_p,
        "   rel dev ", abs(u_f - u_p) / abs(u_p)); flush(stdout)

t_build_vs = NaN
if RUN_VS
    println(">>> PRODUCTION PATH (VolumeSource, same session A/B)"); flush(stdout)
    t_build_vs = @elapsed iface_vs = BI.single_dielectric_box3d_rhs_adaptive(
        L, L, LZ, N_QUAD, svs, 1.0, L_EC, RHS_TOL, EPS_IN, EPS_OUT, Float64;
        max_depth = 12, tkm_kmax = KMAX)
    line("adaptive build (VolumeSource path)", t_build_vs)
    println("  mesh equality: ", BI.num_points(iface_vs) == BI.num_points(iface),
            "  (", BI.num_points(iface_vs), " vs ", BI.num_points(iface), ")")
    println("  SPEEDUP build: ", round(t_build_vs / t_build; digits = 1), "x"); flush(stdout)
end

println("=" ^ 70)
println("  field construction      ", r2(t_field), " s")
println("  build: field ", r2(t_build), " s   vs production ", r2(t_build_vs),
        " s   (field+construction ", r2(t_field + t_build), " s)")
println("  RHS:   field ", r2(t_rhs_f), " s   vs production ", r2(t_rhs_p), " s")
println("  u_int: field ", r2(t_pot_f), " s   vs production ", r2(t_pot_p), " s")
println("  peak RSS ", rss(), " GB")
println("=" ^ 70)

serialize(joinpath(DATA, "raw", "bench_field_path.jls"),
          (; t_load, t_field, t_build, t_rhs_f, t_rhs_p, t_pot_f, t_pot_p, t_build_vs,
             n_points = BI.num_points(iface), nmodes = field.nmodes,
             rss_gb = Sys.maxrss() / 2^30, nthreads = Threads.nthreads()))
println("FIELD PATH BENCH DONE")

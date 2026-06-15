# 6.5e — production-scale validation of the cache_fft mode (commit 3b0a519):
# precomputed fine-grid FFT + native interp-only evaluation for in-box targets.
#
#   JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 \
#     julia --project=/mnt/home/xgao1/codes/BoundaryIntegral.jl \
#       exp65_orbital_bench/scripts/bench_cache_fft.jl

using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
using Serialization

const GRAPHENE = normpath(joinpath(@__DIR__, "..", "..", "..", "graphene"))
include(joinpath(GRAPHENE, "monolayer", "src", "MonolayerOrbitalLoader.jl"))
using .MonolayerOrbitalLoader

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
say(args...) = (println(args...); flush(stdout))

say("threads = ", Threads.nthreads())
say(">>> warm-up")
let g = BI.GaussianVolumeSource((0.0, 0.0, 0.0), 0.3, 8, 1e-3)
    for cf in (false, true)
        f = PrecomputedVolumeField(g; tol = 1e-3, cache_fft = cf)
        BI.single_dielectric_box3d_rhs_adaptive(5.0, 5.0, 1.0, 2, f, 1.0, 0.6, 1e-2, EPS_IN, EPS_OUT, Float64; max_depth = 3)
        volume_field_potential(f, Matrix{Float64}(g.positions[:, 1:50]))
    end
end

say(">>> load")
dg = load_squared_xsf(XSF_1)
sh = MonolayerOrbitalLoader._centering_shift(dg; tol = 1e-3)
dg = shift_datagrid(dg, (sh[1], sh[2], sh[3] + Z_CENTER))
vs1 = BI.VolumeSource(dg, tol = 1e-3)
svs = BI.screened_volume_source(L, L, LZ, vs1, EPS_IN, EPS_OUT, BI.SharpScreening(); tol = RHS_TOL)
tmat = Matrix{Float64}(vs1.positions)
tw = vs1.weights .* vs1.density

say(">>> standard field (cache_fft = false)")
t_con_s = @elapsed f_std = PrecomputedVolumeField(svs; tol = FMM_TOL)
t_pot_s = @elapsed pot_s = volume_field_potential(f_std, tmat)
u_s = dot(tw, pot_s)
say("  construction ", r2(t_con_s), " s   u_int eval ", r2(t_pot_s), " s   u = ", u_s, "   maxrss ", rss(), " GB")
f_std = nothing; GC.gc()

say(">>> cached field (cache_fft = true)")
t_con_c = @elapsed f_c = PrecomputedVolumeField(svs; tol = FMM_TOL, cache_fft = true)
t_pot_c = @elapsed pot_c = volume_field_potential(f_c, tmat)
u_c = dot(tw, pot_c)
say("  construction ", r2(t_con_c), " s   u_int eval ", r2(t_pot_c), " s   u = ", u_c, "   maxrss ", rss(), " GB")
say("  u_int agreement (cached vs standard): ", abs(u_c - u_s) / abs(u_s))
say("  pot max dev: ", maximum(abs, pot_c .- pot_s) / maximum(abs, pot_s))

say(">>> adaptive build with cached field")
t_build = @elapsed iface = BI.single_dielectric_box3d_rhs_adaptive(
    L, L, LZ, N_QUAD, f_c, 1.0, L_EC, RHS_TOL, EPS_IN, EPS_OUT, Float64; max_depth = 12)
say("  build ", r2(t_build), " s   points ", BI.num_points(iface), "   (must be 960768)")

t_rhs = @elapsed rhs_f = rhs_dielectric_box3d_field(iface, f_c, 1.0)
say("  RHS assembly ", r2(t_rhs), " s")
say("  peak RSS ", rss(), " GB")

serialize(joinpath(DATA, "raw", "bench_cache_fft.jls"),
          (; t_con_s, t_pot_s, u_s, t_con_c, t_pot_c, u_c, t_build, t_rhs,
             n_points = BI.num_points(iface), rss_gb = Sys.maxrss() / 2^30,
             nthreads = Threads.nthreads()))
say("CACHE FFT BENCH DONE")

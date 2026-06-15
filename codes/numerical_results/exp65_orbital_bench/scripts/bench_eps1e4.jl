# 6.5i — OLD vs NEW u_int with BOTH at eps = 1e-4, to separate the two
# sources of the field's speedup:
#   (reuse)  the field precomputes the type-1 spectrum once; the old path
#            rebuilds it every call.
#   (tol)    the production legacy path runs ltkm3dc at volume_tol = 1e-3,
#            where each NUFFT is much costlier than at 1e-4 (see
#            bench_tol_isolation.jl).
# By running the OLD ltkm3dc u_int at 1e-4 too (not its production 1e-3), the
# remaining gap vs the field is attributable to reuse alone.
#
# Reports, all warm: old ltkm3dc u_int at 1e-3 (production) and at 1e-4;
# field construction + volume_field_potential at 1e-4; cache_fft u_int at 1e-4.
#
#   JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 \
#     julia --project=/mnt/home/xgao1/codes/BoundaryIntegral.jl-wt-pvf \
#       exp65_orbital_bench/scripts/bench_eps1e4.jl

using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
const TKM = BI.TKM3D
const GRAPHENE = normpath(joinpath(@__DIR__, "..", "..", "..", "graphene"))
include(joinpath(GRAPHENE, "monolayer", "src", "MonolayerOrbitalLoader.jl"))
using .MonolayerOrbitalLoader

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")
const LZ, Z_CENTER, EPS_IN, EPS_OUT, L = 3.35, 3.35 / 2, 3.5, 1.0, 90.0
say(a...) = (println(a...); flush(stdout))
r2(x) = round(x; digits = 3)
tmin(f, n=3) = (f(); minimum(@elapsed(f()) for _ in 1:n))

function main()
    say("threads ", Threads.nthreads())
    dg = load_squared_xsf(XSF_1); sh = MonolayerOrbitalLoader._centering_shift(dg; tol = 1e-3)
    dg = shift_datagrid(dg, (sh[1], sh[2], sh[3] + Z_CENTER))
    vs1 = BI.VolumeSource(dg, tol = 1e-3)
    svs = BI.screened_volume_source(L, L, LZ, vs1, EPS_IN, EPS_OUT, BI.SharpScreening(); tol = 1e-3)
    src, q = BI._volume_source_fmm_sources(svs)
    trg = Matrix{Float64}(vs1.positions); tw = vs1.weights .* vs1.density
    kmax = BI._estimate_tkm3dc_kmax(BI._estimate_source_spacing(svs))
    say("src ", size(src,2), "  trg ", size(trg,2), "  kmax ", round(kmax;digits=2))

    # warm both paths small (avoids cold-FINUFFT crash, removes JIT from timings)
    let gs = BI.GaussianVolumeSource((0.0,0.0,0.0), 0.3, 8, 1e-3)
        gp = Matrix{Float64}(gs.positions); gq = gs.weights .* gs.density
        gk = BI._estimate_tkm3dc_kmax(BI._estimate_source_spacing(gs))
        TKM.ltkm3dc(1e-3, gp; charges=gq, targets=gp[:,1:50], pgt=1, kmax=gk)
        TKM.ltkm3dc(1e-4, gp; charges=gq, targets=gp[:,1:50], pgt=1, kmax=gk)
        f = PrecomputedVolumeField(gs; tol=1e-4); volume_field_potential(f, gp[:,1:50])
        fc = PrecomputedVolumeField(gs; tol=1e-4, cache_fft=true); volume_field_potential(fc, gp[:,1:50])
    end
    say(">>> warm-up done")

    # OLD ltkm3dc u_int at production tol 1e-3
    local u3
    t_old3 = tmin(() -> (global u3 = real.(TKM.ltkm3dc(1e-3, src; charges=q, targets=trg, pgt=1, kmax=kmax).pottarg)))
    say("  OLD ltkm3dc u_int @1e-3 (production):  ", r2(t_old3), " s   u=", dot(tw,u3))
    # OLD ltkm3dc u_int at 1e-4
    local u4
    t_old4 = tmin(() -> (global u4 = real.(TKM.ltkm3dc(1e-4, src; charges=q, targets=trg, pgt=1, kmax=kmax).pottarg)))
    say("  OLD ltkm3dc u_int @1e-4:               ", r2(t_old4), " s   u=", dot(tw,u4))

    # NEW field @1e-4: construction (once) + query
    local fld, uf
    t_con = @elapsed (fld = PrecomputedVolumeField(svs; tol=1e-4))
    t_q = tmin(() -> (global uf = volume_field_potential(fld, trg)))
    say("  NEW field @1e-4: construction ", r2(t_con), " s (once) + query ", r2(t_q), " s   u=", dot(tw,uf))

    # NEW cache_fft @1e-4
    local fldc, ufc
    t_conc = @elapsed (fldc = PrecomputedVolumeField(svs; tol=1e-4, cache_fft=true))
    t_qc = tmin(() -> (global ufc = volume_field_potential(fldc, trg)))
    say("  NEW cache_fft @1e-4: construction ", r2(t_conc), " s (once) + query ", r2(t_qc), " s   u=", dot(tw,ufc))

    say("\n--- attribution (u_int) ---")
    say("  tol effect (old 1e-3 -> 1e-4):   ", r2(t_old3), " -> ", r2(t_old4), " s  (", r2(t_old3/t_old4), "x)")
    say("  reuse effect (old 1e-4 -> field query 1e-4): ", r2(t_old4), " -> ", r2(t_q), " s  (", r2(t_old4/t_q), "x)")
    say("  cache_fft query 1e-4: ", r2(t_qc), " s  (", r2(t_old4/t_qc), "x vs old 1e-4)")
    say("EPS1E4 BENCH DONE")
end

main()

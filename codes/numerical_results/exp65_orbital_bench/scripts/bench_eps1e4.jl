# 6.5i — attribute the field's u_int speedup to its real causes, NOT tol
# (bench_tol_isolation.jl showed tol 1e-3 ~= 1e-4 at fixed grid). Measures,
# all on ONE node, all warm:
#   (A) OLD ltkm3dc u_int at 1e-3 (production) and at 1e-4  -> tol irrelevant?
#   (B) NEW field / cache_fft u_int at 1e-4                 -> reuse benefit
#   (C) a controlled TYPE-2 on the field-box grid vs the ltkm3dc per-call-box
#       grid, same node/tol/upsampfac                       -> grid-dimension effect
#
#   JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 \
#     julia --project=/mnt/home/xgao1/codes/BoundaryIntegral.jl \
#       exp65_orbital_bench/scripts/bench_eps1e4.jl

using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
const TKM = BI.TKM3D
const FN = TKM.FINUFFT
const GRAPHENE = normpath(joinpath(@__DIR__, "..", "..", "..", "graphene"))
include(joinpath(GRAPHENE, "monolayer", "src", "MonolayerOrbitalLoader.jl"))
using .MonolayerOrbitalLoader

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")
const LZ, Z_CENTER, EPS_IN, EPS_OUT, L = 3.35, 3.35 / 2, 3.5, 1.0, 90.0
say(a...) = (println(a...); flush(stdout))
r2(x) = round(x; digits = 3)

# min-of-n warm timing; returns (time, last_value) — no global/local dance
function tmin(f, n = 3)
    v = f()
    t = minimum(@elapsed(f()) for _ in 1:n)
    return t, v
end

# grid (axes/center/dks) for a given box-defining corner set + kmax
function grid_for(src, corners, kmax)
    lengths, center = TKM.combined_box_geometry_3xn(src, corners)
    Lbig = sqrt(sum(abs2, lengths)); dks = ntuple(d -> prevfloat(2π/(lengths[d]+Lbig)), 3)
    kx = TKM.centered_mode_axis(dks[1], kmax); ky = TKM.centered_mode_axis(dks[2], kmax); kz = TKM.centered_mode_axis(dks[3], kmax)
    return (; center, dks, nm = (length(kx),length(ky),length(kz)))
end

function time_type2_on(src, q, trg, G, tol)
    sx = G.dks[1].*(vec(src[1,:]).-G.center[1]); sy = G.dks[2].*(vec(src[2,:]).-G.center[2]); sz = G.dks[3].*(vec(src[3,:]).-G.center[3])
    c0 = FN.nufft3d1(sx, sy, sz, complex.(q), -1, tol, G.nm...); fk = ndims(c0)==4 ? dropdims(c0;dims=4) : c0
    tx = G.dks[1].*(vec(trg[1,:]).-G.center[1]); ty = G.dks[2].*(vec(trg[2,:]).-G.center[2]); tz = G.dks[3].*(vec(trg[3,:]).-G.center[3])
    t, _ = tmin() do
        plan = FN.finufft_makeplan(2, FN.BIGINT[G.nm...], 1, 1, tol; dtype=Float64)
        FN.finufft_setpts!(plan, tx, ty, tz); FN.finufft_exec(plan, fk); FN.finufft_destroy!(plan)
    end
    return t
end

function main()
    say("threads ", Threads.nthreads())
    dg = load_squared_xsf(XSF_1); sh = MonolayerOrbitalLoader._centering_shift(dg; tol = 1e-3)
    dg = shift_datagrid(dg, (sh[1], sh[2], sh[3] + Z_CENTER))
    vs1 = BI.VolumeSource(dg, tol = 1e-3)
    svs = BI.screened_volume_source(L, L, LZ, vs1, EPS_IN, EPS_OUT, BI.SharpScreening(); tol = 1e-3)
    src, q = BI._volume_source_fmm_sources(svs)
    trg = Matrix{Float64}(vs1.positions); tw = vs1.weights .* vs1.density
    h = BI._estimate_source_spacing(svs); kmax = BI._estimate_tkm3dc_kmax(h)
    say("src ", size(src,2), "  trg ", size(trg,2), "  kmax ", round(kmax;digits=2))

    # warm both paths small (avoids cold-FINUFFT crash, removes JIT)
    let gs = BI.GaussianVolumeSource((0.0,0.0,0.0), 0.3, 8, 1e-3)
        gp = Matrix{Float64}(gs.positions); gq = gs.weights .* gs.density
        gk = BI._estimate_tkm3dc_kmax(BI._estimate_source_spacing(gs))
        TKM.ltkm3dc(1e-3, gp; charges=gq, targets=gp[:,1:50], pgt=1, kmax=gk)
        TKM.ltkm3dc(1e-4, gp; charges=gq, targets=gp[:,1:50], pgt=1, kmax=gk)
        volume_field_potential(PrecomputedVolumeField(gs; tol=1e-4), gp[:,1:50])
        volume_field_potential(PrecomputedVolumeField(gs; tol=1e-4, cache_fft=true), gp[:,1:50])
    end
    say(">>> warm-up done")

    # (A) OLD ltkm3dc u_int, production 1e-3 vs 1e-4
    t_old3, u3 = tmin(() -> real.(TKM.ltkm3dc(1e-3, src; charges=q, targets=trg, pgt=1, kmax=kmax).pottarg))
    say("  (A) OLD ltkm3dc u_int @1e-3:  ", r2(t_old3), " s   u=", dot(tw,u3))
    t_old4, u4 = tmin(() -> real.(TKM.ltkm3dc(1e-4, src; charges=q, targets=trg, pgt=1, kmax=kmax).pottarg))
    say("  (A) OLD ltkm3dc u_int @1e-4:  ", r2(t_old4), " s   u=", dot(tw,u4))

    # (B) NEW field / cache_fft @1e-4
    t_con = @elapsed (fld = PrecomputedVolumeField(svs; tol=1e-4))
    t_q, uf = tmin(() -> volume_field_potential(fld, trg))
    say("  (B) NEW field @1e-4: construction ", r2(t_con), " s (once) + query ", r2(t_q), " s   u=", dot(tw,uf))
    t_conc = @elapsed (fldc = PrecomputedVolumeField(svs; tol=1e-4, cache_fft=true))
    t_qc, ufc = tmin(() -> volume_field_potential(fldc, trg))
    say("  (B) NEW cache_fft @1e-4: construction ", r2(t_conc), " s (once) + query ", r2(t_qc), " s   u=", dot(tw,ufc))

    # (C) controlled type-2: field box vs ltkm3dc per-call box, SAME node/tol
    lo = ntuple(d->minimum(view(src,d,:))-5h,3); hi = ntuple(d->maximum(view(src,d,:))+5h,3)
    Gfield = grid_for(src, [lo[1] hi[1]; lo[2] hi[2]; lo[3] hi[3]], kmax)        # source bbox + 5h
    Gcall  = grid_for(src, hcat(vec(minimum(hcat(src,trg);dims=2)), vec(maximum(hcat(src,trg);dims=2))), kmax)  # source∪target, no margin
    say("  (C) field-box grid ", Gfield.nm, "  vs  ltkm3dc-box grid ", Gcall.nm)
    say("      type-2 @1e-4 on field box:   ", r2(time_type2_on(src,q,trg,Gfield,1e-4)), " s")
    say("      type-2 @1e-4 on ltkm3dc box: ", r2(time_type2_on(src,q,trg,Gcall,1e-4)), " s")

    say("\n--- attribution (u_int) ---")
    say("  tol     (old 1e-3 -> 1e-4):              ", r2(t_old3), " -> ", r2(t_old4), " s   (", r2(t_old3/t_old4), "x)")
    say("  reuse   (old 1e-4 -> field query):       ", r2(t_old4), " -> ", r2(t_q),    " s   (", r2(t_old4/t_q), "x)")
    say("  cache   (field query -> cache_fft query):", r2(t_q),    " -> ", r2(t_qc),   " s   (", r2(t_q/t_qc), "x)")
    say("EPS1E4 BENCH DONE")
end

main()

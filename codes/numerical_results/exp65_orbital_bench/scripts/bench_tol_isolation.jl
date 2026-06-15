# 6.5h — ISOLATE the tolerance effect on NUFFT cost, controlled and
# reproducible. Holds the grid FIXED (one box, one kmax, one mode count) and
# varies ONLY the FINUFFT tolerance (1e-3 vs 1e-4) and upsampfac
# (auto / 1.25 / 2.0), for both the type-1 and type-2 transforms. Everything
# is function-wrapped (no top-level-global type instability — the earlier
# "51 s scale loop" was exactly that artifact), warmed, and reported as the
# min of NREP timed runs.
#
# Source: the real screened monolayer density by default; set REPRO_SYNTH=1
# for a self-contained synthetic Gaussian of similar size (portable repro,
# e.g. for a FINUFFT bug report).
#
#   JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 \
#     julia --project=/mnt/home/xgao1/codes/BoundaryIntegral.jl-wt-pvf \
#       exp65_orbital_bench/scripts/bench_tol_isolation.jl

using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
const TKM = BI.TKM3D
const FN = TKM.FINUFFT
const GRAPHENE = normpath(joinpath(@__DIR__, "..", "..", "..", "graphene"))
include(joinpath(GRAPHENE, "monolayer", "src", "MonolayerOrbitalLoader.jl"))
using .MonolayerOrbitalLoader

const SYNTH = get(ENV, "REPRO_SYNTH", "0") == "1"
const NREP = parse(Int, get(ENV, "NREP", "3"))
const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")
const LZ, Z_CENTER, EPS_IN, EPS_OUT, L = 3.35, 3.35 / 2, 3.5, 1.0, 90.0
say(a...) = (println(a...); flush(stdout))
r2(x) = round(x; digits = 2)

# min of NREP warm runs (one warm-up call first)
function tmin(f)
    f()
    minimum(@elapsed(f()) for _ in 1:NREP)
end

function build_source()
    if SYNTH
        n = 71                          # 71^3 = 357911 ~ production 355862
        g = range(-6.5, 7.0; length = n)
        pts = Matrix{Float64}(undef, 3, n^3); k = 0
        for x in g, y in g, z in range(-1.7, 5.0; length = n)
            k += 1; pts[1,k]=x; pts[2,k]=y; pts[3,k]=z
        end
        q = exp.(-(vec(sum(abs2, pts; dims=1))) ./ 8)   # smooth-ish blob
        return pts, q, copy(pts)
    else
        dg = load_squared_xsf(XSF_1); sh = MonolayerOrbitalLoader._centering_shift(dg; tol = 1e-3)
        dg = shift_datagrid(dg, (sh[1], sh[2], sh[3] + Z_CENTER))
        vs1 = BI.VolumeSource(dg, tol = 1e-3)
        svs = BI.screened_volume_source(L, L, LZ, vs1, EPS_IN, EPS_OUT, BI.SharpScreening(); tol = 1e-3)
        src, q = BI._volume_source_fmm_sources(svs)
        return src, q, Matrix{Float64}(vs1.positions)
    end
end

# ONE fixed grid (field-box convention: source bbox + 5h), reused for all tols
function fixed_grid(src)
    h = 0.0781                          # production mean spacing (kmax = pi/h ~ 40.2)
    kmax = BI._estimate_tkm3dc_kmax(h)
    lo = ntuple(d -> minimum(view(src,d,:)) - 5h, 3); hi = ntuple(d -> maximum(view(src,d,:)) + 5h, 3)
    lengths, center = TKM.combined_box_geometry_3xn(src, [lo[1] hi[1]; lo[2] hi[2]; lo[3] hi[3]])
    Lbig = sqrt(sum(abs2, lengths)); dks = ntuple(d -> prevfloat(2π/(lengths[d]+Lbig)), 3)
    kx = TKM.centered_mode_axis(dks[1], kmax); ky = TKM.centered_mode_axis(dks[2], kmax); kz = TKM.centered_mode_axis(dks[3], kmax)
    return (; kmax, center, dks, kx, ky, kz, nm=(length(kx),length(ky),length(kz)), Lbig)
end

usfkw(usf) = usf > 0 ? (; upsampfac = usf) : NamedTuple()

function time_type1(sxn, syn, szn, cq, nm, tol, usf)
    tmin(() -> FN.nufft3d1(sxn, syn, szn, cq, -1, tol, nm...; usfkw(usf)...))
end

function time_type2(txn, tyn, tzn, fk, nm, tol, usf)
    tmin() do
        plan = FN.finufft_makeplan(2, FN.BIGINT[nm...], 1, 1, tol; dtype=Float64, usfkw(usf)...)
        FN.finufft_setpts!(plan, txn, tyn, tzn)
        FN.finufft_exec(plan, fk)
        FN.finufft_destroy!(plan)
    end
end

function main()
    say("threads ", Threads.nthreads(), "   synthetic_source=", SYNTH, "   NREP=", NREP)
    src, q, trg = build_source()
    say("src ", size(src,2), "  trg ", size(trg,2))
    G = fixed_grid(src)
    say("FIXED grid: kmax ", r2(G.kmax), "  modes ", G.nm, " = ", prod(G.nm))
    cq = complex.(q)
    sxn = G.dks[1].*(vec(src[1,:]).-G.center[1]); syn = G.dks[2].*(vec(src[2,:]).-G.center[2]); szn = G.dks[3].*(vec(src[3,:]).-G.center[3])
    txn = G.dks[1].*(vec(trg[1,:]).-G.center[1]); tyn = G.dks[2].*(vec(trg[2,:]).-G.center[2]); tzn = G.dks[3].*(vec(trg[3,:]).-G.center[3])

    say("\n=== TYPE-1 (nufft3d1) on the SAME grid, tol x upsampfac ===")
    say("  tol      auto      usf=1.25   usf=2.0")
    for tol in (1e-3, 1e-4)
        ta = time_type1(sxn,syn,szn,cq,G.nm,tol,0.0); GC.gc()
        t125 = time_type1(sxn,syn,szn,cq,G.nm,tol,1.25); GC.gc()
        t20 = time_type1(sxn,syn,szn,cq,G.nm,tol,2.0); GC.gc()
        say("  ", tol, "   ", r2(ta), " s   ", r2(t125), " s    ", r2(t20), " s   (auto~=", t125<t20 ? "1.25" : "2.0", ")")
    end

    # build a representative coeff for the type-2 test (tol 1e-4, then reuse)
    c0 = FN.nufft3d1(sxn,syn,szn,cq,-1,1e-4,G.nm...); fk = ndims(c0)==4 ? dropdims(c0;dims=4) : c0
    say("\n=== TYPE-2 (makeplan+exec) on the SAME grid + coeff, tol x upsampfac ===")
    say("  tol      auto      usf=1.25   usf=2.0")
    for tol in (1e-3, 1e-4)
        ta = time_type2(txn,tyn,tzn,fk,G.nm,tol,0.0); GC.gc()
        t125 = time_type2(txn,tyn,tzn,fk,G.nm,tol,1.25); GC.gc()
        t20 = time_type2(txn,tyn,tzn,fk,G.nm,tol,2.0); GC.gc()
        say("  ", tol, "   ", r2(ta), " s   ", r2(t125), " s    ", r2(t20), " s   (auto~=", t125<t20 ? "1.25" : "2.0", ")")
    end
    say("\nINTERPRETATION: at FIXED grid + FIXED upsampfac, tighter tol should be")
    say("SLOWER (wider spread kernel). If 1e-4 is not slower than 1e-3, the earlier")
    say("field-vs-ltkm3dc gap was a grid/upsampfac difference, not tol per se.")
    say("TOL ISOLATION DONE")
end

main()

# 6.5j — dig INTO ltkm3dc: instrument a clone of TKM3D._ltkm3dc_eval that
# times each internal step (type-1 / scale / type-2) WITHIN one call, on the
# actual ltkm3dc per-call box grid AND the field box grid, on ONE node, warm.
# Also sweep the FINUFFT thread count (96 / 16 / 1) — leading suspect for the
# unexplained ~42 s full call is a threaded-FFT pathology at these dimensions.
# Cross-check: the instrumented sum must equal the real ltkm3dc(pgt=1) call.
#
#   JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 \
#     julia --project=/mnt/home/xgao1/codes/BoundaryIntegral.jl \
#       exp65_orbital_bench/scripts/diag_ltkm_internal.jl

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
const EPS = 1e-3
say(a...) = (println(a...); flush(stdout))
r2(x) = round(x; digits = 2)

# typed scaling loop (same as inside _ltkm3dc_eval)
function scale!(coeff, kx, ky, kz, kmax, Lbig)
    @inbounds for iz in eachindex(kz), iy in eachindex(ky), ix in eachindex(kx)
        k = sqrt(kx[ix]^2 + ky[iy]^2 + kz[iz]^2)
        coeff[ix, iy, iz] = k <= kmax ? coeff[ix, iy, iz] * TKM.truncated_laplace3d_hat(k, Lbig) :
                                        zero(eltype(coeff))
    end
end

# build the grid for a box defined by `corners` (3x2), like _ltkm3dc_eval does
function grid_for(src, corners, kmax)
    lengths, center = TKM.combined_box_geometry_3xn(src, corners)
    Lbig = sqrt(sum(abs2, lengths)); dks = ntuple(d -> prevfloat(2π/(lengths[d]+Lbig)), 3)
    kx = TKM.centered_mode_axis(dks[1], kmax); ky = TKM.centered_mode_axis(dks[2], kmax); kz = TKM.centered_mode_axis(dks[3], kmax)
    return (; center, dks, kx, ky, kz, nm=(length(kx),length(ky),length(kz)), Lbig)
end

# instrumented _ltkm3dc_eval (pgt=1): returns (t_type1, t_scale, t_type2)
function eval_timed(src, q, trg, G, kmax, eps, nthreads)
    nthkw = nthreads > 0 ? (; nthreads = nthreads) : NamedTuple()
    sx = G.dks[1].*(vec(src[1,:]).-G.center[1]); sy = G.dks[2].*(vec(src[2,:]).-G.center[2]); sz = G.dks[3].*(vec(src[3,:]).-G.center[3])
    tx = G.dks[1].*(vec(trg[1,:]).-G.center[1]); ty = G.dks[2].*(vec(trg[2,:]).-G.center[2]); tz = G.dks[3].*(vec(trg[3,:]).-G.center[3])
    cq = complex.(q)
    t1 = @elapsed c0 = FN.nufft3d1(sx, sy, sz, cq, -1, eps, G.nm...; nthkw...)
    coeff = ndims(c0)==4 ? dropdims(c0;dims=4) : c0
    t2 = @elapsed scale!(coeff, G.kx, G.ky, G.kz, kmax, G.Lbig)
    t3 = @elapsed begin
        plan = FN.finufft_makeplan(2, FN.BIGINT[G.nm...], 1, 1, eps; dtype=Float64, nthkw...)
        FN.finufft_setpts!(plan, tx, ty, tz); FN.finufft_exec(plan, coeff); FN.finufft_destroy!(plan)
    end
    return t1, t2, t3
end

function main()
    say("threads(julia) ", Threads.nthreads(), "   FINUFFT nthreads probed: 96 / 16 / 1   eps=", EPS)
    dg = load_squared_xsf(XSF_1); sh = MonolayerOrbitalLoader._centering_shift(dg; tol = 1e-3)
    dg = shift_datagrid(dg, (sh[1], sh[2], sh[3] + Z_CENTER))
    vs1 = BI.VolumeSource(dg, tol = 1e-3)
    svs = BI.screened_volume_source(L, L, LZ, vs1, EPS_IN, EPS_OUT, BI.SharpScreening(); tol = 1e-3)
    src, q = BI._volume_source_fmm_sources(svs)
    trg = Matrix{Float64}(vs1.positions)
    h = BI._estimate_source_spacing(svs); kmax = BI._estimate_tkm3dc_kmax(h)
    say("src ", size(src,2), "  trg ", size(trg,2), "  kmax ", r2(kmax))

    sb = hcat(vec(minimum(src;dims=2)), vec(maximum(src;dims=2)))
    fieldcorners = hcat(sb[:,1] .- 5h, sb[:,2] .+ 5h)                       # field box (source bbox + 5h)
    callcorners  = hcat(vec(minimum(hcat(src,trg);dims=2)), vec(maximum(hcat(src,trg);dims=2)))  # ltkm3dc box (source∪target)
    Gfield = grid_for(src, fieldcorners, kmax)
    Gcall  = grid_for(src, callcorners,  kmax)
    say("field-box grid   ", Gfield.nm, " = ", prod(Gfield.nm))
    say("ltkm3dc-box grid ", Gcall.nm,  " = ", prod(Gcall.nm))

    # reference: the REAL ltkm3dc full call (warm), default threads
    TKM.ltkm3dc(EPS, src[:,1:2000]; charges=q[1:2000], targets=trg[:,1:50], pgt=1, kmax=kmax)  # warm
    t_real = @elapsed TKM.ltkm3dc(EPS, src; charges=q, targets=trg, pgt=1, kmax=kmax)
    say(">>> REAL ltkm3dc(pgt=1) full call (default threads): ", r2(t_real), " s")

    say("\n  box        nthreads   type-1     scale     type-2     sum")
    for (name, G) in (("field  ", Gfield), ("ltkm3dc", Gcall))
        for nth in (96, 16, 1)
            eval_timed(src, q, trg, G, kmax, EPS, nth)              # warm this (nm, nthreads)
            t1, t2, t3 = eval_timed(src, q, trg, G, kmax, EPS, nth) # timed
            say("  ", name, "    ", lpad(nth,3), "      ", lpad(r2(t1),6), " s  ", lpad(r2(t2),6), " s  ", lpad(r2(t3),6), " s  ", lpad(r2(t1+t2+t3),6), " s")
            GC.gc()
        end
    end
    say("\nIf ltkm3dc-box sum @96 ~= REAL full call, the decomposition is faithful.")
    say("If type-1/type-2 drop sharply at fewer threads, it is a threaded-FFT pathology.")
    say("If field-box << ltkm3dc-box at the same nthreads, it is a grid-dimension effect.")
    say("LTKM INTERNAL DIAG DONE")
end

main()

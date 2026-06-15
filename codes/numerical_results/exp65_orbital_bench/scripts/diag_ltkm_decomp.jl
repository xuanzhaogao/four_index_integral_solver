# 6.5g — DECOMPOSE the warm ltkm3dc u_int call (~56 s) to find where the time
# goes vs the field's 0.41 s type-2. Runs the field construction FIRST so
# FINUFFT/FFTW thread state is initialized (the cold first big transform
# silently killed julia on two nodes). All timings WARM.
#
#   JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 \
#     julia --project=/mnt/home/xgao1/codes/BoundaryIntegral.jl \
#       exp65_orbital_bench/scripts/diag_ltkm_decomp.jl

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
const RHS_TOL = 1e-3          # tol the OLD ltkm3dc u_int call uses (production)
const FMM_TOL = 1e-4          # tol the field uses
say(a...) = (println(a...); flush(stdout))
r2(x) = round(x; digits = 2)

say(">>> load")
dg = load_squared_xsf(XSF_1); sh = MonolayerOrbitalLoader._centering_shift(dg; tol = 1e-3)
dg = shift_datagrid(dg, (sh[1], sh[2], sh[3] + Z_CENTER))
vs1 = BI.VolumeSource(dg, tol = 1e-3)
svs = BI.screened_volume_source(L, L, LZ, vs1, EPS_IN, EPS_OUT, BI.SharpScreening(); tol = RHS_TOL)
src, q = BI._volume_source_fmm_sources(svs)
trg = Matrix{Float64}(vs1.positions)
h = BI._estimate_source_spacing(svs); kmax = BI._estimate_tkm3dc_kmax(h)
say("  src ", size(src,2), "  trg ", size(trg,2), "  kmax ", r2(kmax))

# init FINUFFT/FFTW threads via the field construction (avoids cold-call crash)
say(">>> field construction (also warms FINUFFT)")
t_field = @elapsed fld = PrecomputedVolumeField(svs; tol = FMM_TOL)
say("  construction ", r2(t_field), " s   field grid ", fld.nmodes)
volume_field_potential(fld, trg[:, 1:100])
t_field_u = @elapsed volume_field_potential(fld, trg)
say("  field volume_field_potential (warm, type-2 only): ", round(t_field_u; digits=3), " s")

# ltkm3dc's own grid (combined source+target box, kmax passed)
lengths, center = TKM.combined_box_geometry_3xn(src, trg)
Lbig = sqrt(sum(abs2, lengths))
dks = ntuple(d -> prevfloat(2π/(lengths[d] + Lbig)), 3)
ax(km) = ntuple(d -> TKM.centered_mode_axis(dks[d], km), 3)
kx,ky,kz = ax(kmax); nm = (length(kx),length(ky),length(kz))
say(">>> ltkm3dc grid (kmax=", r2(kmax), ", tol=", RHS_TOL, "): modes ", nm, " = ", prod(nm))

# warm small ltkm3dc then time the full warm call
TKM.ltkm3dc(RHS_TOL, src[:,1:2000]; charges=q[1:2000], targets=trg[:,1:50], pgt=1, kmax=kmax)
t_full = @elapsed TKM.ltkm3dc(RHS_TOL, src; charges=q, targets=trg, pgt=1, kmax=kmax)
say(">>> FULL ltkm3dc(pgt=1) warm: ", r2(t_full), " s")

# --- decompose on ltkm3dc's grid, tol = RHS_TOL ---
sxn = dks[1].*(vec(src[1,:]).-center[1]); syn = dks[2].*(vec(src[2,:]).-center[2]); szn = dks[3].*(vec(src[3,:]).-center[3])
cq = complex.(q)
# type-1 at auto upsampfac (what ltkm3dc uses)
FN.nufft3d1(sxn[1:2000],syn[1:2000],szn[1:2000],cq[1:2000],-1,RHS_TOL,32,32,32)
t_t1_auto = @elapsed c_auto = FN.nufft3d1(sxn,syn,szn,cq,-1,RHS_TOL,nm...)
say("  type-1 nufft3d1 AUTO upsampfac:  ", r2(t_t1_auto), " s")
c_auto = nothing; GC.gc()
for usf in (1.25, 2.0)
    FN.nufft3d1(sxn[1:2000],syn[1:2000],szn[1:2000],cq[1:2000],-1,RHS_TOL,32,32,32; upsampfac=usf)
    t = @elapsed c = FN.nufft3d1(sxn,syn,szn,cq,-1,RHS_TOL,nm...; upsampfac=usf)
    say("  type-1 nufft3d1 upsampfac=", usf, ":   ", r2(t), " s")
    c = nothing; GC.gc()
end
# scale loop cost on this grid
c = FN.nufft3d1(sxn,syn,szn,cq,-1,RHS_TOL,nm...)
c = ndims(c)==4 ? dropdims(c;dims=4) : c
t_scale = @elapsed begin
    @inbounds for iz in eachindex(kz), iy in eachindex(ky), ix in eachindex(kx)
        k = sqrt(kx[ix]^2+ky[iy]^2+kz[iz]^2)
        c[ix,iy,iz] = k<=kmax ? c[ix,iy,iz]*TKM.truncated_laplace3d_hat(k,Lbig) : zero(eltype(c))
    end
end
say("  scale loop (per-mode kernel):    ", r2(t_scale), " s")
# type-2 at auto
txn=dks[1].*(vec(trg[1,:]).-center[1]); tyn=dks[2].*(vec(trg[2,:]).-center[2]); tzn=dks[3].*(vec(trg[3,:]).-center[3])
t_t2 = @elapsed begin
    plan = TKM._finufft_make_type2_plan_3d(txn,tyn,tzn,1,RHS_TOL,nm,1,Float64)
    FN.finufft_exec(plan, c); FN.finufft_destroy!(plan)
end
say("  type-2 makeplan+exec AUTO:       ", r2(t_t2), " s")
say(">>> sum of decomposed (auto): ", r2(t_t1_auto + t_scale + t_t2), " s   vs full ", r2(t_full), " s")
say("DECOMP DONE")

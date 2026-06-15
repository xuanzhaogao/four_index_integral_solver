# 6.5f — DIAGNOSTIC: why was each OLD ltkm3dc call so expensive? Print the
# actual Fourier grid the old per-call path builds vs the field's fixed grid,
# and decompose the old call cold-vs-warm and type-1-vs-type-2.
#
#   JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 \
#     julia --project=/mnt/home/xgao1/codes/BoundaryIntegral.jl \
#       exp65_orbital_bench/scripts/diag_tkm_grid.jl

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
const RHS_TOL = 1e-3
const FMM_TOL = RHS_TOL * 0.1
say(a...) = (println(a...); flush(stdout))

# the Fourier grid that _ltkm3dc_eval / the field constructor builds for a box
function grid_of(src_box, kmax)   # src_box = 3x2 [min max] used to size the box
    lengths, center = TKM.combined_box_geometry_3xn(src_box, src_box)
    Lbig = sqrt(sum(abs2, lengths))
    dks = ntuple(d -> prevfloat(2π / (lengths[d] + Lbig)), 3)
    nk = ntuple(d -> length(TKM.centered_mode_axis(dks[d], kmax)), 3)
    return lengths, Lbig, dks, nk
end
bbox(M) = hcat(vec(minimum(M; dims = 2)), vec(maximum(M; dims = 2)))

say(">>> load")
dg = load_squared_xsf(XSF_1)
sh = MonolayerOrbitalLoader._centering_shift(dg; tol = 1e-3)
dg = shift_datagrid(dg, (sh[1], sh[2], sh[3] + Z_CENTER))
vs1 = BI.VolumeSource(dg, tol = 1e-3)                       # target grid (raw)
svs = BI.screened_volume_source(L, L, LZ, vs1, EPS_IN, EPS_OUT, BI.SharpScreening(); tol = RHS_TOL)
src, q = BI._volume_source_fmm_sources(svs)
trg = Matrix{Float64}(vs1.positions)
h = BI._estimate_source_spacing(svs)
kmax = BI._estimate_tkm3dc_kmax(h)
say("  src ", size(src, 2), "  targets ", size(trg, 2), "  h ", round(h; digits = 4), "  kmax ", round(kmax; digits = 2))

# ---- the two boxes ----
sb = bbox(src); tb = bbox(trg)
say("  source bbox:  x [", round(sb[1,1];digits=2), ", ", round(sb[1,2];digits=2), "]  y [",
    round(sb[2,1];digits=2), ", ", round(sb[2,2];digits=2), "]  z [", round(sb[3,1];digits=2), ", ", round(sb[3,2];digits=2), "]")
say("  target bbox:  x [", round(tb[1,1];digits=2), ", ", round(tb[1,2];digits=2), "]  y [",
    round(tb[2,1];digits=2), ", ", round(tb[2,2];digits=2), "]  z [", round(tb[3,1];digits=2), ", ", round(tb[3,2];digits=2), "]")

# FIELD grid: source bbox + 5h margin
fb = hcat(sb[:,1] .- 5h, sb[:,2] .+ 5h)
lf, Lf, dkf, nkf = grid_of(fb, kmax)
say(">>> FIELD box (source bbox + 5h):  lengths ", round.(lf; digits=2), "  diag ", round(Lf;digits=2),
    "  -> modes ", nkf, " = ", prod(nkf))

# OLD ltkm3dc grid: combined(sources, targets), kmax passed = 40.2
ob = hcat(min.(sb[:,1], tb[:,1]), max.(sb[:,2], tb[:,2]))
lo, Lo, dko, nko = grid_of(ob, kmax)
say(">>> OLD per-call box combined(src,trg), kmax=", round(kmax;digits=2), ":  lengths ", round.(lo;digits=2),
    "  diag ", round(Lo;digits=2), "  -> modes ", nko, " = ", prod(nko))

# (contrast, from earlier bench_meshgen run: estimate_kcut3dc(tol=1e-4) -> kcut=2201,
#  which on the combined box would be ~24000^3 ~ 1e13 modes; production avoids this
#  by passing kmax=40.2 explicitly. Not recomputed here — it OOM'd this script.)

# ---- timing: cold vs warm, decomposed ----
say(">>> timing the OLD u_int call: ltkm3dc(src, targets=trg, kmax=", round(kmax;digits=2), ", pgt=1)")
t_cold = @elapsed v1 = TKM.ltkm3dc(RHS_TOL, src; charges = q, targets = trg, pgt = 1, kmax = kmax)
say("   COLD (first call, incl. JIT compile): ", round(t_cold; digits=2), " s")
t_warm = @elapsed v2 = TKM.ltkm3dc(RHS_TOL, src; charges = q, targets = trg, pgt = 1, kmax = kmax)
say("   WARM (second call):                   ", round(t_warm; digits=2), " s")

# decompose the warm call: type-1 + scale + type-2
say(">>> decomposed (warm) on the old grid ", nko)
sxn = dko[1] .* (vec(src[1,:]) .- (ob[1,1]+ob[1,2])/2)
syn = dko[2] .* (vec(src[2,:]) .- (ob[2,1]+ob[2,2])/2)
szn = dko[3] .* (vec(src[3,:]) .- (ob[3,1]+ob[3,2])/2)
FN.nufft3d1(sxn[1:5000], syn[1:5000], szn[1:5000], complex.(q[1:5000]), -1, RHS_TOL, 32,32,32)  # warm
t_t1 = @elapsed c = FN.nufft3d1(sxn, syn, szn, complex.(q), -1, RHS_TOL, nko...)
say("   type-1 (warm): ", round(t_t1; digits=2), " s")

# field path for reference (warm)
fld = PrecomputedVolumeField(svs; tol = FMM_TOL)
volume_field_potential(fld, trg[:, 1:100])  # warm
t_field_q = @elapsed volume_field_potential(fld, trg)
say(">>> field volume_field_potential (warm, type-2 only): ", round(t_field_q; digits=3), " s")
say("DIAG DONE")

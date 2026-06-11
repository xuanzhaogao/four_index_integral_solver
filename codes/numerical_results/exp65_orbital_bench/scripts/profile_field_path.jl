# 6.5d — step-by-step profiling of the PrecomputedVolumeField pipeline
# (post FMM-far policy, commit cc6a749) + optimization experiments.
#
# Run with the worktree as project:
#   JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 \
#     julia --project=/mnt/home/xgao1/codes/BoundaryIntegral.jl-wt-pvf \
#       exp65_orbital_bench/scripts/profile_field_path.jl
#
# Parts:
#   A  field construction decomposed: type-1 NUFFT / kernel scaling / gradient
#      coefficients — each at FINUFFT upsampfac 2.0 (default) and 1.25
#      (smaller internal FFT grid), with coefficient agreement check.
#   B  type-2 evaluation decomposed (makeplan / setpts / exec) at the three
#      production batch sizes: ~2k (per-depth in-box), 356k (u_int targets),
#      ntrans 1 and 3, upsampfac 2.0 vs 1.25, value agreement check.
#   C  far-path FMM floor at per-depth batch sizes (for the record).
#   D  end-to-end: build + RHS + u_int with the current branch code.
# BI deps only (no Printf).

using BoundaryIntegral
import BoundaryIntegral as BI
using LinearAlgebra
using Serialization

const TKM = BI.TKM3D
const FN = TKM.FINUFFT

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

r3(x) = round(x; digits = 3)
say(args...) = (println(args...); flush(stdout))

say("threads = ", Threads.nthreads())

say(">>> setup")
dg = load_squared_xsf(XSF_1)
sh = MonolayerOrbitalLoader._centering_shift(dg; tol = 1e-3)
dg = shift_datagrid(dg, (sh[1], sh[2], sh[3] + Z_CENTER))
vs1 = BI.VolumeSource(dg, tol = 1e-3)
svs = BI.screened_volume_source(L, L, LZ, vs1, EPS_IN, EPS_OUT, BI.SharpScreening(); tol = RHS_TOL)
src, q = BI._volume_source_fmm_sources(svs)
h = BI._estimate_source_spacing(svs)
km = BI._estimate_tkm3dc_kmax(h)
say("  src ", length(q), "  kmax ", round(km; digits = 2))

# Fourier box exactly as the constructor builds it
lo = ntuple(d -> minimum(view(src, d, :)) - 5.0 * h, 3)
hi = ntuple(d -> maximum(view(src, d, :)) + 5.0 * h, 3)
corners = [lo[1] hi[1]; lo[2] hi[2]; lo[3] hi[3]]
lengths, center = TKM.combined_box_geometry_3xn(src, corners)
Lbig = sqrt(sum(abs2, lengths))
dks = ntuple(d -> prevfloat(2π / (lengths[d] + Lbig)), 3)
kx = TKM.centered_mode_axis(dks[1], km); ky = TKM.centered_mode_axis(dks[2], km)
kz = TKM.centered_mode_axis(dks[3], km)
nm = (length(kx), length(ky), length(kz))
say("  modes ", nm, " = ", prod(nm))
sxn = dks[1] .* (vec(src[1, :]) .- center[1])
syn = dks[2] .* (vec(src[2, :]) .- center[2])
szn = dks[3] .* (vec(src[3, :]) .- center[3])
cq = complex.(q)

scale!(c) = (@inbounds for iz in eachindex(kz), iy in eachindex(ky), ix in eachindex(kx)
    k = sqrt(kx[ix]^2 + ky[iy]^2 + kz[iz]^2)
    c[ix, iy, iz] = k <= km ? c[ix, iy, iz] * TKM.truncated_laplace3d_hat(k, Lbig) :
                              zero(eltype(c))
end; c)

say(">>> PART A: construction stages, upsampfac 2.0 vs 1.25")
say("  (warm-up type-1 at each upsampfac)")
FN.nufft3d1(sxn[1:5000], syn[1:5000], szn[1:5000], cq[1:5000], -1, FMM_TOL, 64, 64, 64)
FN.nufft3d1(sxn[1:5000], syn[1:5000], szn[1:5000], cq[1:5000], -1, FMM_TOL, 64, 64, 64; upsampfac = 1.25)
coeffs = Dict{Float64, Array{ComplexF64, 3}}()
for usf in (2.0, 1.25)
    t1 = @elapsed c0 = FN.nufft3d1(sxn, syn, szn, cq, -1, FMM_TOL, nm...; upsampfac = usf)
    c = ndims(c0) == 4 ? dropdims(c0; dims = 4) : c0
    t2 = @elapsed scale!(c)
    t3 = @elapsed gc_ = TKM._spectral_gradient_coeffs_3d(c, kx, ky, kz)
    say("  usf=", usf, "  type-1 ", r3(t1), " s   scale ", r3(t2), " s   gradcoeff ", r3(t3), " s")
    coeffs[usf] = c
    gc_ = nothing; GC.gc()
end
let a = coeffs[2.0], b = coeffs[1.25]
    dev = maximum(abs, a .- b) / maximum(abs, a)
    say("  type-1 coeff agreement (rel, max): ", dev)
end
coeff = coeffs[2.0]
gradc = TKM._spectral_gradient_coeffs_3d(coeff, kx, ky, kz)
coeffs = nothing; GC.gc()

say(">>> PART B: type-2 decomposed (makeplan/setpts/exec), usf 2.0 vs 1.25")
prefactor = prod(dks) / (2π)^3
trg_in(n) = begin   # points inside the box, near the density
    t = src[:, 1:n] .+ 0.5 .* h .* (rand(3, n) .- 0.5)
    t
end
for nt in (2000, 356000), (lbl, fk, ntrans) in (("pot", coeff, 1), ("grad", gradc, 3))
    trg = trg_in(min(nt, size(src, 2)))
    txn = dks[1] .* (vec(trg[1, :]) .- center[1])
    tyn = dks[2] .* (vec(trg[2, :]) .- center[2])
    tzn = dks[3] .* (vec(trg[3, :]) .- center[3])
    vals_ref = nothing
    for usf in (2.0, 1.25)
        tplan = @elapsed plan = FN.finufft_makeplan(2, FN.BIGINT[nm...], 1, ntrans, FMM_TOL;
                                                    dtype = Float64, upsampfac = usf)
        tset = @elapsed FN.finufft_setpts!(plan, txn, tyn, tzn)
        texec = @elapsed v = FN.finufft_exec(plan, fk)
        FN.finufft_destroy!(plan)
        if usf == 2.0
            vals_ref = v
        else
            dev = maximum(abs, v .- vals_ref) / maximum(abs, vals_ref)
            say("    [usf 1.25 vs 2.0 value agreement: ", dev, "]")
        end
        say("  n=", size(trg, 2), " ", lbl, " (ntrans=", ntrans, ") usf=", usf,
            "  plan ", r3(tplan), " s   setpts ", r3(tset), " s   exec ", r3(texec), " s")
    end
end
gradc = nothing; GC.gc()

say(">>> PART C: far-path FMM floor (per-depth batch sizes)")
trg_far(n) = hcat([[lo[1] + rand() * (hi[1] - lo[1]),
                    lo[2] + rand() * (hi[2] - lo[2]),
                    hi[3] + 1.0 + rand() * 30.0] for _ in 1:n]...)
BI.lfmm3d(FMM_TOL, src; charges = q, targets = trg_far(64), pgt = 2)
for n in (1000, 4000, 10000)
    t = @elapsed BI.lfmm3d(FMM_TOL, src; charges = q, targets = trg_far(n), pgt = 2)
    say("  lfmm3d pgt=2  ", n, " targets: ", r3(t), " s")
end

say(">>> PART D: end-to-end with current branch code")
t_field = @elapsed field = PrecomputedVolumeField(svs; tol = FMM_TOL)
say("  field construction: ", r3(t_field), " s")
t_build = @elapsed iface = BI.single_dielectric_box3d_rhs_adaptive(
    L, L, LZ, N_QUAD, field, 1.0, L_EC, RHS_TOL, EPS_IN, EPS_OUT, Float64; max_depth = 12)
say("  adaptive build (FMM far policy): ", r3(t_build), " s   points ", BI.num_points(iface))
t_rhs = @elapsed rhs_f = rhs_dielectric_box3d_field(iface, field, 1.0)
say("  RHS assembly (field): ", r3(t_rhs), " s")
tmat = Matrix{Float64}(vs1.positions)
t_pot = @elapsed pot = volume_field_potential(field, tmat)
say("  u_int potential (field): ", r3(t_pot), " s   u = ", dot(vs1.weights .* vs1.density, pot))
say("  peak RSS ", round(Sys.maxrss() / 2^30; digits = 2), " GB")

serialize(joinpath(DATA, "raw", "profile_field_path.jls"),
          (; nm, t_field, t_build, t_rhs, t_pot, nthreads = Threads.nthreads(),
             rss_gb = Sys.maxrss() / 2^30))
say("PROFILE DONE")

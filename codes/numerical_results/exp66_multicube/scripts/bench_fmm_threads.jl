# Isolate the 96-thread evaluation collapse: is it FMM3D or FINUFFT?
#
# The eval splits 7.17M targets 88.9% far (lfmm3d, full OpenMP pool) / 11.1% near
# (PrecomputedVolumeField -> FINUFFT, capped at 16 threads by volume_field.jl). Since the NUFFT
# path is capped it is identical at 64 and 96 threads and cannot explain a 4.6x difference --
# this times the two paths separately at one thread count to confirm that directly.
#
# Run one thread count per process (OpenMP reads its pool at startup):
#   JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 julia --project=. exp66_multicube/scripts/bench_fmm_threads.jl
using BoundaryIntegral, Serialization, Printf
import BoundaryIntegral as BI
using FMM3D

const REF = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
TEMPL = [BI.read_xsf(joinpath(REF, "graphene_0000$(i).xsf"))[2] for i in 1:2]
CEN = [ntuple(d -> Float64(BI.density_centroid(dg)[d]), 3) for dg in TEMPL]
b = BI.assemble_lattice_batch(TEMPL, Dict(1 => BI.OrbitalInstance(1, 1,
        BI.snap_orbital(TEMPL[1], CEN[1], CEN[1]))), [(1,1)]; support_rtol = 1e-4)
src = BI.batch_volume_sources(b)[1]
L = 270.0; tz = 9.0; cx,cy,cz = CEN[1]; h = L/2; czi = cz - tz/2 - h
boxes = BI.BoxGeom[(center=(cx-h,cy,czi),Lx=L,Ly=L,Lz=L),(center=(cx+h,cy,czi),Lx=L,Ly=L,Lz=L),
                   (center=(cx,cy,cz),Lx=42.95357656,Ly=44.06001031,Lz=tz)]
scr = BI.screened_volume_source(boxes, Float64[11.9,3.9,2.4], 1.0, src, BI.SharpScreening())

tc = deserialize(joinpath(@__DIR__, "..", "data_sec53evalpilot", "targets_model.jls"))
T = tc.positions
geom = BI.near_field_geometry(scr; c_pad = 5.0)
near_idx = findall(i -> BI.in_near_region(geom, T, i), 1:size(T,2))
far_idx  = setdiff(1:size(T,2), near_idx)
@printf("threads: julia %d  OMP %s   near %d  far %d\n",
        Threads.nthreads(), get(ENV,"OMP_NUM_THREADS","unset"), length(near_idx), length(far_idx))

q = scr.weights .* scr.density
tol = 1e-3
# --- FAR: the uncapped FMM3D call, 88.9% of the targets
far_t = T[:, far_idx]
lfmm3d(tol, scr.positions; charges = reshape(q, 1, :), targets = far_t[:, 1:1000], pgt = 1, nd = 1)  # warm
t_fmm = @elapsed lfmm3d(tol, scr.positions; charges = reshape(q, 1, :), targets = far_t, pgt = 1, nd = 1)
@printf("  lfmm3d  (FMM3D, full pool)      %8.2f s   over %d targets\n", t_fmm, length(far_idx))

# --- NEAR: the FINUFFT path, capped at 16 threads inside BI
near_t = T[:, near_idx]
BI.PrecomputedVolumeField(scr; tol = tol, c_pad = 5.0, compute_grad = false)                       # warm
t_nufft = @elapsed begin
    fld = BI.PrecomputedVolumeField(scr; tol = tol, c_pad = 5.0, compute_grad = false)
    BI.volume_field_potential(fld, near_t)
end
@printf("  PrecomputedVolumeField (NUFFT)  %8.2f s   over %d targets\n", t_nufft, length(near_idx))
@printf("  FMM share of the two: %.1f%%\n", 100*t_fmm/(t_fmm + t_nufft))

# --- pottrg: the corrected layer-potential map over ALL targets. NOT part of the two calls
# above, and by elimination the dominant term: lfmm3d + NUFFT come to ~2 s against a measured
# eval of 21.6 s (64 threads) and 99 s (96). It is FMM3D over the interface PLUS adaptive
# hcubature corrections for every near target-panel pair, and the correction count scales with
# the target count -- 52x more targets here than the row protocol's 136,580.
const P = (n_quad = 6, edge_level = 4, rhs_tol = 1e-3, lhs_tol = 1e-5, max_depth = 12)
env = BI.envelope_volume_source(b)
t_itf = @elapsed interface = BI.multi_dielectric_box3d_rhs_adaptive(
    P.n_quad, tz / 2.0^P.edge_level * 1.01, boxes, Float64[11.9, 3.9, 2.4], env, P.rhs_tol;
    eps_out = 1.0, max_depth = P.max_depth, tkm_kmax = BI._estimate_tkm3dc_kmax(env))
@printf("  interface build                 %8.2f s   N = %d\n", t_itf, BI.num_points(interface))
t_pot = @elapsed pottrg = BI.laplace3d_pottrg_fmm3d_corrected_hcubature(
    interface, T, P.lhs_tol, P.lhs_tol, 5.0)
@printf("  pottrg build (FMM + hcubature)  %8.2f s   over ALL %d targets\n", t_pot, size(T,2))
sig = randn(BI.num_points(interface))
t_mv = @elapsed pottrg * sig
@printf("  pottrg matvec (per column)      %8.2f s\n", t_mv)
@printf("  ==> accounted: fmm %.1f + nufft %.1f + pottrg %.1f + matvec %.1f = %.1f s\n",
        t_fmm, t_nufft, t_pot, t_mv, t_fmm + t_nufft + t_pot + t_mv)

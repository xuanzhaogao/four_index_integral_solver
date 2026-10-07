# Before/after for the outer-loop threading fix in laplace3d_pottrg_corrections_hcubature.
#
# Times the pottrg build -- the corrected layer-potential map over the fixed 198-orbital target
# set (N_p = 7,168,390) -- which isolation measurements pinned as the source of the evaluation's
# 96-thread collapse: lfmm3d (0.95 s) and the NUFFT path (1.25 s) were both flat across thread
# counts, leaving pottrg to account for the 21.6 s (64 threads) vs 99.1 s (96) difference.
#
# Run the SAME script against each BI checkout, one thread count per process:
#   --project=/mnt/home/xgao1/codes/BoundaryIntegral.jl                    (live, unfixed)
#   --project=/mnt/home/xgao1/codes/BoundaryIntegral-hcub-outer-threads    (worktree, fixed)
#
# No Printf: it is not a dependency of BI's own project, and these run under that project so
# that the worktree's Manifest (with its absolute dev paths) resolves TKM3D/FINUFFT.
using BoundaryIntegral, Serialization
import BoundaryIntegral as BI

const REF = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
TEMPL = [BI.read_xsf(joinpath(REF, "graphene_0000$(i).xsf"))[2] for i in 1:2]
CEN = [ntuple(d -> Float64(BI.density_centroid(dg)[d]), 3) for dg in TEMPL]
b = BI.assemble_lattice_batch(TEMPL, Dict(1 => BI.OrbitalInstance(1, 1,
        BI.snap_orbital(TEMPL[1], CEN[1], CEN[1]))), [(1,1)]; support_rtol = 1e-4)
L = 270.0; tz = 9.0; cx, cy, cz = CEN[1]; h = L/2; czi = cz - tz/2 - h
boxes = BI.BoxGeom[(center=(cx-h,cy,czi),Lx=L,Ly=L,Lz=L),(center=(cx+h,cy,czi),Lx=L,Ly=L,Lz=L),
                   (center=(cx,cy,cz),Lx=42.95357656,Ly=44.06001031,Lz=tz)]
env = BI.envelope_volume_source(b)
itf = BI.multi_dielectric_box3d_rhs_adaptive(6, tz/2.0^4*1.01, boxes, Float64[11.9,3.9,2.4],
        env, 1e-3; eps_out = 1.0, max_depth = 12, tkm_kmax = BI._estimate_tkm3dc_kmax(env))

tc = deserialize("/mnt/home/xgao1/work/four_index_integral_solver/codes/article/" *
                 "sec52_performance/data_sec53evalpilot/targets_model.jls")
T = tc.positions

println("threads=", Threads.nthreads(), "  OMP=", get(ENV,"OMP_NUM_THREADS","unset"),
        "  N=", BI.num_points(itf), "  N_p=", size(T,2))
# warm on a small slice so the timed call is not paying compilation
BI.laplace3d_pottrg_fmm3d_corrected_hcubature(itf, T[:, 1:2000], 1e-5, 1e-5, 5.0)
t = @elapsed BI.laplace3d_pottrg_fmm3d_corrected_hcubature(itf, T, 1e-5, 1e-5, 5.0)
println("POTTRG_BUILD_SECONDS ", round(t; digits = 2))

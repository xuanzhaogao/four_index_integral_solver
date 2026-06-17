# Decisive lever test: time the multi/envelope build at a given thread config.
# If the build time is ~unchanged when FINUFFT(OMP)/BLAS threads are capped, the cost
# is the SERIAL Julia adaptive-refinement work (resolution-check interpolation loop),
# NOT the NUFFT/FMM transforms or any threaded library. Build only, no solve.
#
# Run e.g.:  OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 julia -t 16 --project=<BI> diag_build_times.jl

using BoundaryIntegral
import BoundaryIntegral as BI
using Printf, LinearAlgebra

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const TEMPL = [BI.read_xsf(joinpath(REF_DIR, "graphene_00001.xsf"))[2],
               BI.read_xsf(joinpath(REF_DIR, "graphene_00002.xsf"))[2]]
const CENTROID = [ntuple(d -> Float64(BI.density_centroid(dg)[d]), 3) for dg in TEMPL]
const A1 = (2.465, 0.0, 0.0)
const A2 = (-1.2325, 2.1347526, 0.0)
addv(a, b) = (a[1]+b[1], a[2]+b[2], a[3]+b[3])
scalev(s, a) = (s*a[1], s*a[2], s*a[3])
distv(a, b) = sqrt((a[1]-b[1])^2 + (a[2]-b[2])^2 + (a[3]-b[3])^2)

const LZ, EPS_IN, EPS_OUT, L = 3.35, 3.5, 1.0, 90.0
const N_QUAD, EDGE, RHS_TOL, MAXDEPTH, SUPPORT_RTOL = 6, 4, 1e-3, 12, 1e-4
const L_EC = LZ / 2.0^EDGE * 1.01
const BOXES = BI.BoxGeom[(center = (0.0, 0.0, CENTROID[1][3]), Lx = L, Ly = L, Lz = LZ)]
const EPSES = Float64[EPS_IN]

function batch(cutoff; nrange = 3)
    center = CENTROID[1]
    sites = Tuple{Int,NTuple{3,Float64}}[(1, center)]
    for t in 1:2, n1 in -nrange:nrange, n2 in -nrange:nrange
        (t == 1 && n1 == 0 && n2 == 0) && continue
        pos = addv(CENTROID[t], addv(scalev(Float64(n1), A1), scalev(Float64(n2), A2)))
        distv(pos, center) <= cutoff && push!(sites, (t, pos))
    end
    insts = Dict(id => BI.OrbitalInstance(id, t, BI.snap_orbital(TEMPL[t], CENTROID[t], pos))
                 for (id, (t, pos)) in enumerate(sites))
    return BI.assemble_lattice_batch(TEMPL, insts, [(1, id) for id in 1:length(sites)];
                                     support_rtol = SUPPORT_RTOL)
end

function build(cutoff)
    b = batch(cutoff)
    env = BI.envelope_volume_source(b)
    kmax = BI._estimate_tkm3dc_kmax(env)
    return BI.multi_dielectric_box3d_rhs_adaptive(N_QUAD, L_EC, BOXES, EPSES, env, RHS_TOL;
        eps_out = EPS_OUT, max_depth = MAXDEPTH, tkm_kmax = kmax)
end

@printf("julia -t %d   OMP=%s   OPENBLAS=%s   BLAS.get_num_threads=%d\n",
        Threads.nthreads(), get(ENV, "OMP_NUM_THREADS", "?"),
        get(ENV, "OPENBLAS_NUM_THREADS", "?"), BLAS.get_num_threads())
build(0.0)  # warm-up / compile
for cutoff in (1.5, 2.5)
    GC.gc()
    t = @elapsed itf = build(cutoff)
    @printf("cutoff %.1f  build %.1f s  pts %d\n", cutoff, t, BI.num_points(itf))
end
println("DIAG BUILD TIMES DONE")

# 6.5 multi-RHS — walltime + peak-RAM benchmark of the multi-RHS production pipeline on
# the REAL graphene monolayer: one central Wannier orbital (sublattice A) plus its
# lattice neighbors, paired into K densities rho = phi_center * phi_neighbor solved on
# ONE shared interface. K grows with the neighbor cutoff (onsite -> nn -> nnn shells).
#
# Mirrors run_single_rhs.jl's instrumentation (warm-up, inline per-stage walltime +
# Sys.maxrss snapshots, smoke mode). Drives the BI.jl lattice multi-RHS primitives
# directly so each stage is timed; builds the batched operator EXPLICITLY with
# correct_edges=true / max_order=64 (production), not the solve_dielectric_box3d_block
# defaults (correct_edges=false / max_order=8). Eval uses four_index_matrix.
#
# Geometry: lattice_scale/demo_2x2 convention (box center cz = A-centroid z ~ 7.5,
# orbital at box MIDPLANE) — differs from run_single_rhs.jl (top face); see spec caveat.
#
# Env knobs:
#   ORBBENCH_SMOKE=1     coarse tolerances + reduced cutoff list
#   CORRECT_EDGES=0      disable edge correction (default on, production)
#   ORBBENCH_CUTOFFS     comma-separated cutoff override, e.g. "0.0,1.5,2.5"
#   ORBBENCH_GEOM_ONLY=1 print K per cutoff and exit (no solve)
#
# Run (production, cluster): JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 julia --project=<BI> run_multi_rhs.jl
# Smoke (dev):               ORBBENCH_SMOKE=1 julia -t 8 --project=<BI> run_multi_rhs.jl

using BoundaryIntegral
import BoundaryIntegral as BI
using Krylov
using LinearAlgebra
using Printf
using Serialization

const GRAPHENE = normpath(joinpath(@__DIR__, "..", "..", "..", "graphene"))
include(joinpath(GRAPHENE, "bilayer_slab", "src", "ScreenedOrbitalSolve.jl"))
using .ScreenedOrbitalSolve

const SMOKE = get(ENV, "ORBBENCH_SMOKE", "0") == "1"
const CORRECT_EDGES = get(ENV, "CORRECT_EDGES", "1") == "1"
const GEOM_ONLY = get(ENV, "ORBBENCH_GEOM_ONLY", "0") == "1"

const DATA = joinpath(@__DIR__, "..", "data", SMOKE ? "smoke" : "")
mkpath(joinpath(DATA, "raw"))

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")   # sublattice A (type 1)
const XSF_2 = joinpath(REF_DIR, "graphene_00002.xsf")   # sublattice B (type 2)

# graphene hexagonal lattice (bohr), demo_2x2 convention
const A1 = (2.465, 0.0, 0.0)
const A2 = (-1.2325, 2.1347526, 0.0)

const LZ, EPS_IN, EPS_OUT, L = 3.35, 3.5, 1.0, 90.0
const P = SMOKE ?
    (n_quad = 6, edge_level = 2, rhs_tol = 1e-2, lhs_tol = 1e-3,
     gmres_atol = 1e-3, gmres_rtol = 1e-3, max_order = 64, max_depth = 12,
     support_rtol = 1e-3) :
    (n_quad = 6, edge_level = 4, rhs_tol = 1e-3, lhs_tol = 1e-5,
     gmres_atol = 1e-5, gmres_rtol = 1e-5, max_order = 64, max_depth = 12,
     support_rtol = 1e-4)
const VOLUME_TOL = P.rhs_tol

const CUTOFFS = let env = get(ENV, "ORBBENCH_CUTOFFS", "")
    !isempty(env) ? parse.(Float64, split(env, ",")) :
        (SMOKE ? [0.0, 1.5] : [0.0, 1.5, 2.5, 2.9, 3.9])
end

# raw phi datagrids (sublattice A, B) and their Cartesian density centroids
const TEMPL = [BI.read_xsf(XSF_1)[2], BI.read_xsf(XSF_2)[2]]
const CENTROID = [ntuple(d -> Float64(BI.density_centroid(dg)[d]), 3) for dg in TEMPL]

addv(a, b) = (a[1] + b[1], a[2] + b[2], a[3] + b[3])
scalev(s, a) = (s * a[1], s * a[2], s * a[3])
distv(a, b) = sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2 + (a[3] - b[3])^2)

gb(x) = x / 2^30
rss_gb() = gb(Sys.maxrss())
live_gb() = gb(Base.gc_live_bytes())

# Build the K pair densities rho = phi_center * phi_neighbor for the central A orbital
# (id 1) and every A/B orbital within `cutoff` of it. Returns (LatticeBatch, pairs).
function build_geometry(cutoff::Float64; nrange::Int = 3)
    center = CENTROID[1]                          # sublattice A, cell (0,0)
    sites = Tuple{Int, NTuple{3, Float64}}[(1, center)]   # (template type, pos); center first
    for t in 1:2, n1 in -nrange:nrange, n2 in -nrange:nrange
        (t == 1 && n1 == 0 && n2 == 0) && continue        # skip center duplicate
        pos = addv(CENTROID[t], addv(scalev(Float64(n1), A1), scalev(Float64(n2), A2)))
        distv(pos, center) <= cutoff && push!(sites, (t, pos))
    end
    insts = Dict{Int, BI.OrbitalInstance}()
    for (id, (t, pos)) in enumerate(sites)
        steps = BI.snap_orbital(TEMPL[t], CENTROID[t], pos)
        insts[id] = BI.OrbitalInstance(id, t, steps)
    end
    pairs = [(1, id) for id in 1:length(sites)]           # onsite (1,1) + (1,neighbor)
    b = BI.assemble_lattice_batch(TEMPL, insts, pairs; support_rtol = P.support_rtol)
    return b, pairs
end

if GEOM_ONLY
    @printf("smoke=%s  centroid_A=%s  centroid_B=%s\n", SMOKE, CENTROID[1], CENTROID[2])
    for cutoff in CUTOFFS
        b, pairs = build_geometry(cutoff)
        @printf("cutoff %.2f  ->  K = %2d pairs   support points = %d\n",
                cutoff, length(pairs), size(b.densities, 1))
    end
    println("GEOM ONLY DONE")
    exit(0)
end

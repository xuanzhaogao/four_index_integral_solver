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

# Per-K instrumented pipeline. `stage!` records (label, t, maxrss, live) when record=true.
function run_pipeline(b, pairs, stages::Vector; record::Bool)
    note = record ? (l, t) -> begin
        push!(stages, (; label = l, t, maxrss_gb = rss_gb(), live_gb = live_gb()))
        @printf("    %-34s %9.2f s   maxrss %6.2f GB   live %6.2f GB\n", l, t, rss_gb(), live_gb())
        flush(stdout)
    end : (l, t) -> nothing

    K = length(pairs)
    l_ec = LZ / 2.0^P.edge_level * 1.01
    boxes = BI.BoxGeom[(center = (0.0, 0.0, CENTROID[1][3]), Lx = L, Ly = L, Lz = LZ)]
    epses = Float64[EPS_IN]

    # --- precompute (RHS-independent across the K columns) ---
    local env, kmax
    t = @elapsed begin
        env = BI.envelope_volume_source(b)
        kmax = BI._estimate_tkm3dc_kmax(env)
    end
    note("envelope + tkm kmax", t)
    local interface
    t = @elapsed interface = BI.multi_dielectric_box3d_rhs_adaptive(
        P.n_quad, l_ec, boxes, epses, env, P.rhs_tol;
        eps_out = EPS_OUT, max_depth = P.max_depth, tkm_kmax = kmax)
    note("interface build (envelope)", t)
    local op
    t = @elapsed op = BI.batched_lhs_dielectric_box3d_fmm3d_corrected(
        interface, P.lhs_tol, P.lhs_tol, P.max_order; correct_edges = CORRECT_EDGES)
    note("batched LHS operator", t)

    sources = BI.batch_volume_sources(b)

    # --- solve (scales with K) ---
    local F
    t = @elapsed F = BI.rhs_dielectric_box3d_fmm3d(interface, sources, P.rhs_tol)
    note("RHS assembly (batched nd=K)", t)
    local sigma_block, bstats
    t = @elapsed begin
        sigma_block, bstats = Krylov.block_gmres(op, F;
            rtol = P.gmres_rtol, atol = P.gmres_atol, itmax = 500)
    end
    note("block GMRES", t)
    block_resid = norm(op * sigma_block - F) / max(norm(F), eps(Float64))

    # --- sequential baseline: same interface, single-RHS operator built ONCE, K gmres ---
    local lhs_seq
    t_seq_op = @elapsed lhs_seq = BI.lhs_dielectric_box3d_fmm3d_corrected(
        interface, P.lhs_tol, P.lhs_tol, P.max_order; correct_edges = CORRECT_EDGES)
    sigma_seq = Matrix{Float64}(undef, size(sigma_block)...)
    t_seq_solve = @elapsed for k in 1:K
        rhs_k = BI.rhs_dielectric_box3d_fmm3d(interface, sources[k], P.rhs_tol)
        sk, _ = Krylov.gmres(lhs_seq, rhs_k; atol = P.gmres_atol, rtol = P.gmres_rtol)
        sigma_seq[:, k] = sk
    end
    note("sequential LHS operator", t_seq_op)
    note("sequential K gmres", t_seq_solve)
    seq_agree = norm(sigma_block - sigma_seq) / max(norm(sigma_block), eps(Float64))

    # --- eval + contraction into V[a,b] (TKM u_inc; 96-thread pathology by design) ---
    local V
    t = @elapsed V = BI.four_index_matrix(interface, sources, sigma_block;
        lhs_tol = P.lhs_tol, volume_tol = VOLUME_TOL, range_factor = 5.0)
    note("eval + contract (four_index_matrix)", t)

    scaleV = maximum(abs.(V))
    max_rel_asym = scaleV > 0 ? maximum(abs.(V - V')) / scaleV : 0.0
    onsite_N = sum(sources[1].weights .* sources[1].density)
    u_onsite_ev = to_eV(V[1, 1], onsite_N, onsite_N)

    return (; K, n_points = BI.num_points(interface), n_src = size(b.densities, 1),
            niter = bstats.niter, block_resid, seq_agree, max_rel_asym,
            v11_raw = V[1, 1], u_onsite_ev, V, pairs)
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

@printf("threads = %d   correct_edges = %s   smoke = %s   cutoffs = %s\n",
        Threads.nthreads(), CORRECT_EDGES, SMOKE, CUTOFFS)
flush(stdout)

# warm-up: compile every stage on the smallest (K=1) real batch, un-recorded.
println(">>> warm-up (K=1 batch through every stage)"); flush(stdout)
let (bw, pw) = build_geometry(0.0)
    tw = @elapsed run_pipeline(bw, pw, NamedTuple[]; record = false)
    @printf("  warm-up: %.1f s   baseline maxrss %.2f GB\n", tw, rss_gb())
end
const RSS_BASELINE = rss_gb()

const RESULTS = NamedTuple[]
for cutoff in CUTOFFS
    println("=" ^ 72)
    @printf(">>> cutoff = %.2f bohr  (L=%g Lz=%g eps_in=%g)\n", cutoff, L, LZ, EPS_IN)
    flush(stdout)
    b, pairs = build_geometry(cutoff)
    stages = NamedTuple[]
    t_total = @elapsed res = run_pipeline(b, pairs, stages; record = true)

    tof(lbl) = first(s.t for s in stages if s.label == lbl)
    t_precompute = tof("envelope + tkm kmax") + tof("interface build (envelope)") +
                   tof("batched LHS operator")
    t_solve_block = tof("RHS assembly (batched nd=K)") + tof("block GMRES")
    t_seq_op = tof("sequential LHS operator")
    # end-to-end solve cost vs t_solve_block; both exclude the 1-time LHS op build.
    # t_solve_seq = K RHS builds + K gmres;  t_solve_block = 1 batched RHS + block gmres
    # (so the speedup reflects both RHS batching and block GMRES).
    t_solve_seq = tof("sequential K gmres")
    t_eval = tof("eval + contract (four_index_matrix)")

    @printf("  K=%d  interface points %d  src %d  niter %d  block_resid %.2e\n",
            res.K, res.n_points, res.n_src, res.niter, res.block_resid)
    @printf("  seq_agree %.2e  max_rel_asym %.2e  u_onsite = %.4f eV\n",
            res.seq_agree, res.max_rel_asym, res.u_onsite_ev)
    @printf("  precompute %.1f s | block solve %.1f s (vs seq %.1f s, %.2fx) | eval %.1f s | total %.1f s\n",
            t_precompute, t_solve_block, t_solve_seq,
            t_solve_seq / max(t_solve_block, eps()), t_eval, t_total)
    @printf("  per-RHS marginal (block solve+eval)/K = %.1f s   peak RSS %.2f GB\n",
            (t_solve_block + t_eval) / res.K, rss_gb())

    out = (; smoke = SMOKE, correct_edges = CORRECT_EDGES, cutoff,
           L, Lz = LZ, eps_in = EPS_IN, eps_out = EPS_OUT,
           K = res.K, pairs = res.pairs, n_points = res.n_points, n_src = res.n_src,
           niter = res.niter, block_resid = res.block_resid, seq_agree = res.seq_agree,
           max_rel_asym = res.max_rel_asym, v11_raw = res.v11_raw,
           u_onsite_ev = res.u_onsite_ev, V = res.V,
           t_precompute, t_solve_block, t_seq_op, t_solve_seq, t_eval, t_total,
           stages = copy(stages), rss_baseline_gb = RSS_BASELINE, rss_peak_gb = rss_gb(),
           nthreads = Threads.nthreads(), hostname = gethostname())
    serialize(joinpath(DATA, "raw", "multi_rhs_K$(res.K)_cut$(cutoff).jls"), out)
    push!(RESULTS, out)
end

# summary table + CSV
println("=" ^ 72)
@printf("%-4s %-8s %-9s %-7s %-7s %-7s %-7s %-7s %-9s %-9s\n",
        "K", "n_pts", "precomp", "blk", "seq", "spdup", "eval", "tot", "perRHS", "peakGB")
for r in RESULTS
    @printf("%-4d %-8d %-9.1f %-7.1f %-7.1f %-7.2f %-7.1f %-7.1f %-9.1f %-9.2f\n",
            r.K, r.n_points, r.t_precompute, r.t_solve_block, r.t_solve_seq,
            r.t_solve_seq / max(r.t_solve_block, eps()), r.t_eval, r.t_total,
            (r.t_solve_block + r.t_eval) / r.K, r.rss_peak_gb)
end
println("=" ^ 72)

let csv = joinpath(DATA, "multi_rhs.csv")
    newfile = !isfile(csv)
    open(csv, "a") do io
        newfile && println(io, join(["hostname", "nthreads", "smoke", "correct_edges",
            "cutoff", "K", "n_points", "n_src", "niter", "block_resid", "seq_agree",
            "max_rel_asym", "t_precompute", "t_solve_block", "t_seq_op", "t_solve_seq",
            "t_eval", "t_total", "rss_peak_gb", "u_onsite_ev"], ","))
        for r in RESULTS
            println(io, join(string.([r.hostname, r.nthreads, r.smoke, r.correct_edges,
                r.cutoff, r.K, r.n_points, r.n_src, r.niter, r.block_resid, r.seq_agree,
                r.max_rel_asym, r.t_precompute, r.t_solve_block, r.t_seq_op, r.t_solve_seq,
                r.t_eval, r.t_total, r.rss_peak_gb, r.u_onsite_ev]), ","))
        end
    end
    println("saved CSV -> $csv")
end
println("MULTI RHS BENCH DONE")

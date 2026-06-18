# 6.6 multi-cube multi-RHS — the exp65 multi-RHS pipeline on the two-cube heterojunction.
# One central graphene Wannier p_z orbital (sublattice A) plus its lattice neighbors,
# paired into K densities rho = phi_center * phi_neighbor, solved on ONE shared interface.
# K grows with the neighbor cutoff (onsite -> nn -> nnn shells). The SUBSTRATE is the
# multicube: two L x L x L cubes (eps1 | eps2, Si/SiO2) sharing the buried face x = c_x,
# with an L x L x (L/10) slab (eps_slab) on top, anchored on the central orbital
# (slab center = orbital centroid, junction directly beneath it).
#
# Mirrors exp65/run_multi_rhs.jl exactly (build_geometry, batched LHS, block GMRES, eval
# of the central row V[rho_11, rho_b]); only the boxes/epses differ (multicube vs single
# graphene slab). Runs at 96 threads: TKM3D's in-function FINUFFT nthreads cap (16) avoids
# the FFTW 96-thread plan pathology, so no DUCC override is needed.
#
# Env knobs:
#   MULTICUBE_CUTOFF     single cutoff per process -> one Slurm-array task per K
#   MULTICUBE_SMOKE=1    coarse tolerances + reduced cutoff list
#   MULTICUBE_CUTOFFS    comma-separated cutoff override, e.g. "0.0,1.5,2.5"
#   CORRECT_EDGES=0      disable edge correction (default on, production)
#   MULTICUBE_GEOM_ONLY=1  print K per cutoff and exit (no solve)
#
# Run (one cutoff): MULTICUBE_CUTOFF=1.5 JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 \
#                     julia --project=. exp66_multicube/scripts/run_multicube.jl
# Smoke (dev):      MULTICUBE_SMOKE=1 MULTICUBE_CUTOFF=1.5 julia -t 8 --project=. ... run_multicube.jl

using BoundaryIntegral
import BoundaryIntegral as BI
using Krylov
using LinearAlgebra
using Printf
using Serialization

const SMOKE = get(ENV, "MULTICUBE_SMOKE", "0") == "1"
const CORRECT_EDGES = get(ENV, "CORRECT_EDGES", "1") == "1"
const GEOM_ONLY = get(ENV, "MULTICUBE_GEOM_ONLY", "0") == "1"

const DATA = joinpath(@__DIR__, "..", "data", SMOKE ? "smoke" : "")
mkpath(joinpath(DATA, "raw"))

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")   # sublattice A (type 1)
const XSF_2 = joinpath(REF_DIR, "graphene_00002.xsf")   # sublattice B (type 2)

# graphene hexagonal lattice (bohr), demo_2x2 convention (as exp65)
const A1 = (2.465, 0.0, 0.0)
const A2 = (-1.2325, 2.1347526, 0.0)

# multicube substrate: two L x L x L cubes (eps1 | eps2) sharing x = c_x, slab on top.
const L = 90.0
const SLAB_THICK = L / 10            # 9: thick enough to fully contain the orbital plane (cf. smoke)
const EPS1, EPS2, EPS_SLAB, EPS_OUT = 11.9, 3.9, 10.0, 1.0

const P = SMOKE ?
    (n_quad = 6, edge_level = 2, rhs_tol = 1e-2, lhs_tol = 1e-3,
     gmres_atol = 1e-3, gmres_rtol = 1e-3, max_order = 64, max_depth = 12,
     support_rtol = 1e-3) :
    (n_quad = 6, edge_level = 4, rhs_tol = 1e-3, lhs_tol = 1e-5,
     gmres_atol = 1e-5, gmres_rtol = 1e-5, max_order = 64, max_depth = 12,
     support_rtol = 1e-4)
const VOLUME_TOL = P.rhs_tol

const CUTOFFS = let env = get(ENV, "MULTICUBE_CUTOFFS", "")
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

# K pair densities rho = phi_center * phi_neighbor for the central A orbital (id 1) and
# every A/B orbital within `cutoff` of it. Identical to exp65. Returns (LatticeBatch, pairs).
function build_geometry(cutoff::Float64; nrange::Int = 3)
    center = CENTROID[1]                          # sublattice A, cell (0,0)
    sites = Tuple{Int, NTuple{3, Float64}}[(1, center)]
    for t in 1:2, n1 in -nrange:nrange, n2 in -nrange:nrange
        (t == 1 && n1 == 0 && n2 == 0) && continue
        pos = addv(CENTROID[t], addv(scalev(Float64(n1), A1), scalev(Float64(n2), A2)))
        distv(pos, center) <= cutoff && push!(sites, (t, pos))
    end
    insts = Dict{Int, BI.OrbitalInstance}()
    for (id, (t, pos)) in enumerate(sites)
        steps = BI.snap_orbital(TEMPL[t], CENTROID[t], pos)
        insts[id] = BI.OrbitalInstance(id, t, steps)
    end
    pairs = [(1, id) for id in 1:length(sites)]
    b = BI.assemble_lattice_batch(TEMPL, insts, pairs; support_rtol = P.support_rtol)
    return b, pairs
end

# multicube substrate boxes, anchored on the central orbital (slab center = centroid).
function multicube_boxes()
    cx, cy, cz = CENTROID[1]
    h = L / 2; tz = SLAB_THICK; czi = cz - tz / 2 - h     # cube tops at the slab bottom
    boxes = BI.BoxGeom[
        (center = (cx - h, cy, czi), Lx = L, Ly = L, Lz = L),   # Omega_1, x in [cx-L, cx]  (eps1)
        (center = (cx + h, cy, czi), Lx = L, Ly = L, Lz = L),   # Omega_2, x in [cx, cx+L]  (eps2)
        (center = (cx, cy, cz),      Lx = L, Ly = L, Lz = tz)]  # slab over the junction    (eps_slab)
    return boxes, Float64[EPS1, EPS2, EPS_SLAB], tz
end

function run_pipeline(b, pairs, stages::Vector; record::Bool)
    note = record ? (l, t) -> begin
        push!(stages, (; label = l, t, maxrss_gb = rss_gb(), live_gb = live_gb()))
        @printf("    %-34s %9.2f s   maxrss %6.2f GB   live %6.2f GB\n", l, t, rss_gb(), live_gb())
        flush(stdout)
    end : (l, t) -> nothing

    K = length(pairs)
    boxes, epses, tz = multicube_boxes()
    l_ec = tz / 2.0^P.edge_level * 1.01

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
    # Box-based multi-region screening: the interface spans 3 eps regions, so the
    # interface-based screened_volume_source exp65 used (uniform eps_in only) can't apply
    # here — screen each source by which box it sits in (all orbitals lie in the slab).
    screened = [BI.screened_volume_source(boxes, epses, EPS_OUT, sources[k], BI.SharpScreening()) for k in 1:K]

    # --- solve (scales with K) ---
    # RHS per source (the batched rhs_dielectric_box3d_fmm3d re-screens via the interface;
    # for >1 distinct-position source it loops the single-source path anyway, so we loop the
    # 4-arg single-source form with eps_src=1 on the already-screened sources).
    local F
    t = @elapsed begin
        F = Matrix{Float64}(undef, BI.num_points(interface), K)
        for k in 1:K
            F[:, k] = BI.rhs_dielectric_box3d_fmm3d(interface, screened[k], 1.0, P.rhs_tol)
        end
    end
    note("RHS assembly (per-source, multi-region)", t)
    local sigma_block, bstats
    t = @elapsed begin
        sigma_block, bstats = Krylov.block_gmres(op, F;
            rtol = P.gmres_rtol, atol = P.gmres_atol, itmax = 500)
    end
    note("block GMRES", t)
    block_resid = norm(op * sigma_block - F) / max(norm(F), eps(Float64))

    # --- eval: central row V[rho_11, rho_b], target = onsite rho_11 = phi_1^2 (exp65 protocol).
    onsite = sources[1]
    nz = findall(>=(P.support_rtol * maximum(abs, onsite.density)), abs.(onsite.density))
    tgt = Matrix{Float64}(onsite.positions[:, nz])
    tw = onsite.weights[nz] .* onsite.density[nz]
    # eval-side precompute: corrected layer-potential map over the fixed rho_11 targets.
    local pottrg
    t = @elapsed pottrg = BI.laplace3d_pottrg_fmm3d_corrected_hcubature(interface, tgt, P.lhs_tol, P.lhs_tol, 5.0)
    note("eval pottrg build", t)
    local Vrow
    t = @elapsed begin
        Vrow = Vector{Float64}(undef, K)
        for bcol in 1:K
            sb = screened[bcol]                      # box-based multi-region screened source
            vals = BI.TKM3D.ltkm3dc(VOLUME_TOL, sb.positions; charges = sb.weights .* sb.density,
                                    targets = tgt, pgt = 1, kmax = BI._estimate_tkm3dc_kmax(sb))
            Vrow[bcol] = dot(tw, real.(vals.pottarg) .+ (pottrg * sigma_block[:, bcol]))
        end
    end
    note("eval (onsite-row, K field evals)", t)

    # bare (interface-free) onsite self-energy reference: int rho_11 * TKM[rho_11 unscreened]
    q1 = onsite.weights .* onsite.density
    vac = BI.TKM3D.ltkm3dc(VOLUME_TOL, onsite.positions; charges = q1, targets = tgt,
                           pgt = 1, kmax = BI._estimate_tkm3dc_kmax(onsite))
    v11_vac = dot(tw, real.(vac.pottarg))

    return (; K, n_points = BI.num_points(interface), n_src = size(b.densities, 1),
            niter = bstats.niter, block_resid,
            v11_raw = Vrow[1], v11_vac, screen_ratio = v11_vac / Vrow[1], Vrow, pairs)
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

# Single cutoff per invocation (MULTICUBE_CUTOFF) -> one Slurm-array task per cutoff/K.
const ONE_CUTOFF = let v = get(ENV, "MULTICUBE_CUTOFF", ""); isempty(v) ? nothing : parse(Float64, v) end
const SWEEP = ONE_CUTOFF === nothing ? CUTOFFS : [ONE_CUTOFF]

@printf("threads = %d   correct_edges = %s   smoke = %s   cutoffs = %s\n",
        Threads.nthreads(), CORRECT_EDGES, SMOKE, SWEEP)
@printf("substrate: two %gx%gx%g cubes eps %g|%g, slab %gx%gx%g eps %g, eps_out %g\n",
        L, L, L, EPS1, EPS2, L, L, SLAB_THICK, EPS_SLAB, EPS_OUT)
flush(stdout)

# warm-up: compile every stage on the smallest (K=1) real batch, un-recorded.
println(">>> warm-up (K=1 batch through every stage)"); flush(stdout)
let (bw, pw) = build_geometry(0.0)
    tw = @elapsed run_pipeline(bw, pw, NamedTuple[]; record = false)
    @printf("  warm-up: %.1f s   baseline maxrss %.2f GB\n", tw, rss_gb())
end
const RSS_BASELINE = rss_gb()

const RESULTS = NamedTuple[]
for cutoff in SWEEP
    println("=" ^ 72)
    @printf(">>> cutoff = %.2f bohr\n", cutoff); flush(stdout)
    b, pairs = build_geometry(cutoff)
    stages = NamedTuple[]
    t_total = @elapsed res = run_pipeline(b, pairs, stages; record = true)

    tof(lbl) = first(s.t for s in stages if s.label == lbl)
    t_precompute = tof("envelope + tkm kmax") + tof("interface build (envelope)") +
                   tof("batched LHS operator")
    t_solve_block = tof("RHS assembly (per-source, multi-region)") + tof("block GMRES")
    t_pottrg = tof("eval pottrg build")
    t_eval = tof("eval (onsite-row, K field evals)")

    @printf("  K=%d  interface points %d  src %d  niter %d  block_resid %.2e\n",
            res.K, res.n_points, res.n_src, res.niter, res.block_resid)
    @printf("  V[11] = %.6e   V_vac[11] = %.6e   V_vac/V = %.4f\n",
            res.v11_raw, res.v11_vac, res.screen_ratio)
    @printf("  precompute %.1f s | pottrg %.1f s | block solve %.1f s | eval %.1f s | total %.1f s   peak RSS %.2f GB\n",
            t_precompute, t_pottrg, t_solve_block, t_eval, t_total, rss_gb())
    flush(stdout)

    out = (; smoke = SMOKE, correct_edges = CORRECT_EDGES, cutoff,
           L, slab_thick = SLAB_THICK, eps1 = EPS1, eps2 = EPS2, eps_slab = EPS_SLAB, eps_out = EPS_OUT,
           K = res.K, pairs = res.pairs, n_points = res.n_points, n_src = res.n_src,
           niter = res.niter, block_resid = res.block_resid,
           v11_raw = res.v11_raw, v11_vac = res.v11_vac, screen_ratio = res.screen_ratio, Vrow = res.Vrow,
           t_precompute, t_pottrg, t_solve_block, t_eval, t_total,
           stages = copy(stages), rss_baseline_gb = RSS_BASELINE, rss_peak_gb = rss_gb(),
           nthreads = Threads.nthreads(), hostname = gethostname())
    serialize(joinpath(DATA, "raw", "multicube_K$(res.K)_cut$(cutoff).jls"), out)
    push!(RESULTS, out)
end

# Full-sweep mode (single process): print summary + write the combined CSV. In per-cutoff
# (Slurm-array) mode each task writes only its .jls; gather them afterward.
if ONE_CUTOFF === nothing
    println("=" ^ 72)
    @printf("%-4s %-9s %-12s %-8s %-9s %-9s %-9s %-9s\n",
            "K", "n_pts", "V[11]", "Vvac/V", "precomp", "pottrg", "block", "eval")
    for r in RESULTS
        @printf("%-4d %-9d %-12.5e %-8.3f %-9.1f %-8.1f %-9.1f %-9.1f\n",
                r.K, r.n_points, r.v11_raw, r.screen_ratio, r.t_precompute, r.t_pottrg,
                r.t_solve_block, r.t_eval)
    end
    println("=" ^ 72)
    open(joinpath(DATA, "multicube.csv"), "w") do io
        println(io, join(["hostname", "nthreads", "smoke", "correct_edges", "cutoff",
            "K", "n_points", "n_src", "niter", "block_resid", "v11_raw", "v11_vac", "screen_ratio",
            "t_precompute", "t_pottrg", "t_solve_block", "t_eval", "t_total", "rss_peak_gb"], ","))
        for r in RESULTS
            println(io, join(string.([r.hostname, r.nthreads, r.smoke, r.correct_edges,
                r.cutoff, r.K, r.n_points, r.n_src, r.niter, r.block_resid,
                r.v11_raw, r.v11_vac, r.screen_ratio,
                r.t_precompute, r.t_pottrg, r.t_solve_block, r.t_eval, r.t_total, r.rss_peak_gb]), ","))
        end
    end
    println("saved CSV -> ", joinpath(DATA, "multicube.csv"))
end
println("MULTICUBE MULTI-RHS DONE")

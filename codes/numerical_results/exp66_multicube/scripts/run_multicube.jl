# 6.6 multi-cube multi-RHS — the exp65 multi-RHS pipeline on the two-cube heterojunction.
# One central graphene Wannier p_z orbital (sublattice A) plus its lattice neighbors,
# paired into K densities rho = phi_center * phi_neighbor, solved on ONE shared interface.
# K grows with the neighbor cutoff (onsite -> nn -> nnn shells). The SUBSTRATE is the
# multicube: two L x L x L cubes (eps1 | eps2, Si/SiO2) sharing the buried face x = c_x,
# with a slab (eps_slab, 9 A thick) on top, anchored on the central orbital
# (slab center = orbital centroid, junction directly beneath it). By default this is the
# SAME system the article's lattice section (Sec. 5.3) reports: L = 270 and eps_slab = 2.4
# (see MULTICUBE_GEOM below; MULTICUBE_GEOM=published recovers the L = 90, eps_slab = 10
# geometry behind the submitted tables).
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
#   MULTICUBE_GEOM       substrate variant: "sec53" (default, = article Sec. 5.3) or
#                          "published" (L = 90, eps_slab = 10; the submitted tables)
#   MULTICUBE_EVAL       evaluation path: "both" (default; times the multi-RHS batched eval
#                          AND the original per-column loop on the same sigma, and checks they
#                          agree), "batched", or "percolumn"
#   CORRECT_EDGES=0      disable edge correction (default on, production)
#   MULTICUBE_GEOM_ONLY=1  print K per cutoff and exit (no solve)
#   MULTICUBE_GMRES_ITMAX  cap block-GMRES iterations (default 500 = converges). Set to 1
#                          for a TIMING-ONLY run: measures one iteration's cost so the full
#                          block-solve time can be reconstructed offline as
#                          t_1iter * niter_ref (niter from a real converged run), without
#                          paying for a full multi-hour low-thread-count solve. The residual
#                          and physics outputs (V, screen_ratio, ...) are NOT meaningful when
#                          itmax=1 — only t_precompute/t_pottrg/t_solve_block/t_eval are used.
#   MULTICUBE_RUNTAG      tag (e.g. "nt8") that routes output to data/raw_threads/ under a
#                          distinct filename, so a thread-count sweep can never collide with,
#                          or be swept up by the gather regex for, the production K-sweep
#                          records in data/raw/. Leave unset for normal (K-sweep) runs.
#
# Run (one cutoff): MULTICUBE_CUTOFF=1.5 JULIA_NUM_THREADS=96 OMP_NUM_THREADS=96 \
#                     julia --project=. exp66_multicube/scripts/run_multicube.jl
# Smoke (dev):      MULTICUBE_SMOKE=1 MULTICUBE_CUTOFF=1.5 julia -t 8 --project=. ... run_multicube.jl
# Thread-timing (one thread count): MULTICUBE_CUTOFF=0.0 MULTICUBE_GMRES_ITMAX=1 \
#     MULTICUBE_RUNTAG=nt8 JULIA_NUM_THREADS=8 OMP_NUM_THREADS=8 julia --project=. ... run_multicube.jl

using BoundaryIntegral
import BoundaryIntegral as BI
using Krylov
using LinearAlgebra
using Printf
using Serialization
using Random

const SMOKE = get(ENV, "MULTICUBE_SMOKE", "0") == "1"
const CORRECT_EDGES = get(ENV, "CORRECT_EDGES", "1") == "1"
# Right-diagonal preconditioner (BI design note 2026-07-31). This script calls
# Krylov.block_gmres directly, so it does NOT inherit BI's new
# `precondition = true` default on solve_dielectric_box3d_block.
const PRECONDITION = get(ENV, "PRECONDITION", "1") == "1"
const GEOM_ONLY = get(ENV, "MULTICUBE_GEOM_ONLY", "0") == "1"
# Thread-scaling knobs. GMRES_ITMAX=1 runs a single block-GMRES iteration (timing
# only; the solve does NOT converge and the physics outputs are meaningless): the
# per-iteration time is multiplied offline by niter(K=1 @ 96 threads) to recover
# the full block-solve time cheaply at low thread counts. RUNTAG tags the output
# filename so a thread sweep does not collide with the K-sweep records.
const GMRES_ITMAX = parse(Int, get(ENV, "MULTICUBE_GMRES_ITMAX", "500"))
# Evaluation path: "percolumn" (default), "batched" (multi-RHS), or "both" (run each on
# identical sigma and compare). See the eval block below for why percolumn is the default:
# this benchmark's targets are ALL in the near region, so the batched path has no far field
# to batch and is measurably slower.
# MULTICUBE_EVAL_ONLY=1 benchmarks the EVALUATION alone: the interface is still built (pottrg
# and the near corrections need it, and N sets the matvec cost), but the corrected LHS operator,
# the RHS assembly and block GMRES are all skipped and sigma is filled with a fixed-seed random
# field. Evaluation cost depends on sigma's SHAPE, not its values, so the timing is exact while
# the expensive solve is not paid for. Every physics output of such a record (V, screen_ratio,
# niter) is meaningless and is written as NaN.
const EVAL_ONLY = get(ENV, "MULTICUBE_EVAL_ONLY", "0") == "1"
const EVAL_MODE = get(ENV, "MULTICUBE_EVAL", "block")
EVAL_MODE in ("batched", "percolumn", "both", "block", "model") ||
    error("MULTICUBE_EVAL must be model|block|batched|percolumn|both, got $(EVAL_MODE)")
# "model" mode reads its FIXED target set (the union of all 198 Sec. 5.3 orbitals' quadrature
# points) from this cache, built by scripts/prep_model_targets.jl.
const TARGET_CACHE = get(ENV, "MULTICUBE_TARGET_CACHE",
                         joinpath(@__DIR__, "..", "data" * get(ENV, "RERUN_TAG", ""), "targets_model.jls"))
const RUNTAG = get(ENV, "MULTICUBE_RUNTAG", "")

# RERUN_TAG: appends a suffix to this experiment's output directory so a rerun
# never overwrites the data behind the submitted manuscript. Empty = original paths.
const TAG = get(ENV, "RERUN_TAG", "")
const DATA = joinpath(@__DIR__, "..", "data" * TAG, SMOKE ? "smoke" : "")
mkpath(joinpath(DATA, "raw")); mkpath(joinpath(DATA, "raw_threads"))
mkpath(joinpath(DATA, "raw"))

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1 = joinpath(REF_DIR, "graphene_00001.xsf")   # sublattice A (type 1)
const XSF_2 = joinpath(REF_DIR, "graphene_00002.xsf")   # sublattice B (type 2)

# graphene hexagonal lattice (bohr), demo_2x2 convention (as exp65)
const A1 = (2.465, 0.0, 0.0)
const A2 = (-1.2325, 2.1347526, 0.0)

# multicube substrate: two L x L x L cubes (eps1 | eps2) sharing x = c_x, slab on top.
#
# MULTICUBE_GEOM selects the substrate (default "sec53"):
#   "sec53"     -- the SAME system the article's lattice section (Sec. 5.3) reports:
#                  x3 converged Si|SiO2 cubes (L = 270) and the graphene slab at its cRPA
#                  eps = 2.4, slab lateral extent as in campaigns/lattice_conv_l3_eps2.4.toml.
#                  The eps = 10 slab of the published tables was a generic demo value, and
#                  L = 90 leaves a finite-size boundary artifact in the on-site U
#                  (see lattice_scale: 90 A gives a bowl, x3 is converged).
#   "published" -- L = 90, eps_slab = 10: the geometry behind the submitted Tables 2/3.
# Only the substrate boxes/epses change; the pipeline, discretization (n_quad, edge_level,
# tolerances) and the K sweep are untouched, so the timings stay comparable.
const GEOM = get(ENV, "MULTICUBE_GEOM", "sec53")
GEOM in ("sec53", "published") || error("MULTICUBE_GEOM must be \"sec53\" or \"published\", got $(GEOM)")
const SEC53 = GEOM == "sec53"
const L = SEC53 ? 270.0 : 90.0
# Slab thickness is a property of the graphene sheet's z-support, NOT of the cube edge: the
# published `L / 10` happened to equal 9 only because L was 90. Sec. 5.3 fixes it at 9 (the
# slab spans z in [3, 12] with the orbital plane at 7.5), so it is a constant here.
const SLAB_THICK = 9.0
# Slab lateral extent. Sec. 5.3 sizes it as (orbital-block extent + 2 x 11 A margin) so the
# phi^2 support is fully inside the slab; those exact dimensions are reused here. The largest
# cluster below (cutoff 5 A) has a 5 A position half-extent, so with the ~10 A in-plane phi^2
# support it sits well inside this slab's 21.5 A half-width. Fixed across K on purpose: an
# interface that grew with K would confound the multi-RHS scaling study.
const SLAB_LX, SLAB_LY = SEC53 ? (42.95357656, 44.06001031) : (L, L)
const EPS1, EPS2, EPS_OUT = 11.9, 3.9, 1.0
const EPS_SLAB = SEC53 ? 2.4 : 10.0

# MULTICUBE_ACCURACY selects the discretization and tolerances.
#
#   "sec53" (default) -- IDENTICAL to the lattice campaign (campaigns/lattice_conv_l3_eps2.4
#       .toml): edge_refine_level 3, rhs_tol = volume_tol = lhs_tol = gmres_rtol = 1e-5,
#       max_order 64, max_depth 128. The two sections then benchmark the same system at the
#       same accuracy, so their on-site U values are directly comparable.
#   "published" -- what the submitted tables used: edge_level 4 (a FINER l_ec, 0.568 vs 1.136
#       A) but rhs_tol = volume_tol = 1e-3 (a LOOSER source/interface refinement) and
#       max_depth 12. Neither uniformly tighter nor looser than the lattice campaign, which is
#       why the cross-validated on-site U differed by 0.17%.
#
# max_order is 64 in every set: p_up is clamped to it with no adaptive fallback for
# non-touching pairs and no warning, and on a real conv_l3 interface max_order = 8 left 92.5%
# of upsample pairs pinned at the cap.
const ACC = get(ENV, "MULTICUBE_ACCURACY", "sec53")
ACC in ("sec53", "published") ||
    error("MULTICUBE_ACCURACY must be sec53|published, got $(ACC)")

const P = SMOKE ?
    (n_quad = 6, edge_level = 2, rhs_tol = 1e-2, lhs_tol = 1e-3,
     gmres_atol = 1e-3, gmres_rtol = 1e-3, max_order = 64, max_depth = 12,
     support_rtol = 1e-3) :
    ACC == "sec53" ?
    (n_quad = 6, edge_level = 3, rhs_tol = 1e-5, lhs_tol = 1e-5,
     gmres_atol = 1e-5, gmres_rtol = 1e-5, max_order = 64, max_depth = 128,
     support_rtol = 1e-4) :
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
# nrange bounds the lattice-index search for neighbours. A fixed 3 silently TRUNCATES once the
# cutoff exceeds ~7.0 A: at cutoff 7.4 it finds 66 sites where the true count is 67, at 7.7
# it finds 70 of 73, at 8.0 it finds 74 of 79. The dropped sites are real neighbours, so both K
# and the physical cluster would be wrong with no error raised. Derive it from the cutoff
# instead (|A1| = |A2| = 2.465 A, plus margin for the non-orthogonal basis). Identical for every
# cutoff <= 7.0, where nrange = 3 was already sufficient -- so no existing point moves.
_nrange_for(cutoff) = max(3, ceil(Int, cutoff / 2.465) + 3)

function build_geometry(cutoff::Float64; nrange::Int = _nrange_for(cutoff))
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
        # slab over the junction (eps_slab). Sec. 5.3 anchors the same system on the junction
        # at absolute (5.547, 10.318); here it stays anchored on the central orbital, which is
        # the same configuration translated (junction directly beneath the orbital, cube tops
        # at the slab bottom z = cz - 4.5 = 3, slab z in [3, 12] as there).
        (center = (cx, cy, cz),      Lx = SLAB_LX, Ly = SLAB_LY, Lz = tz)]
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
    if !EVAL_ONLY
        t = @elapsed op = BI.batched_lhs_dielectric_box3d_fmm3d_corrected(
            interface, P.lhs_tol, P.lhs_tol, P.max_order; correct_edges = CORRECT_EDGES)
        note("batched LHS operator", t)
    else
        note("batched LHS operator", 0.0)     # skipped: EVAL_ONLY
    end

    sources = BI.batch_volume_sources(b)
    # Box-based multi-region screening: the interface spans 3 eps regions, so the
    # interface-based screened_volume_source exp65 used (uniform eps_in only) can't apply
    # here — screen each source by which box it sits in (all orbitals lie in the slab).
    # Only the RHS assembly and the per-column eval consume these; evaluate_batch_potential
    # screens internally. Under EVAL_ONLY both are skipped, so screening K sources here would
    # be pure overhead charged to a benchmark that is not measuring it.
    screened = EVAL_ONLY ? BI.VolumeSource{Float64, 3}[] :
        [BI.screened_volume_source(boxes, epses, EPS_OUT, sources[k], BI.SharpScreening()) for k in 1:K]

    # --- solve (scales with K) ---
    # RHS per source (the batched rhs_dielectric_box3d_fmm3d re-screens via the interface;
    # for >1 distinct-position source it loops the single-source path anyway, so we loop the
    # 4-arg single-source form with eps_src=1 on the already-screened sources).
    local F
    local sigma_block, bstats
    if EVAL_ONLY
        note("RHS assembly (per-source, multi-region)", 0.0)   # skipped: EVAL_ONLY
        # Fixed seed so the benchmark is reproducible; values are irrelevant to the timing.
        sigma_block = randn(Random.MersenneTwister(20260908), BI.num_points(interface), K)
        bstats = (; niter = 0)
        note("block GMRES", 0.0)                               # skipped: EVAL_ONLY
    else
        t = @elapsed begin
            F = Matrix{Float64}(undef, BI.num_points(interface), K)
            for k in 1:K
                F[:, k] = BI.rhs_dielectric_box3d_fmm3d(interface, screened[k], 1.0, P.rhs_tol)
            end
        end
        note("RHS assembly (per-source, multi-region)", t)
        t = @elapsed begin
            Nprec = PRECONDITION ?
                Diagonal(BI.dielectric_diagonal_scaling(interface)) :
                LinearAlgebra.I
            sigma_block, bstats = Krylov.block_gmres(op, F; N = Nprec,
                rtol = P.gmres_rtol, atol = P.gmres_atol, itmax = GMRES_ITMAX)
        end
        note("block GMRES", t)
    end
    # itmax=1 runs are timing-only (see GMRES_ITMAX doc above): the residual is not
    # meaningful, so skip the (also costly) matvec used only for the printed diagnostic.
    block_resid = (EVAL_ONLY || GMRES_ITMAX <= 1) ? NaN :
        norm(op * sigma_block - F) / max(norm(F), eps(Float64))

    # --- eval: central row V[rho_11, rho_b], target = onsite rho_11 = phi_1^2 (exp65 protocol).
    onsite = sources[1]
    nz = findall(>=(P.support_rtol * maximum(abs, onsite.density)), abs.(onsite.density))
    tgt = Matrix{Float64}(onsite.positions[:, nz])
    tw = onsite.weights[nz] .* onsite.density[nz]

    # Two evaluation paths, selected by MULTICUBE_EVAL:
    #
    #  "batched"   -- BI.evaluate_batch_potential: the MULTI-RHS eval. One corrected pottrg map
    #                 shared by all K columns, the Section-3 near/far split applied once, and
    #                 ALL far targets collapsed into a single nd=K point-charge FMM
    #                 (lattice_batch.jl). This is the path the distributed pipeline of the
    #                 article's lattice section runs.
    #
    #                 MEASURED: it does NOT help here, and cannot. Every pair density rho_1j
    #                 lies within ~6.5 A of the others with ~10 A support, so with c_pad = 5
    #                 every orbital's support is inside every other's near region: over the
    #                 whole K = 46 K x K block, 0 of 703 million (target, source) point pairs
    #                 are far-field, in all 2116 (target set, source) combinations. far_idx is
    #                 empty, the nd=K FMM never runs, and the path degenerates to K per-column
    #                 PrecomputedVolumeField builds -- 4% SLOWER than the loop below, from the
    #                 extra bookkeeping and from losing ltkm3dc's own kmax tuning. Batched eval
    #                 would only pay for well-separated pairs, which a localized Wannier basis
    #                 within one neighbour cutoff does not contain.
    #  "percolumn" -- the original loop: one ltkm3dc per column over the whole target set, no
    #                 near/far split, no batching. Retained because it produced the submitted
    #                 evaluation timings.
    #  "both"      -- run both on the SAME interface and sigma, time each, and report the
    #                 largest relative difference in Vrow. Comparing across separate jobs would
    #                 confound the eval comparison with run-to-run variation in the solve, so
    #                 the amortization claim for the eval stage is measured this way.
    #
    # t_pottrg/t_eval below always report the SELECTED path, so the tables keep their meaning;
    # t_eval_batched / t_eval_percolumn carry both when "both" ran. Note that in "both" mode the
    # wall-clock t_total necessarily includes the discarded pass.
    run_model     = EVAL_MODE == "model"
    run_block     = EVAL_MODE == "block"

    # ---- "model": K sources evaluated at a FIXED target set -- the union of the quadrature
    # points of ALL orbitals of the Sec. 5.3 model (N_p, independent of K).
    #
    # This is the benchmark that answers "what does one source cost to evaluate?". The K x N_p
    # potential matrix is what the ERI assembly needs: every V[rho_a, rho_b] for any target
    # pair a in the model is a contraction of a column of Phi against rho_a, so Phi IS the
    # evaluation work and t_eval/K is its per-source cost.
    #
    # With N_p FIXED, t_eval/K isolates the amortization: the corrected pottrg map over the
    # N_p targets is one K-independent cost shared by all K columns, so per-source time should
    # fall roughly as pottrg/K + (per-source field evaluation) and flatten to a floor. The
    # earlier "block" mode used the batch's own grid, which grows with K, so its flat per-source
    # time was two effects cancelling and measured nothing.
    local Phi_model
    t_eval_model = NaN
    n_tgt_model = 0
    if run_model
        isfile(TARGET_CACHE) || error("no target cache at $(TARGET_CACHE); run scripts/prep_model_targets.jl first")
        tc = deserialize(TARGET_CACHE)
        @printf("    model targets: N_p = %d from %d orbitals (%s)%s\n",
                size(tc.positions, 2), tc.n_orbitals, basename(tc.campaign),
                record ? "" : "  [warm-up: using the batch grid instead]")
        flush(stdout)
        # The warm-up pass exists only to force JIT; it compiles the identical method on the
        # batch's own (much smaller) grid. Paying the full N_p there would cost ~93 s at 96
        # threads and ~38 min at 1 thread, per task, measuring nothing.
        warm_tgt = record ? tc.positions : sources[1].positions
        t = @elapsed Phi_model = BI.evaluate_batch_potential(interface, sigma_block, sources,
                warm_tgt; lhs_tol = P.lhs_tol, volume_tol = VOLUME_TOL, c_pad = 5.0,
                screen_boxes = boxes, screen_epses = epses, screen_eps_out = EPS_OUT)
        n_tgt_model = size(warm_tgt, 2)
        t_eval_model = t
        note("eval model (K columns at fixed N_p targets)", t)
        @printf("    model eval: K = %d sources at N_p = %d targets in %.1f s  ->  %.2f s PER SOURCE  (%.3f us per K*N_p value)\n",
                K, n_tgt_model, t, t / K, 1e6 * t / (K * n_tgt_model))
        flush(stdout)
    end
    run_batched   = !EVAL_ONLY && EVAL_MODE in ("batched", "both")
    # model mode also runs the cheap central row, so v11_raw / screen_ratio stay comparable
    # with every other record in the sweep -- pointless under EVAL_ONLY, where sigma is random.
    run_percolumn = !EVAL_ONLY && EVAL_MODE in ("percolumn", "both", "model")

    # ---- "block": the FULL K x K tensor block on the union of all orbital quadrature points.
    #
    # This is the protocol the lattice section uses (lec_conv_single.jl passes the batch grid
    # as both source positions and targets), and it is where the multi-RHS evaluation actually
    # pays. Targets are the batch's shared grid -- the union of every orbital's quad points --
    # so ONE evaluate_batch_potential call gives Phi[g, b] for all K columns, and every
    # V[a,b] = sum_g w_g rho_a(g) Phi_b(g) then falls out as a single K x K matrix product.
    #
    # Cost: one pottrg build + K field evaluations + K matvecs for K^2 tensor entries, i.e.
    # O(K) work per O(K^2) entries -- the evaluation analogue of sharing one interface across
    # K right-hand sides. The per-column path below instead pays one field evaluation PER
    # ENTRY, because it evaluates only over rho_11's support and so cannot reuse a field.
    #
    # It is also more accurate: the row protocol truncates its target set at support_rtol,
    # which biases the contraction low. Here every grid point is a target, as in Sec. 5.3.
    local Vblock
    t_eval_block = NaN
    n_tgt_block = 0
    if run_block
        grid = sources[1].positions          # shared-positions contract: the union grid
        t = @elapsed begin
            Phi = BI.evaluate_batch_potential(interface, sigma_block, sources, grid;
                lhs_tol = P.lhs_tol, volume_tol = VOLUME_TOL, c_pad = 5.0,
                screen_boxes = boxes, screen_epses = epses, screen_eps_out = EPS_OUT)
            # target-side densities are UNSCREENED (evaluate_batch_potential screens the
            # sources internally), matching lattice_scale's onsite_U convention.
            Wrho = Matrix{Float64}(undef, size(Phi, 1), K)
            for a in 1:K
                Wrho[:, a] .= sources[a].weights .* sources[a].density
            end
            Vblock = Wrho' * Phi                       # V[a, b], the full K x K block
            n_tgt_block = size(Phi, 1)
        end
        t_eval_block = t
        note("eval block (K x K on union grid)", t)
        @printf("    block eval: %d x %d = %d tensor entries in %.1f s  (%.4f s/entry, %d targets)\n",
                K, K, K * K, t, t / K^2, n_tgt_block)
        # Exchange symmetry V[a,b] = V[b,a] is free once the whole block is in hand, and is the
        # same end-to-end check the lattice section quotes for the distributed assembly.
        if K > 1
            asym = maximum(abs.(Vblock .- transpose(Vblock))) / maximum(abs, Vblock)
            @printf("    block eval: max relative exchange asymmetry %.3e\n", asym)
        end
        flush(stdout)
    end

    local Vrow_b, Vrow_p, t_pot_p
    t_eval_b = t_eval_p = t_pot_p = NaN

    if run_batched
        t = @elapsed begin
            # NOTE: pass the UNSCREENED sources -- evaluate_batch_potential screens internally
            # (box-based, via screen_boxes) exactly as the RHS assembly above does.
            Phi = BI.evaluate_batch_potential(interface, sigma_block, sources, tgt;
                lhs_tol = P.lhs_tol, volume_tol = VOLUME_TOL, c_pad = 5.0,
                screen_boxes = boxes, screen_epses = epses, screen_eps_out = EPS_OUT)
            Vrow_b = [dot(tw, view(Phi, :, bcol)) for bcol in 1:K]
        end
        t_eval_b = t
        note("eval batched (shared pottrg + nd=K far FMM)", t)
    end

    if run_percolumn
        local pottrg
        t = @elapsed pottrg = BI.laplace3d_pottrg_fmm3d_corrected_hcubature(interface, tgt, P.lhs_tol, P.lhs_tol, 5.0)
        t_pot_p = t
        note("eval pottrg build (per-column path)", t)
        t = @elapsed begin
            Vrow_p = Vector{Float64}(undef, K)
            for bcol in 1:K
                sb = screened[bcol]                      # box-based multi-region screened source
                vals = BI.TKM3D.ltkm3dc(VOLUME_TOL, sb.positions; charges = sb.weights .* sb.density,
                                        targets = tgt, pgt = 1, kmax = BI._estimate_tkm3dc_kmax(sb))
                Vrow_p[bcol] = dot(tw, real.(vals.pottarg) .+ (pottrg * sigma_block[:, bcol]))
            end
        end
        t_eval_p = t
        note("eval per-column (K x ltkm3dc)", t)
    end

    # The two paths must agree: same operator, same sigma, different near/far bookkeeping.
    eval_rel_diff = NaN
    if run_batched && run_percolumn
        eval_rel_diff = maximum(abs.(Vrow_b .- Vrow_p) ./ max.(abs.(Vrow_p), eps(Float64)))
        @printf("    eval paths agree to %.3e (max rel. diff over K=%d entries)   speedup %.2fx\n",
                eval_rel_diff, K, (t_pot_p + t_eval_p) / t_eval_b)
        flush(stdout)
    end

    # under EVAL_ONLY none of the row paths ran, so there is no Vrow to select
    Vrow      = EVAL_ONLY ? Float64[] :
                run_block ? Vblock[1, :] : (run_percolumn ? Vrow_p : Vrow_b)
    t_pottrg  = run_percolumn ? t_pot_p : 0.0                # the others build pottrg internally
    t_eval_sel = run_model ? t_eval_model :
                 run_block ? t_eval_block :
                 (EVAL_MODE == "percolumn" ? t_eval_p : t_eval_b)

    # bare (interface-free) onsite self-energy reference: int rho_11 * TKM[rho_11 unscreened]
    v11_vac = NaN
    if !EVAL_ONLY
        q1 = onsite.weights .* onsite.density
        vac = BI.TKM3D.ltkm3dc(VOLUME_TOL, onsite.positions; charges = q1, targets = tgt,
                               pgt = 1, kmax = BI._estimate_tkm3dc_kmax(onsite))
        v11_vac = dot(tw, real.(vac.pottarg))
    end

    return (; K, n_points = BI.num_points(interface), n_src = size(b.densities, 1),
            niter = bstats.niter, block_resid,
            v11_raw = EVAL_ONLY ? NaN : Vrow[1], v11_vac,
            screen_ratio = EVAL_ONLY ? NaN : v11_vac / Vrow[1],
            Vrow = EVAL_ONLY ? Float64[] : Vrow, pairs, eval_only = EVAL_ONLY,
            # eval timings are returned rather than looked up by stage label, because which
            # stages exist now depends on EVAL_MODE
            eval_mode = EVAL_MODE, t_pottrg = t_pottrg, t_eval = t_eval_sel,
            t_eval_batched = t_eval_b, t_eval_percolumn = t_eval_p,
            t_eval_block = t_eval_block, t_eval_model = t_eval_model,
            n_tgt_model = n_tgt_model,
            n_entries = run_model ? K * n_tgt_model : (run_block ? K * K : K),
            Vblock = run_block ? Vblock : nothing,
            t_pottrg_percolumn = t_pot_p, eval_rel_diff)
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

@printf("threads = %d   correct_edges = %s   smoke = %s   eval = %s   accuracy = %s   cutoffs = %s\n",
        Threads.nthreads(), CORRECT_EDGES, SMOKE, EVAL_MODE, ACC, SWEEP)
@printf("accuracy set: n_quad %d  edge_level %d (l_ec %.4f)  rhs/volume_tol %.0e  lhs_tol %.0e  gmres_rtol %.0e  max_order %d  max_depth %d  support_rtol %.0e\n",
        P.n_quad, P.edge_level, SLAB_THICK / 2.0^P.edge_level * 1.01, P.rhs_tol, P.lhs_tol,
        P.gmres_rtol, P.max_order, P.max_depth, P.support_rtol)
@printf("geom = %s   substrate: two %gx%gx%g cubes eps %g|%g, slab %gx%gx%g eps %g, eps_out %g\n",
        GEOM, L, L, L, EPS1, EPS2, SLAB_LX, SLAB_LY, SLAB_THICK, EPS_SLAB, EPS_OUT)
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
    t_pottrg = res.t_pottrg
    t_eval = res.t_eval

    @printf("  K=%d  interface points %d  src %d  niter %d  block_resid %.2e\n",
            res.K, res.n_points, res.n_src, res.niter, res.block_resid)
    @printf("  V[11] = %.6e   V_vac[11] = %.6e   V_vac/V = %.4f\n",
            res.v11_raw, res.v11_vac, res.screen_ratio)
    @printf("  precompute %.1f s | pottrg %.1f s | block solve %.1f s | eval %.1f s | total %.1f s   peak RSS %.2f GB\n",
            t_precompute, t_pottrg, t_solve_block, t_eval, t_total, rss_gb())
    flush(stdout)

    out = (; smoke = SMOKE, correct_edges = CORRECT_EDGES, precondition = PRECONDITION, cutoff,
           accuracy = ACC, params = P, l_ec_used = SLAB_THICK / 2.0^P.edge_level * 1.01,
           geom = GEOM, L, slab_thick = SLAB_THICK, slab_lx = SLAB_LX, slab_ly = SLAB_LY,
           eps1 = EPS1, eps2 = EPS2, eps_slab = EPS_SLAB, eps_out = EPS_OUT,
           K = res.K, pairs = res.pairs, n_points = res.n_points, n_src = res.n_src,
           niter = res.niter, block_resid = res.block_resid,
           v11_raw = res.v11_raw, v11_vac = res.v11_vac, screen_ratio = res.screen_ratio, Vrow = res.Vrow,
           eval_only = res.eval_only, eval_mode = res.eval_mode, t_eval_batched = res.t_eval_batched,
           t_eval_percolumn = res.t_eval_percolumn, t_pottrg_percolumn = res.t_pottrg_percolumn,
           t_eval_block = res.t_eval_block, t_eval_model = res.t_eval_model,
           n_tgt_model = res.n_tgt_model, n_entries = res.n_entries, Vblock = res.Vblock,
           eval_rel_diff = res.eval_rel_diff,
           t_precompute, t_pottrg, t_solve_block, t_eval, t_total,
           stages = copy(stages), rss_baseline_gb = RSS_BASELINE, rss_peak_gb = rss_gb(),
           nthreads = Threads.nthreads(), gmres_itmax = GMRES_ITMAX, hostname = gethostname())
    # RUNTAG (e.g. a thread-count tag) routes output to its own raw_threads/ directory so a
    # timing-only itmax=1 thread sweep can never collide with, or be picked up by the gather
    # regex for, the production K-sweep records in raw/.
    if isempty(RUNTAG)
        serialize(joinpath(DATA, "raw", "multicube_K$(res.K)_cut$(cutoff).jls"), out)
    else
        mkpath(joinpath(DATA, "raw_threads"))
        serialize(joinpath(DATA, "raw_threads", "multicube_$(RUNTAG).jls"), out)
    end
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
        println(io, join(["hostname", "nthreads", "smoke", "correct_edges", "geom", "eps_slab", "cutoff",
            "K", "n_points", "n_src", "niter", "block_resid", "v11_raw", "v11_vac", "screen_ratio",
            "t_precompute", "t_pottrg", "t_solve_block", "t_eval", "t_total", "rss_peak_gb"], ","))
        for r in RESULTS
            println(io, join(string.([r.hostname, r.nthreads, r.smoke, r.correct_edges,
                r.geom, r.eps_slab, r.cutoff, r.K, r.n_points, r.n_src, r.niter, r.block_resid,
                r.v11_raw, r.v11_vac, r.screen_ratio,
                r.t_precompute, r.t_pottrg, r.t_solve_block, r.t_eval, r.t_total, r.rss_peak_gb]), ","))
        end
    end
    println("saved CSV -> ", joinpath(DATA, "multicube.csv"))
end
println("MULTICUBE MULTI-RHS DONE")

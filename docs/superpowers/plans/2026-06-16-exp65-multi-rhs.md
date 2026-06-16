# exp65 multi-RHS extension — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an instrumented multi-RHS benchmark driver to `exp65_orbital_bench` that solves a central graphene orbital plus its lattice neighbors (K grows with a neighbor cutoff) on one shared interface, measuring precompute amortization, block-vs-sequential GMRES, and peak RAM vs K.

**Architecture:** A single Julia script `scripts/run_multi_rhs.jl` mirroring the existing `run_single_rhs.jl` (warm-up, inline per-stage walltime + `Sys.maxrss()` snapshots, smoke mode, env knobs). It drives the revised `BoundaryIntegral.jl` lattice multi-RHS primitives directly — `assemble_lattice_batch` → explicit `batched_lhs_dielectric_box3d_fmm3d_corrected` + `Krylov.block_gmres` → `four_index_matrix` — plus a sequential single-RHS baseline on the same interface. An sbatch template runs the production sweep on one exclusive genoa node (user submits).

**Tech Stack:** Julia (juliaup), `BoundaryIntegral.jl` (live checkout), `Krylov.jl`, FINUFFT/FMM3D/TKM3D (vendored in the BI.jl env). Run via `/mnt/home/xgao1/.juliaup/bin/julia --project=/mnt/home/xgao1/codes/BoundaryIntegral.jl`.

**Spec:** `docs/superpowers/specs/2026-06-16-exp65-multi-rhs-design.md`

---

## Conventions / shared facts (read once)

- **Paths (absolute):**
  - Script dir: `/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/scripts`
  - BI.jl project: `/mnt/home/xgao1/codes/BoundaryIntegral.jl`
  - juliaup julia: `/mnt/home/xgao1/.juliaup/bin/julia`
  - Orbitals: `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15/graphene_0000{1,2}.xsf`
- **Smoke dev command** (use 8 threads in dev — 96 threads triggers the FINUFFT/FFTW pathology and is only for the production cluster run):
  ```bash
  ORBBENCH_SMOKE=1 /mnt/home/xgao1/.juliaup/bin/julia -t 8 \
    --project=/mnt/home/xgao1/codes/BoundaryIntegral.jl \
    /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/scripts/run_multi_rhs.jl
  ```
- **API facts (verified against the live checkout):**
  - `BI.read_xsf(path)` → tuple; `[2]` is the `DataGrid` (raw φ values).
  - `BI.density_centroid(dg)` → indexable length-3 (Cartesian centroid).
  - `BI.snap_orbital(dg, centroid::NTuple{3,Float64}, pos::NTuple{3,Float64})` → `NTuple{3,Int}` grid-step offset.
  - `BI.OrbitalInstance(id::Int, template_id::Int, steps::NTuple{3,Int})`.
  - `BI.assemble_lattice_batch(templates::Vector{DataGrid}, instances::Dict{Int,OrbitalInstance}, pairs::Vector{Tuple{Int,Int}}; support_rtol)` → `LatticeBatch` (all K columns share grid positions; templates must share grid geometry — asserted inside).
  - `BI.envelope_volume_source(b)` → `VolumeSource` (rss across columns).
  - `BI.batch_volume_sources(b)` → `Vector{VolumeSource{Float64,3}}` (length K, shared positions).
  - `BI.multi_dielectric_box3d_rhs_adaptive(n_quad::Int, l_ec, boxes::Vector{BoxGeom}, epses::Vector{Float64}, vs::VolumeSource, rhs_atol; eps_out, max_depth, tkm_kmax)` → `DielectricInterface`. `BoxGeom = (center=(cx,cy,cz), Lx=, Ly=, Lz=)`.
  - `BI.batched_lhs_dielectric_box3d_fmm3d_corrected(interface, fmm_tol, up_tol, max_order; correct_edges)` → `BatchedDielectricOperator`. `op * X::Matrix` works; `Krylov.block_gmres(op, F::Matrix; rtol, atol, itmax)` → `(Σ, stats)`.
  - `BI.rhs_dielectric_box3d_fmm3d(interface, vss::Vector{VolumeSource}, thresh)` → N×K (batched nd=K FMM; screens each source).
  - `BI.rhs_dielectric_box3d_fmm3d(interface, vs::VolumeSource, thresh)` → length-N (single-source; screens).
  - `BI.lhs_dielectric_box3d_fmm3d_corrected(interface, fmm_tol, up_tol, max_order; correct_edges)` → `LinearMap`; `Krylov.gmres(lhs, rhs; atol, rtol)` → `(σ, stats)`.
  - `BI.four_index_matrix(interface, sources::Vector{VolumeSource}, Σ; lhs_tol, volume_tol, range_factor)` → K×K `V[a,b]` (uses TKM3D for u_inc — hits the 96-thread pathology by design).
  - `BI._estimate_tkm3dc_kmax(vs::VolumeSource)` → Float64.
  - `to_eV(raw, Ni, Nj)` from `ScreenedOrbitalSolve` (same include as `run_single_rhs.jl`).
- **Graphene lattice (bohr):** `a1 = (2.465, 0, 0)`, `a2 = (-1.2325, 2.1347526, 0)` (from `lattice_scale/campaigns/demo_2x2.toml`). Sublattice A = `graphene_00001.xsf` (type 1), B = `graphene_00002.xsf` (type 2).
- **Box / params:** `L = 90`, `Lz = 3.35`, `eps_in = 3.5`, `eps_out = 1.0`; box center xy `(0,0)`, cz = the A centroid's z (≈ 7.5) → orbital at box **midplane** (demo_2x2 convention; differs from single-RHS top-face — see spec validation caveat).

---

## File Structure

- **Create** `exp65_orbital_bench/scripts/run_multi_rhs.jl` — the whole driver (geometry + pipeline + sequential baseline + eval/contract + instrumentation + sweep). Single cohesive script (~230 lines), matching the `run_single_rhs.jl` pattern.
- **Create** `exp65_orbital_bench/slurm/pvf_multi_rhs.sbatch` — exclusive-genoa sbatch (copy of `pvf_field.sbatch`, 96 threads, `--time=01:00:00`).
- **Writes at runtime:** `data/raw/multi_rhs_K<k>.jls` (per K), `data/multi_rhs.csv` (one row per K); smoke variants under `data/smoke/`.
- **Modify after the production run** `exp65_orbital_bench/REPORT.md` — addendum table + observations.

---

### Task 1: Driver skeleton + geometry helper (verify K progression)

**Files:**
- Create: `exp65_orbital_bench/scripts/run_multi_rhs.jl`

- [ ] **Step 1: Write the skeleton (header, imports, constants, geometry helper, geometry-only main)**

Create `run_multi_rhs.jl` with exactly this content:

```julia
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
```

- [ ] **Step 2: Verify the geometry helper produces the expected K progression**

Run:
```bash
ORBBENCH_GEOM_ONLY=1 /mnt/home/xgao1/.juliaup/bin/julia -t 4 \
  --project=/mnt/home/xgao1/codes/BoundaryIntegral.jl \
  /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/scripts/run_multi_rhs.jl
```
Expected: prints centroids near `(~0, ~0, ~7.5)` and a K progression `cutoff 0.00 -> K = 1`, `1.50 -> K = 4`, `2.50 -> K = 10`, `2.90 -> K = 13`, `3.90 -> K = 19` (exact K may differ by ±a couple if a shell sits near a cutoff — confirm it is monotonically increasing and K=1 at cutoff 0). `support points` should be > 0 (and, because far-shell pair products nearly vanish, may stay roughly flat across the larger cutoffs even as K grows). Ends with `GEOM ONLY DONE`.

If K at cutoff 0 is not 1, or the run errors in `assemble_lattice_batch` with a grid-incompatibility assertion, STOP — the two templates do not share a grid (revisit before continuing).

- [ ] **Step 3: Commit**

```bash
cd /mnt/home/xgao1/work/four_index_integral_solver
git add codes/numerical_results/exp65_orbital_bench/scripts/run_multi_rhs.jl
git commit -m "exp65 multi-RHS: driver skeleton + center+neighbors geometry helper"
```

---

### Task 2: Pipeline — precompute, block solve, sequential baseline, eval/contract

**Files:**
- Modify: `exp65_orbital_bench/scripts/run_multi_rhs.jl` (insert before the `if GEOM_ONLY` block)

- [ ] **Step 1: Add the instrumented per-K pipeline function**

Insert this `run_pipeline` function immediately after `build_geometry` and before the `if GEOM_ONLY` block:

```julia
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
    env = BI.envelope_volume_source(b)
    t = @elapsed kmax = BI._estimate_tkm3dc_kmax(env)
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
    note("sequential (1 op + K gmres)", t_seq_op + t_seq_solve)
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
```

- [ ] **Step 2: Verify one cutoff runs end-to-end in smoke mode**

Temporarily test a single cutoff by overriding the list. Run:
```bash
ORBBENCH_SMOKE=1 ORBBENCH_CUTOFFS=1.5 /mnt/home/xgao1/.juliaup/bin/julia -t 8 \
  --project=/mnt/home/xgao1/codes/BoundaryIntegral.jl -e '
  include("/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/scripts/run_multi_rhs.jl")' 2>&1 | tail -5 || true
```
This will error at the bottom because `main` is not written yet (the `if GEOM_ONLY` is false and nothing follows). That is expected. Instead verify `run_pipeline` directly:
```bash
ORBBENCH_SMOKE=1 /mnt/home/xgao1/.juliaup/bin/julia -t 8 \
  --project=/mnt/home/xgao1/codes/BoundaryIntegral.jl -e '
  ENV["ORBBENCH_GEOM_ONLY"]="0";
  include("/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/scripts/run_multi_rhs.jl");
  b, pairs = build_geometry(1.5);
  st = NamedTuple[];
  r = run_pipeline(b, pairs, st; record=true);
  println("K=", r.K, " niter=", r.niter, " block_resid=", r.block_resid,
          " seq_agree=", r.seq_agree, " max_rel_asym=", r.max_rel_asym,
          " u_onsite_ev=", r.u_onsite_ev)'
```
Expected: prints the per-stage table (envelope, interface build, batched LHS, RHS, block GMRES, sequential, eval) with positive times, then a summary line. Sanity gates:
- `block_resid` ≲ `gmres_rtol` (≲ 1e-2 in smoke),
- `seq_agree` small (≲ 1e-2 — block and sequential solve the same system),
- `max_rel_asym` small (≲ 1e-2 — V is ~symmetric),
- `u_onsite_ev` in the ~5–15 eV range (same order as the single-RHS ~10 eV; not expected to match exactly — midplane vs top-face convention).

If `seq_agree` or `max_rel_asym` are O(1), STOP and debug before proceeding (likely an operator/threshold mismatch between the block and sequential paths).

- [ ] **Step 3: Commit**

```bash
cd /mnt/home/xgao1/work/four_index_integral_solver
git add codes/numerical_results/exp65_orbital_bench/scripts/run_multi_rhs.jl
git commit -m "exp65 multi-RHS: precompute + block/sequential solve + four-index eval"
```

---

### Task 3: Warm-up, sweep loop, serialization, CSV, summary table

**Files:**
- Modify: `exp65_orbital_bench/scripts/run_multi_rhs.jl` (append after the `if GEOM_ONLY ... end` block)

- [ ] **Step 1: Append the main driver (warm-up + sweep + output)**

Append this to the end of the file (after the `if GEOM_ONLY` block):

```julia
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
    t_solve_seq = tof("sequential (1 op + K gmres)")
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
           t_precompute, t_solve_block, t_solve_seq, t_eval, t_total,
           stages = copy(stages), rss_baseline_gb = RSS_BASELINE, rss_peak_gb = rss_gb(),
           nthreads = Threads.nthreads(), hostname = gethostname())
    serialize(joinpath(DATA, "raw", "multi_rhs_K$(res.K).jls"), out)
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
            "max_rel_asym", "t_precompute", "t_solve_block", "t_solve_seq", "t_eval",
            "t_total", "rss_peak_gb", "u_onsite_ev"], ","))
        for r in RESULTS
            println(io, join(string.([r.hostname, r.nthreads, r.smoke, r.correct_edges,
                r.cutoff, r.K, r.n_points, r.n_src, r.niter, r.block_resid, r.seq_agree,
                r.max_rel_asym, r.t_precompute, r.t_solve_block, r.t_solve_seq, r.t_eval,
                r.t_total, r.rss_peak_gb, r.u_onsite_ev]), ","))
        end
    end
    println("saved CSV -> $csv")
end
println("MULTI RHS BENCH DONE")
```

- [ ] **Step 2: Verify the full smoke sweep**

Run the smoke command:
```bash
ORBBENCH_SMOKE=1 /mnt/home/xgao1/.juliaup/bin/julia -t 8 \
  --project=/mnt/home/xgao1/codes/BoundaryIntegral.jl \
  /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/scripts/run_multi_rhs.jl
```
Expected: warm-up line, then one block per cutoff (`0.0`, `1.5`) with the per-stage table and summary lines, then the summary table (rows for K=1 and K=4) and `MULTI RHS BENCH DONE`. Verify:
- `data/smoke/multi_rhs.csv` exists with a header + 2 rows,
- `data/smoke/raw/multi_rhs_K1.jls` and `multi_rhs_K4.jls` exist,
- across rows `seq_agree` and `max_rel_asym` are small (≲ 1e-2), and the speedup column (`spdup`) is reported (≥ ~1; the batched path should not be slower).

- [ ] **Step 3: Commit**

```bash
cd /mnt/home/xgao1/work/four_index_integral_solver
git add codes/numerical_results/exp65_orbital_bench/scripts/run_multi_rhs.jl \
        codes/numerical_results/exp65_orbital_bench/data/smoke/
git commit -m "exp65 multi-RHS: warm-up + cutoff sweep + CSV/jls output + summary table"
```

---

### Task 4: Production sbatch template

**Files:**
- Create: `exp65_orbital_bench/slurm/pvf_multi_rhs.sbatch`

- [ ] **Step 1: Write the sbatch (copy of `pvf_field.sbatch`, multi-RHS script, 1 h)**

Create `exp65_orbital_bench/slurm/pvf_multi_rhs.sbatch`:

```bash
#!/bin/bash
# exp65 multi-RHS benchmark — center orbital + lattice-neighbor pair densities,
# block solve vs sequential, full cutoff sweep (K up to ~19). One exclusive 96-core
# genoa node. Eval runs at 96 threads ON PURPOSE (documents the FINUFFT/FFTW thread
# pathology); four_index_matrix is K separate ~50 s TKM calls, so the largest-K eval is
# ~16 min and the full 5-cutoff sweep is ~60-90 min -> 2 h walltime.
# Submit yourself:  sbatch slurm/pvf_multi_rhs.sbatch
#SBATCH --job-name=pvf_multi_rhs
#SBATCH --partition=gen
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_multi_rhs_%j.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_multi_rhs_%j.err

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
WT=/mnt/home/xgao1/codes/BoundaryIntegral.jl
J=/mnt/home/xgao1/.juliaup/bin/julia          # juliaup binary (never module load julia)
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96

echo "=== node $(hostname)  cores $(nproc)  $(date) ==="
$J --project=$WT $NR/exp65_orbital_bench/scripts/run_multi_rhs.jl
echo "=== DONE $(date) ==="
```

- [ ] **Step 2: Verify the script is well-formed (syntax + headers)**

Run:
```bash
bash -n /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/slurm/pvf_multi_rhs.sbatch && echo "sbatch syntax OK"
grep -E "constraint=genoa|exclusive|time=02:00:00|run_multi_rhs.jl|JULIA_NUM_THREADS=96" \
  /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/slurm/pvf_multi_rhs.sbatch
```
Expected: `sbatch syntax OK` and all five grep lines present. (Do NOT submit — the user submits sbatch.)

- [ ] **Step 3: Commit**

```bash
cd /mnt/home/xgao1/work/four_index_integral_solver
git add codes/numerical_results/exp65_orbital_bench/slurm/pvf_multi_rhs.sbatch
git commit -m "exp65 multi-RHS: exclusive-genoa sbatch (96 threads, 1h, documents pathology)"
```

---

### Task 5: Production run (user-submitted) + REPORT addendum

This task runs AFTER tasks 1–4 are merged and smoke-validated. The production sweep is long (~30–40 min) and the org policy is that the **user submits `sbatch`** — the implementer does not.

**Files:**
- Modify: `exp65_orbital_bench/REPORT.md` (append an addendum)

- [ ] **Step 1: Ask the user to submit the production job**

Tell the user:
> Smoke-validated and committed. To get the production numbers, submit from `codes/lattice_scale`-style cwd:
> `cd codes/numerical_results/exp65_orbital_bench && sbatch slurm/pvf_multi_rhs.sbatch`
> I'll check status with `squeue -u $USER` / `sacct` and write the report addendum once `logs/slurm_multi_rhs_<jobid>.out` shows `MULTI RHS BENCH DONE`.

- [ ] **Step 2: After the job completes, read the results**

Read `exp65_orbital_bench/data/multi_rhs.csv` and the log `logs/slurm_multi_rhs_<jobid>.out`. Confirm `MULTI RHS BENCH DONE`, one CSV row per cutoff, `seq_agree`/`max_rel_asym` small, and the K=1 row's `u_onsite_ev` in the ~10 eV range.

- [ ] **Step 3: Append the REPORT addendum**

Append to `exp65_orbital_bench/REPORT.md` a section with this structure (fill the bracketed numbers from the CSV/log — these are real measured values, not placeholders, written after the run):

```markdown
# Addendum (2026-06-16): multi-RHS — central orbital + lattice neighbors

Script: `scripts/run_multi_rhs.jl` (commit <hash>), raw `data/raw/multi_rhs_K*.jls`,
log `logs/slurm_multi_rhs_<jobid>.out`. One exclusive genoa node, 96 threads.

Central graphene A-sublattice Wannier orbital + its neighbors (pairs
rho = phi_center * phi_neighbor) on ONE shared interface; K grows with the neighbor
cutoff. Box at the demo_2x2 midplane convention. Eval (four_index_matrix) runs at 96
threads — the TKM u_inc hits the documented FFTW thread pathology (THREAD_PATHOLOGY.md).

| K | interface pts | precompute (s) | block solve (s) | seq solve (s) | block speedup | eval (s) | total (s) | per-RHS (s) | peak GB |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | … | … | … | … | … | … | … | … | … |
| 4 | … | … | … | … | … | … | … | … | … |
| 10 | … | … | … | … | … | … | … | … | … |
| 13 | … | … | … | … | … | … | … | … | … |
| 19 | … | … | … | … | … | … | … | … | … |

Correctness: block vs sequential Sigma agree to <…> (rel); V symmetry
max|V - V'|/max|V| = <…>; K=1 onsite u = <…> eV (same order as the single-RHS
9.88 eV; differs by the midplane-vs-top-face convention — see the z-placement caveat).

Observations: <amortization — fixed precompute vs per-RHS marginal cost>;
<block vs sequential speedup and its source>; <peak RAM vs K>;
<eval/TKM 96-thread cost at this scale>.
```

- [ ] **Step 4: Commit**

```bash
cd /mnt/home/xgao1/work/four_index_integral_solver
git add codes/numerical_results/exp65_orbital_bench/REPORT.md \
        codes/numerical_results/exp65_orbital_bench/data/multi_rhs.csv \
        codes/numerical_results/exp65_orbital_bench/data/raw/multi_rhs_K*.jls
git commit -m "exp65 multi-RHS: production sweep results + REPORT addendum"
```

---

## Self-Review

**Spec coverage:**
- Real orbital pairs (center + lattice neighbors), cutoff sweep → Task 1 (`build_geometry`) ✓
- Shared interface + batched LHS (precompute amortization) → Task 2 precompute stages ✓
- Block-GMRES vs K sequential solves on the same interface → Task 2 (`run_pipeline`) ✓
- Four-index V contraction + symmetry → Task 2 (`four_index_matrix`, `max_rel_asym`) ✓
- Eval at 96 threads, document the pathology → Task 4 sbatch (96 threads) + Task 5 report ✓
- Phase split / per-RHS marginal / RAM vs K → Task 3 sweep + summary + CSV ✓
- K=1 single-RHS sanity anchor (with z-convention caveat) → Task 2/3 (`u_onsite_ev`), Task 5 ✓
- sbatch on exclusive genoa, juliaup, user submits → Task 4 ✓
- Smoke mode + env knobs → Task 1 constants, Task 3 verify ✓

**Placeholder scan:** No TBD/TODO in the script tasks. The REPORT addendum (Task 5) is filled with real measured numbers after the run — its bracketed `…` are explicitly post-run data, not unspecified plan content.

**Type/name consistency:** `build_geometry` returns `(b, pairs)`, consumed identically in Tasks 2/3; `run_pipeline(b, pairs, stages; record)` signature matches all call sites (warm-up, sweep, Step-2 verify); stage labels in `note(...)` exactly match the `tof("...")` lookups in Task 3 (`"envelope + tkm kmax"`, `"interface build (envelope)"`, `"batched LHS operator"`, `"RHS assembly (batched nd=K)"`, `"block GMRES"`, `"sequential (1 op + K gmres)"`, `"eval + contract (four_index_matrix)"`); `P` NamedTuple fields (`n_quad, edge_level, rhs_tol, lhs_tol, gmres_atol, gmres_rtol, max_order, max_depth, support_rtol`) are all referenced consistently.
</content>

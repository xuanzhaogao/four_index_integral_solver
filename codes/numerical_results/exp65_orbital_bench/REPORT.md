# Walltime and peak-RAM benchmark — production pipeline on the real monolayer system

Date: 2026-06-11. Hardware: worker7018 (96-core AMD Genoa, 1.5 TB), 96 Julia
threads + 96 OpenMP threads. Script: `scripts/run_single_rhs.jl` (commit
`15bf8fc`), raw record `data/raw/single_rhs_edges1.jls`, log
`logs/single_rhs.log`.

## What was measured

One full production four-index solve on the **real system**: graphene
monolayer modeled as a dielectric slab, with the **real Wannier pz orbital**
(`k_323201_nb_144_c_15/graphene_00001.xsf`, 57 MB) as the source density.
Single RHS (source = |phi_1|^2), single evaluation (onsite channel: target
density = source density).

The benchmark inlines `ScreenedOrbitalSolve.solve_screened_mode` stage by
stage — identical BI calls, identical calibrated parameters as
`screened_monolayer_hund_calibrated.jl`:

| parameter | value |
|---|---|
| slab | L x L x Lz = 90 x 90 x 3.35, eps_in = 3.5 (calibrated), vacuum outside |
| orbital placement | centroid of \|phi_1\|^2 at (0, 0, Lz/2), production convention (see caveat below) |
| source_tol | 1e-3 (-> 355,862 truncated grid points, N_phi1 = 78.416205) |
| panel order / edges | n_quad = 6, edge_refine_level = 4, **correct_edges = true** |
| tolerances | rhs_tol 1e-3, lhs_tol 1e-5, gmres atol = rtol = 1e-5 |
| caps | max_order 64, max_depth 12 |

A tiny Gaussian warm-up solve runs through every stage first, so all reported
times are warm (JIT-free). RAM is tracked two ways, and they agree:
`Sys.maxrss()` snapshots after each stage (process high-water mark — the
per-stage *increment* attributes new peak memory to that stage), and an
OS-level `/usr/bin/time -v` wrapper (18,027,040 kB = 17.19 GiB).

## Headline result

Problem size: **960,768 interface points** (reproduces the 960,552 of the
production records), 355,862 source = target points. GMRES: 8 iterations,
relative residual 7.1e-6. Result: u_onsite = 9.8849 eV at eps_in = 3.5
(production value at eps_in = 2.0 is 12.77 eV — consistent ordering).

| phase | wall | share | breakdown |
|---|---:|---:|---|
| load | 3.2 s | 1% | XSF read + tol-truncation (one-time, I/O) |
| **precompute** | **218.7 s** | **74%** | screened source 0.0 + TKM kmax 0.4 + **interface build 201.4** + LHS near-correction assembly 16.9 |
| **solve** | **25.3 s** | **9%** | RHS assembly (FMM) 2.2 + GMRES 23.1 |
| **eval** | **48.6 s** | **16%** | volume potential (TKM) 44.2 + scattered potential (FMM + hcubature) 4.5 |
| total | 295.8 s (~5 min) | | end-to-end wall 297.8 s; outer wall incl. startup + warm-up 5:25 |

**Peak RAM: 17.2 GB.** Build-up: 1.9 GB baseline after warm-up → 2.0 GB after
orbital load → **11.0 GB during the interface build** → flat through
LHS/RHS assembly → **17.2 GB during GMRES** (FMM workspaces over the 0.96M-point
operator) → flat through evaluation.

Per-stage record (maxrss = high-water mark at end of stage):

| stage | wall (s) | maxrss (GB) | live heap (GB) |
|---|---:|---:|---:|
| load (XSF read + truncation) | 3.18 | 2.02 | 0.47 |
| screened_volume_source | 0.00 | 2.02 | 0.51 |
| TKM kmax estimate | 0.37 | 2.36 | 0.61 |
| interface build (RHS-adaptive) | 201.40 | 10.99 | 0.28 |
| LHS assembly (near corrections) | 16.90 | 10.99 | 3.40 |
| RHS assembly (FMM) | 2.19 | 10.99 | 3.72 |
| GMRES solve (8 iters) | 23.11 | 17.19 | 2.28 |
| eval: volume potential (TKM) | 44.18 | 17.19 | 3.01 |
| eval: scattered potential (FMM+hcub) | 4.46 | 17.19 | 0.86 |

## Observations

1. **Precompute dominates (74%) and is almost entirely the RHS-adaptive
   interface build** (201 s of 219 s). This is exactly the RHS-independent
   part: with the mesh and operator reused, the marginal cost of each
   additional RHS is solve + eval ≈ 74 s — a ~4x amortization already at the
   second right-hand side. This is the quantitative anchor for the multi-RHS
   amortization story.
2. **Evaluation cost is the bare interaction, not the screening.** The TKM
   volume potential (u_int, present even in vacuum) costs 44 s; the BIE
   scattered-potential evaluation adds only 4.5 s. The dielectric correction
   is nearly free at evaluation time.
3. **The near-correction assembly is subdominant** — 16.9 s (~6%) even with
   edge correction enabled at this extreme slab aspect ratio (90:90:3.35).
   (6.1 showed edge correction is load-bearing for accuracy; here it is also
   cheap.)
4. **Parallel efficiency is the optimization target.** Average CPU
   utilization was 2326% of 9600% (~24 of 96 cores) — the interface build's
   panelwise adaptive recursion is the serial bottleneck. Halving precompute
   would nearly halve total wall.
5. **Memory is modest**: 17.2 GB peak for a ~1M-DOF production solve, far
   below the worker's 1.5 TB; the two RSS measures (in-process and OS) agree
   to 3 digits.

## Caveat: orbital z-placement in the production geometry

All monolayer production scripts shift the orbital centroid to
`z_center = Lz/2` and the docstring describes this as "slab midplane" — but
both `screened_volume_source(Lx, Ly, Lz, ...)` and
`single_dielectric_box3d_rhs_adaptive` build the slab **centered at the
origin** (z in [-Lz/2, +Lz/2]). The carbon plane therefore sits on the slab's
**top face**, with roughly half the pz density outside the dielectric. The
eps_in = 3.5 calibration was performed in this same geometry, so production
results are internally consistent; but if midplane was the intent,
`z_center` should be 0 and the calibration would shift. This benchmark
mirrors production as-is — the cost numbers are unaffected by the choice.

## Validation

- Smoke configuration (coarse tolerances, same geometry + real orbital):
  216,000 interface points, 11.4 s total, 3.5 GB peak, u_onsite = 10.37 eV —
  consistent with the production-tolerance result.
- N_phi1 = 78.416205 and the interface size reproduce the production CSV
  records exactly / to 0.02%.
- GMRES residual 7.1e-6 at the 1e-5 target; u_onsite ordering vs the
  eps_in = 2.0 production value is physically correct.

## Follow-ups (not yet run)

- `CORRECT_EDGES=0` single row for the legacy-config comparison.
- Multi-RHS sweep reusing interface + LHS: total time vs number of RHS
  (the phase split here already provides the per-RHS marginal cost).

---

# Addendum (2026-06-11): mesh-generation deep dive and the precompute-once fix

Script: `scripts/bench_meshgen.jl` (raw record `data/raw/bench_meshgen.jls`,
log `logs/bench_meshgen.log`). Three parts: (A) instrumented replication of
the production adaptive build, (B) component microbenchmarks, (C) prototype
of the proposed strategy — fix the near-evaluation region from the density,
precompute type-1 NUFFT + diagonal kernel scaling once, type-2 per round.

## Where the 201 s actually goes (part A)

Part A reproduces the production build exactly: 200.7 s vs 201.4 s, same
final mesh (26,688 panels = 960,768 points), 9 refinement depths.

| component | total | notes |
|---|---:|---|
| **ltkm3dc near calls** | **188 s (94%)** | 8 calls x ~23.5 s, **flat whether 14 or 4,352 targets are evaluated** |
| lfmm3d far calls | 12.5 s | also per-call setup dominated (B1: ~1.4-1.7 s/call nearly independent of target count) |
| KDTree classify | 0.5 s | negligible (build over 356k points: 0.03 s) |
| interpolation checks | ~0 s | |

The hypothesis is confirmed: the cost is the repeated per-call setup of the
volume-field evaluation over the full 356k-point density — dominated by the
TKM path (type-1 NUFFT + kernel scaling + gradient coefficients + FFT
planning, redone from scratch at every depth because the Fourier box depends
on the per-call target set).

Decomposed single-call stages on a FIXED box (B4, kmax = 40.21, modes
(441, 401, 353) = 6.2e7, coeff 0.93 GB): type-1 3.5 s + scale 0.4 s + grad
coeffs 0.3 s + type-2 (ntrans=3, plan+exec) ~1.0 s ≈ 5.2 s — vs 23.5 s per
ltkm3dc call. The ~18 s/call gap is per-call overhead inside ltkm3dc;
the likely mechanism (inference, consistent with all measurements): the
target-dependent box gives different FFT grid dimensions every call, so
FINUFFT/FFTW re-plans and re-allocates the ~5e8-point upsampled grid each
time (including one plan that is created and never executed on the
gradient-only path), whereas the fixed-box prototype plans the same
dimensions every round (FFTW wisdom reuse) and pays ~1 s. Worth a profile
when upstreaming.

Two side findings:
- `estimate_kcut3dc(1e-4)` finds **no usable spectral cutoff** (kcut ≈
  Nyquist): the screened density is discontinuous at the slab faces (the
  orbital straddles the top face in this geometry) and has nuclear cusps, so
  its spectrum has no tail below the grid Nyquist. The production
  kmax = 40.21 (grid-spacing Nyquist) is the right operating point — and
  part C validates it to 1.4e-5.
- The threaded direct sum runs at 7.3e10 pairs/s (B3): 356k sources x 10k
  targets in 0.05 s. At these batch sizes the far path needs no FMM at all.

## The proposed strategy works (part C)

Fixed near box B = density bbox + 5h margin; type-1 + truncated-kernel
scaling + spectral gradient coefficients precomputed ONCE on B's
target-independent Fourier box (4.1 s, 6.2e7 modes); per depth: type-2 NUFFT
at the targets inside B (~0.9 s/depth) + direct threaded sum outside
(~0.02 s/depth). Replaying the identical refinement path:

|  | current (part A) | precompute-once (part C) |
|---|---:|---:|
| adaptive-build evaluation | 200.7 s | **12.4 s** (4.1 precompute + 8.3 eval) |
| max RHS deviation | — | 8.0e-5 abs = **1.4e-5 rel** (budget: rhs_atol 1e-3) |
| refinement-decision flips | — | **0** (identical mesh) |

**16x on the dominant phase.** Projected single-RHS totals if upstreamed:
precompute 219 s -> ~30 s (interface 201 -> 12; LHS assembly 17 s unchanged),
end-to-end 296 s -> ~107 s. The same precomputed coefficients can also serve
the evaluation phase's volume potential (currently a separate 44 s ltkm3dc
call at 356k targets -> one type-2 exec), taking the projected total to
~70 s — a ~4x end-to-end speedup, and the per-RHS marginal cost drops
accordingly for the multi-RHS story.

## Upstream implementation: measured (2026-06-11, branch feature/precomputed-volume-field)

Implemented in a BI.jl worktree (`~/codes/BoundaryIntegral.jl-wt-pvf`, 11
commits, insertion-only, 33 new tests; plan
`docs/superpowers/plans/2026-06-11-precomputed-volume-field.md`): exported
`PrecomputedVolumeField`, `volume_field_potential`, `volume_field_gradient`,
`rhs_dielectric_box3d_field`, plus a field overload of
`single_dielectric_box3d_rhs_adaptive`. Production-scale benchmark of the
REAL branch code (`scripts/bench_field_path.jl`, worker7018, 96 threads,
same-session A/B, raw record `data/raw/bench_field_path.jls`):

|  | production path | field path | speedup |
|---|---:|---:|---:|
| field construction (once) | — | 4.6 s | |
| adaptive mesh build | 200.2 s | **9.0 s** | **22.2x** (14.7x incl. construction) |
| RHS assembly | 3.3 s (FMM) | 5.0 s | 0.7x (see note) |
| u_int volume potential (356k targets) | 43.5 s (ltkm3dc) | **0.39 s** | **111x** |

- Mesh equality: exact — 960,768 points on both paths.
- u_int agreement: 4.0e-5 relative. Peak RSS 15.0 GB (incl. both paths).
- RHS note: the field assembly is slightly slower at 960k interface targets
  (most fall outside the near box -> large direct sum) and deviates 4.7e-3
  from the FMM reference — but that deviation is dominated by the PRODUCTION
  path's own error: `Rhs_dielectric_box3d_fmm3d` treats the density as point
  charges even for panels touching the support, while the field's spectral
  near evaluation resolves the continuous density (the same effect measured
  in the unit tests, where the point-sum reference is 2-5% off near the
  support). The field RHS is the more accurate of the two.
- Projected end-to-end single RHS with the field path: load 3.1 + field 4.6
  + build 9.0 + LHS 16.9 + RHS 3.3 + GMRES 23.1 + u_int 0.4 + u_scatter 4.5
  ≈ **65 s vs 296 s (~4.5x)**.

Status: branch ready to merge (final review passed; live checkout untouched);
merging is the user's call.

## Evaluation-routing policy (decision note, 2026-06-11)

Per user decision (commit `cc6a749`), `PrecomputedVolumeField` evaluation
routes targets as follows, with no size-based switching:

- **Inside the field box (near the density): always the TKM spectral path**
  (type-2 NUFFT on the stored, truncated-kernel-scaled coefficients). This is
  a correctness requirement, not a performance choice: near the support, any
  point-charge representation of the density (direct summation or FMM) carries
  the quadrature error of the source grid — measured 2-5% on the test
  Gaussian, and the dominant part of the 4.7e-3 RHS deviation observed against
  `Rhs_dielectric_box3d_fmm3d` at production scale. Only the spectral
  evaluation resolves the continuous density there.
- **Outside the box: always FMM** (`lfmm3d` at the field tolerance). For
  well-separated targets the exact point sum and FMM agree within tolerance,
  so correctness is equivalent; FMM is asymptotically scalable in the batch
  size (the direct sum used previously is faster below ~3e5 targets/batch
  because each lfmm3d call rebuilds the source tree, ~1.4 s at 356k sources,
  but becomes quadratic beyond). Cost consequence at production scale: the
  adaptive build pays the FMM setup floor on the depths that have far targets
  (~1.4 s x ~6 depths), and the 960k-target RHS assembly gets faster
  (FMM 3.3 s vs direct 5.0 s).

The deleted direct-sum kernels remain available in git history if a
small-batch fast path is ever wanted again.

## Profiling + cache_fft optimization (2026-06-11, commits 3b0a519 / 810d736)

Step-by-step profiling (`scripts/profile_field_path.jl`) localized the
remaining cost: every FINUFFT type-2 exec re-does the FFT of the FIXED
coefficient grid (0.4 s pot / 0.85 s grad per call, independent of target
count — 2k and 356k targets identical), and FINUFFT auto-selects the fast
small-upsampling FFT at tol 1e-4 already (forcing upsampfac 2.0 is ~20x
WORSE: 73 s vs 3.7 s type-1), so the FFT-per-exec was the only lever left.

**cache_fft mode** (optional, default off): at construction, deconvolve the
coefficients by the ES-kernel Fourier factors, bfft once onto a padded
(x1.25) fine grid, store the grids (coefficient arrays dropped); per
evaluation, in-box targets are interpolated natively in threaded Julia using
TKM3D's spread-only kernel machinery. (FINUFFT's own spreadinterponly exec
was validated as numerically correct but rejected: finufft 2.5.x
value-initializes its nf-sized internal workspace on every exec, making it
grid-bound anyway. The native interpolation matches FINUFFT's
spreadinterponly exec to 1.4e-8.) Review-hardened: interpolation kernel
parameters frozen in the struct at construction; the TKM3D private-API
surface (8 underscore functions + TKM3D.FFTW) is documented in-code and
should be promoted to a public TKM3D spread-interp API before wide use.

Production-scale results (`scripts/bench_cache_fft.jl`, idle worker, 96
threads; mesh identical at 960,768 points; u_int agreement 5.0e-7):

| stage | standard field | cache_fft field | vs ORIGINAL production |
|---|---:|---:|---:|
| field construction (once) | 3.8 s | 7.3 s | — |
| adaptive mesh build | (9-24 s) | **2.4 s** | 200.2 s -> **85x** |
| RHS assembly | 3.2 s | 0.85 s | |
| u_int volume potential | 0.42 s | **0.03 s** | 43.5 s -> **1450x** |

Memory: cached grids ~2x the coefficient arrays (~8 GB at production);
construction transient ~16 GB (documented in the docstring).

Updated end-to-end single-RHS projection with cache_fft: load 3.1 + field
7.3 + build 2.4 + LHS 16.9 + RHS 0.9 + GMRES 23.1 + u_int 0.03 + u_scatter
4.5 ≈ **58 s vs 296 s (~5x)** — and the RHS-independent precompute drops
from 219 s to ~10 s, so the multi-RHS marginal cost is now GMRES-dominated.

Earlier profile-run caveats: the 23.7 s build / erratic FMM-floor numbers in
the profiling run were contaminated by ~10 cores of external load on the
worker; the 18 s "scale" reading was a script artifact (untyped global
closure), the real constructor scaling loop is 0.4 s.

## Slurm-validated numbers + warm ltkm3dc decomposition (2026-06-15)

Three jobs on dedicated **exclusive** genoa nodes (`slurm/pvf_{field,cache,
decomp}.sbatch`), all warm. These supersede the interactive worker numbers
(which drifted with node load):

**Field A/B (worker7020, exclusive):** field construction 3.9 s; adaptive
build 20.7 s (field) vs 195.5 s (VolumeSource) = **9.5x**, mesh identical
(960768); RHS 2.5 s (field) vs 0.7 s (FMM); **u_int 0.40 s (field) vs 42.5 s
(old ltkm3dc) = 106x**; u_int agreement 4.0e-5; peak 14 GB.

**cache_fft (worker7093, exclusive):** standard u_int 0.42 s -> **cached
u_int 0.10 s**; cached adaptive build **4.1 s**; cached vs standard agreement
5.0e-7; peak 27 GB.

**Decomposition of the warm old ltkm3dc u_int (worker7095, exclusive):**
full call **42.6 s**, and it splits cleanly as **type-1 NUFFT ~21.5 s +
type-2 NUFFT ~21.3 s** (forward build of the spectrum + inverse read-out to
targets), scaling negligible. So the old per-call cost is exactly TWO full
NUFFT passes over the ~5e7-mode grid. The field does the type-1 ONCE at
construction and replaces the per-call type-2 with a 0.40 s stored-coefficient
read-out.

Two honesty notes on the decomposition script (`diag_ltkm_decomp.jl`):
1. Its standalone "scale loop" (51 s) is a TYPE-UNSTABLE top-level-global
   artifact, not a real cost — inside the typed `_ltkm3dc_eval` the scaling is
   ~1-2 s, which is why type-1 + type-2 alone (~42.8 s) already equals the full
   call (42.6 s). The "sum of decomposed (94 s) > full (42.6 s)" line is
   entirely this artifact; it is NOT contention (it reproduced on a clean
   exclusive node).
2. **Unresolved, flagged honestly:** the per-NUFFT cost is strongly
   tolerance-dependent in a counterintuitive direction. At the production
   u_int tol (volume_tol = rhs_tol = **1e-3**) each NUFFT is ~21 s; at the
   field's tol (**1e-4**) the type-1 is 3.9 s and the type-2 is 0.4 s — i.e.
   the TIGHTER tolerance is ~5-50x FASTER here, at the same auto upsampfac
   (1.25; forced 2.0 is 64 s). This means part of the field's u_int advantage
   is that it operates at 1e-4 while the legacy `evaluate_volume_potential`
   path runs ltkm3dc at 1e-3. Actionable corollary: calling the OLD u_int at
   tol 1e-4 should itself be several-fold faster (~8 s vs 42 s) — worth a
   one-line change in the legacy path independent of the field work. The root
   cause of the tol-1e-3-slower-than-1e-4 behavior (a FINUFFT spread/FFT
   heuristic quirk near 1e-3) is not yet isolated.

Net, grid-independent and node-independent conclusion: the old runs were slow
because every call rebuilt the source spectrum (type-1) and read it out
(type-2) — two full NUFFT passes on a grid that is the SAME size as (slightly
smaller than) the field's. The field pays the forward pass once; cache_fft
additionally pre-FFTs so even the read-out is interpolation-only (0.10 s).

## ROOT CAUSE of the slow ltkm3dc: a FFTW 96-thread pathology (2026-06-15)

The "why is each old call ~50 s" investigation finally resolved by
instrumenting `_ltkm3dc_eval` internally and sweeping the FINUFFT thread count
on one idle node (`diag_ltkm_internal.jl`, worker7001). The decomposition is
now internally consistent (sum ~= the real full call):

```
REAL ltkm3dc(pgt=1) full call, 96 threads:  54.3 s
box        nthreads   type-1    scale    type-2    sum
field        96        3.61      0.44     0.43      4.47
field         1        2.77      0.44     2.26      5.48
ltkm3dc      96       27.51      0.35    27.23     55.09   <- the 54 s
ltkm3dc      16        0.48      0.36     0.42      1.26    <- 43x faster
ltkm3dc       1        1.96      0.35     1.90      4.22
```

**It is a FFTW multithreading pathology, not tol, grid size, or even grid
dimensions per se.** On the ltkm3dc per-call box's mode grid (413,375,325) the
FFT at 96 threads costs 27 s per transform but only 0.48 s at 16 threads — a
~57x slowdown from using MORE threads (FFTW picks a pathological plan for
those upsampled dimensions at high thread count). The field box (441,401,353)
does not trigger it (4.5 s at any thread count). The full 54 s call is exactly
type-1 + type-2, both hitting the pathology.

**Corrections to earlier claims in this report (the headline factors were
inflated by this bug):**
- The earlier "tol effect" paragraph is WRONG — `bench_tol_isolation.jl`
  showed tol 1e-3 ~= 1e-4 at fixed grid; tol is irrelevant. The real lever is
  the FFTW thread count.
- The u_int "1450x / 137x / 111x" figures were dominated by this 96-thread
  pathology in the OLD path. A thread-capped old ltkm3dc u_int is ~1.3 s, not
  ~50 s. The field's TRUE u_int advantage over a correct old path is the reuse
  factor (type-1 done once) ~= 3x, plus cache_fft's interp-only read-out.
- The 9.5-22x adaptive-build speedups are likely ALSO partly this bug (the
  VolumeSource build's per-depth ltkm3dc near-evals can hit pathological
  per-depth boxes) and should be re-measured with a thread cap before being
  quoted. The field/cache_fft absolute build times (21 s / 4 s) stand; the
  *ratio* to the old path needs the thread-capped baseline.

**Actionable fix (independent of PVF; do in a TKM3D worktree per the workflow
rule):** cap the FINUFFT `nthreads` (~16) inside `ltkm3dc`/`_ltkm3dc_eval`
(or choose mode-grid dimensions that thread well). This turns the 54 s call
into ~1.3 s (43x) and speeds up the ENTIRE existing pipeline, not just the
field path. The field path remains valuable for its structural reuse and
because its box dimensions avoid the pathology by construction — but the
dramatic speedup numbers must be re-baselined against a thread-fixed ltkm3dc.

---

# Addendum (2026-06-16): multi-RHS — central orbital + lattice neighbors

Script: `scripts/run_multi_rhs.jl` (commit `548343e`), sbatch
`slurm/pvf_multi_rhs.sbatch`, raw `data/raw/multi_rhs_K*_cut*.jls`, CSV
`data/multi_rhs.csv`, log `logs/slurm_multi_rhs_6519476.out`. One **exclusive
genoa node** (worker7001, 96 Julia + 96 OpenMP threads); job 6519476,
warm + 01:07:47 wall. Same calibrated production parameters as the single-RHS
benchmark (`n_quad = 6`, edge_refine_level 4, rhs_tol 1e-3, lhs_tol 1e-5,
gmres 1e-5, max_order 64, **correct_edges = true**).

## What was measured

The multi-RHS path of the revised BI.jl: a central graphene A-sublattice
Wannier orbital plus its lattice neighbors form **K pair densities**
`rho = phi_center * phi_neighbor` (onsite included) on **one shared interface**;
K grows with a neighbor cutoff (graphene shells). Each K is run through the
production pipeline stage by stage — `assemble_lattice_batch` →
`multi_dielectric_box3d_rhs_adaptive` (envelope) →
`batched_lhs_dielectric_box3d_fmm3d_corrected` → `Krylov.block_gmres` →
`four_index_matrix` (the K×K Coulomb matrix V) — alongside a **sequential
baseline** (the single-RHS operator built once on the same interface, then K
independent `Krylov.gmres` solves). Geometry follows the
`lattice_scale/demo_2x2` convention (box center cz = 7.5, orbital at the slab
**midplane**), which differs from `run_single_rhs.jl` (top face).

The interface is genuinely production scale — **~960,768 points at every
cutoff** (reproduces the single-RHS 960,768), independent of K; only the source
support (n_src 136k → 438k) and K (1 → 19) grow.

## Headline numbers

`block`/`seq` are the **solve** phases (batched RHS + block-GMRES vs K single
RHS + K gmres), both excluding the one-time LHS operator build (which is
counted in `precompute`); `spdup = seq/block`. `prod` is the **production
end-to-end** block path = precompute + block + eval (the sequential baseline is
a benchmark-only extra and is excluded). Times in seconds.

| cutoff | K | precompute | block | seq | spdup | eval | prod total | per-RHS marginal | peak RAM |
|---:|--:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 1 | 38.9 | 17.5 | 19.6 | 1.12x | 64.5 | 120.9 | 82.0 | 18.5 GB |
| 1.5 | 4 | 35.6 | 31.4 | 60.7 | 1.93x | 94.6 | 161.6 | 31.5 | 40.6 GB |
| 2.5 | 10 | 149.2 | 62.3 | 163.0 | 2.61x | 456.8 | 668.3 | 51.9 | 89.6 GB |
| 2.9 | 13 | 145.3 | 75.5 | 233.9 | 3.10x | 599.1 | 819.9 | 51.9 | 114.8 GB |
| 3.9 | 19 | 162.2 | 109.9 | 325.7 | 2.96x | 888.2 | 1160.3 | 52.5 | 164.0 GB |

`per-RHS marginal = (block + eval)/K` (cost of one more RHS on a solved batch).

## Observations

1. **Block solve beats K sequential solves, and the win grows with K** —
   1.12x (K=1, trivially ≈1) → 1.93x → 2.61x → **3.10x** (K=13) → 2.96x (K=19).
   The batched `nd=K` FMM matvec inside block-GMRES (one tree over the shared
   interface points serving all K columns) plus the batched RHS assembly
   amortize the FMM setup that the sequential path pays K times. This is the
   quantitative payoff of the multi-RHS API: the four-index V row costs ~3x less
   to solve as a block than column by column.

2. **Precompute is RHS-independent and amortizes across K.** The interface
   build + LHS operator (~39 s at K=1 up to ~162 s at K=19 — it grows with the
   source support, not with K per se) is built **once** per batch and serves all
   K right-hand sides. Run as K independent single-RHS jobs you would pay it K
   times; here the K=19 precompute (162 s) is amortized to ~8.5 s/RHS.

3. **Evaluation dominates at 96 threads — the documented FFTW pathology.**
   `four_index_matrix` builds the corrected layer-potential map once, then loops
   one `TKM3D.ltkm3dc` u_inc call per column; at 96 threads each call hits the
   ~47 s FFTW pathology (see the 2026-06-15 root-cause section above), so eval
   scales ~K×47 s and is **77 % of the K=19 production cost** (888 of 1160 s).
   This was measured as-shipped on purpose. Actionable corollary: capping the
   FINUFFT thread count (~16) in `ltkm3dc` would cut each call to ~1.3 s, taking
   the K=19 eval from ~888 s to ~25 s and the production total from ~1160 s to
   ~300 s — the single highest-leverage fix for the whole multi-RHS pipeline,
   independent of this benchmark.

4. **Peak RAM grows with K — the block-vs-sequential trade-off.** Peak RSS
   climbs 18.5 → 164 GB across the sweep, and the per-stage increments pin the
   jump to **block GMRES** (the K-blocked Krylov subspace over the ~960k-point
   operator). The sequential path stays near the precompute baseline. 164 GB at
   K=19 is comfortable on the 1.5 TB node, but the block speedup is bought with
   memory — relevant if K or the interface grows much larger.

5. **Correctness holds at every K.** Block vs sequential Σ agree to 1e-15 (K=1,
   identical system) and ~1e-5 (K≥4, at the gmres tol); the four-index matrix is
   symmetric to `max|V − Vᵀ|/max|V| ≤ 4.3e-6`; block-GMRES residual ~5e-6 at the
   1e-5 target in 7–8 iterations (no penalty vs single-RHS). The full K×K V for
   each cutoff is in the `.jls` records.

## Caveats

- **Onsite value differs from the single-RHS 9.88 eV by the z-convention.**
  `u_onsite` is 7.05–7.09 eV across all cutoffs (stable, as it should be — the
  onsite element is cutoff-independent). The ~7 eV vs the single-RHS 9.88 eV is
  the **midplane (here) vs top-face (`run_single_rhs.jl`) orbital placement**;
  the eps_in = 3.5 calibration shifts with it (see the single-RHS caveat). The
  *cost* numbers are directly comparable to the single-RHS REPORT (same params,
  same ~960k interface); the eV value is not, by design.
- The CSV/jls also record `t_total` (full `run_pipeline` wall, which *includes*
  the sequential baseline) — the production figures above use `precompute +
  block + eval`, not `t_total`.
- Re-running appends to `data/multi_rhs.csv`; delete it first to re-measure.

## UPDATE (2026-06-16): the precompute/eval cost was the FINUFFT 96-thread FFTW pathology — fixed in-function

The OMP=96 numbers above are **inflated by the FFTW 96-thread plan pathology** (the same
one the single-RHS report found in the eval). Investigation (`scripts/diag_*.jl`):

- **The interface build is NOT algorithmic growth — it's the FFTW thread pathology.**
  Build of the cutoff-2.5 envelope: **171 s @96 threads → 31 s @64 → 17 s @16**
  (`diag_build_times.jl`). The lever is `OMP_NUM_THREADS` (FINUFFT/FFTW) alone — Julia and
  OpenBLAS thread counts are irrelevant. A sampling profile mislabels it (FINUFFT's FFT runs
  in C/OpenMP worker threads the Julia profiler can't attribute), but the thread-lever test
  is decisive.
- **Direct FINUFFT micro-benchmark** (`diag_finufft_fftw.jl`, type-1 on the pathological
  413×375×325 grid): exec **0.70 s @16 threads vs 27.66 s @96** with `FFTW_ESTIMATE`
  (FINUFFT's default). `FFTW_MEASURE` fixes the exec (0.51 s @96) but its *planning* costs
  25 s @16 / 54 s @96 — a net loss unless wisdom is persisted, and useless for the build's
  varying mode grids. A good plan execs in ~0.5 s at **any** thread count, so capping the
  FFT threads loses nothing.
- **`n_src` growth is a truncation effect, not new support** (`diag_support.jl`): the union
  of pair supports is always ⊆ supp(φ_center); n_src grows 137k→438k because the rss
  envelope retains φ_center's faint tail (where the neighbor pairs, linear in φ_center, are
  ~36× the quadratic onsite) — bounded by |supp(φ_center)|. Genuine build growth at a sane
  thread count is only ~6 s.

**Fix (TKM3D, branch `finufft-nthreads-cap`): cap FINUFFT's `nthreads` IN-FUNCTION inside
`ltkm3dc`** (default 16, env `TKM3D_FINUFFT_NTHREADS`), threaded into the type-1 `nufft3d1`
and type-2 plans. This caps *only* the FFT; the FMM (`lfmm3d`, a separate library) keeps the
full OpenMP pool — strictly better than a global `OMP` cap, which throttles the FMM and (at
64) still leaves the eval crippled.

Re-run at the production env (JULIA=96, OMP=96) with the in-function cap (worker7122,
`data/multi_rhs.csv`; OMP=96 baseline preserved as `data/multi_rhs_omp96.csv`):

| K | precompute 96 → fix | eval 96 → fix | block solve 96 → fix | prod total (pre+blk+eval) 96 → fix |
|--:|--:|--:|--:|--:|
| 1 | 38.9 → 40.2 | 64.5 → 68.0 | 17.5 → 18.1 | 121 → 126 |
| 10 | 149.2 → **38.4** | 456.8 → **172.6** | 62.3 → 62.6 | 668 → **274** |
| 13 | 145.3 → **38.1** | 599.1 → **197.3** | 75.5 → 81.6 | 820 → **317** |
| 19 | 162.2 → **38.5** | 888.2 → **277.5** | 110 → 114 | 1160 → **430** |

- **Build flat at ~21 s** across all K (was 131–145 s) — precompute now K-independent (~38 s).
- **Per-K eval 46 → 12 s/K** (the ~34 s/call ltkm3dc pathology removed). The residual eval is
  now **FMM-bound** (the corrected `pottrg` map build + per-column apply), not FFT — next
  optimization is a batched `nd=K` corrected `pottrg`.
- **FMM-based stages unchanged** (block solve 110→114 s) — the cap hit only the FFT.
- **Correctness identical**: seq_agree 1e-15→1e-5, V-symmetry ≤4.3e-6, u_onsite 7.05–7.09 eV.
- End-to-end block path at K=19: **2.7× faster** (1160→430 s); precompute 4.2×, eval 3.2×.

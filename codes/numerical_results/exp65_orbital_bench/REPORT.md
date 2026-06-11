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

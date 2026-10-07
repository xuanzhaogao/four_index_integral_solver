# Experiment 6.1 — convergence study: results and findings

Date: 2026-06-10. Hardware: worker7018 (96-core AMD Genoa, 1.5 TB), 96 threads.
All solves run by hand on the worker; scripts under `scripts/`, raw data under
`data/raw/`, figures under `figs/`. Everything below uses the production
pipeline (per-point screened source -> RHS-adaptive multi-box interface ->
FMM operator with sparse near corrections -> GMRES -> TKM + corrected
evaluation), with **all component tolerances tied to a single eps**.

## Systems

| | slab | Fig.-1 system |
|---|---|---|
| geometry | 10 x 10 x 1 slab, eps 10, vacuum outside | 10 x 10 x 1 slab (eps 10) resting on two 10 x 10 x 10 cubes (eps 4 / eps 12) sharing the face x = 0 |
| source | Gaussian s = 0.05 at slab center, support fully inside | same, at slab center (0, 0, 0.5) |
| V target | identical Gaussian displaced (0.2, 0, 0) | identical, displaced laterally 0.2 |
| phi zones | 200 pts each: 0.01 above top face / source-support ball / far sphere | same construction |
| special features | edge/corner singularities only | + two slab-cube contact faces, buried cube-cube face, triple-junction lines (slab-cube-vacuum, slab-Omega1-Omega2) |

Sweep: eps = 1e-4 fixed; p in {2, 4, 6} x edge depth r in 1..5 (slab also r = 6).
References (differ in every discretization parameter, edge-corrected):
p = 8, eps = 1e-6, r = 6, finer source/TKM grid (margin 1.4).
Slab: V_ref = 0.071804515761 (N = 3.29M, 16 iters, 149 s).
Fig.-1: V_ref = 0.055678122198 (N = 11.84M, 33 iters, 1134 s).

## Finding 1 (central): edge-pair correction is load-bearing

With the production default `correct_edges = false`, edge panels are excluded
from near-field correction (both in the operator and in the target-evaluation
corrections). Consequences, measured on the slab:

- V-error decays only ~2x per edge-refinement level (first order in l_ec),
  identically for every p — the edge-region quadrature error dominates all
  metrics (V, every phi zone).
- The p-ordering INVERTS: p = 6 is worst (0.20 rel. error at r = 1!), p = 2
  "best" — at p = 2 the massive upsampled near-pair set (~3e5 pairs)
  incidentally covers edge neighborhoods, while at p = 6 the classifier finds
  nothing to upsample and edge pairs go fully uncorrected.
- The original reference run (p = 6, r = 6, eps = 1e-6, no edge correction)
  was the LEAST converged number in the suite: V = 0.072404 vs true 0.071805
  (+0.83%). All first-pass "errors" were measuring the reference's own error.
- Enabling `correct_edges = true` (adaptive quadtree on edge-touching pairs):
  30-80x error reduction at equal mesh (p = 6, r = 3: 4.4e-2 -> 8.4e-4), and
  GMRES iterations drop ~2x (22-25 -> 10-12) — it improves conditioning too.

Verification of the converged value: edge-corrected r-series extrapolates to
0.071802; uncorrected p = 2 series descends from above to ~0.07183 at r = 7;
the independent p = 8 edge-corrected reference lands at 0.0718045. Tightening
eps 1e-4 -> 1e-6 at fixed mesh moves V by only 2e-5, so the tied tolerance is
not the binding error in this regime — r is.

**Action taken:** `correct_edges = true` is now the harness default for the
operator AND the evaluator's target corrections (the stock BI evaluator
excludes edge panels; our `eval_scatter` includes them). The no-correction
sweep is archived in `data/raw/noedges_run1/` and plotted as the dashed
contrast in the convergence figure.

## Finding 2: clean geometric convergence with a universal error-vs-DOF curve

With edge correction on, all p-series converge monotonically at ~0.6x per
edge level, and the curves collapse: p = 4 at level r matches p = 6 at level
r-1 to ~3% on both systems. Interpretation: the kernel quadrature is fully
corrected; the remaining error is the RESOLUTION of the sigma edge/corner
singularity, which is controlled by the graded-mesh depth r, with p buying
exactly one level per doubling. Headline numbers (V rel. error):

- slab: p=2: 1.0e-2 -> 1.5e-3 (r 1..5); p=4: 4.1e-3 -> 4.4e-4; p=6: 2.5e-3 -> 2.4e-4
- Fig.-1: p=2: 4.0e-3 -> 6.3e-4; p=4: 1.9e-3 -> 2.0e-4; p=6: 1.2e-3 -> **9.9e-5
  (reaches the tied eps = 1e-4 target, N = 3.3M, 258 s)**

phi-zone errors track V with fixed ratios (support ~0.5x, near-interface
~3-5x, far ~4-6x) and the same decay rate.

## Finding 3: no degradation on the three-material contact geometry

The Fig.-1 system — contact faces, buried junction, multiple coexisting
contrasts (gamma = 0.6, 0.85, -0.5, +0.43, -0.27, +0.82 across its six
interface types) — converges with the same rate and constant structure as the
single slab. N_iter: 16-24 across the sweep (33 at the 11.8M reference) vs
7-12 for the slab; growth with N is ~logarithmic in both.

## Finding 4: cost ordering favors high p

p = 6 runs are simultaneously the most accurate and the fastest (slab: 17-46 s
vs 67-84 s at p = 2). p = 2 is a bad operating point: n_quad = 2 forces
upsampling corrections on ~3e5 panel pairs, so near-block assembly dominates
its wall time while accuracy per DOF is no better than p = 6.

## Finding 5: correctness checks (Fig.-1 system, eps = 1e-4, p = 6, r = 4)

- Geometry: census symmetric (cube|vac 14076 each; contacts 2204 each;
  cube-cube face 2968, emitted once); zero duplicate panels.
- Transmission: [phi] continuity 4e-9 (x = 0 face) down to 2.5e-12 (contact
  faces); [eps dn phi] flux condition 2.8e-4 .. 2.2e-3 relative to the
  face-wide flux scale (the 4-node FD probe floors near eps/delta = 2e-2;
  measured values sit an order BELOW that).
- Conservation: global monopole int rho_scr + sum sigma = 0.998387
  (dev 1.6e-3, consistent with the r = 4 edge-resolution level since it sums
  raw sigma quadrature on singular edge panels); Gauss flux through an
  enclosing box = -1.000003 (**dev 3.1e-6**, evaluated through the corrected
  phi — the honest measure of field quality).
- Per-body sum(sigma dS) is NOT an invariant for touching bodies (bound charge
  lives on shared interfaces); reported as diagnostics only. The isolated-body
  version of the check remains valid for gap geometries (6.2's box).

## Technical findings (pipeline-level, affect later experiments)

1. **TKM far-target blowup:** TKM's Fourier domain covers sources AND targets;
   asking it for far-field points inflated the grid ~1150^3 and hung. Fix in
   `Harness.eval_incident`: targets outside the source-support bounding box
   (+3 spacings) use direct threaded summation (exact to quadrature accuracy
   there); only near/inside-support targets go through TKM.
2. **Batched evaluation:** V-grid + all zone sets are evaluated through ONE
   post-refined operator per run (one refinement + one FMM + one hcubature
   assembly) — ~4x cheaper post-eval.
3. FMM3D accepts tol down to 1e-14 (saturates ~3e-14 rel.) — no floor issue
   for any planned tolerance.
4. Julia footgun caught in the FD checks: `4f2` parses as Float32 4e2;
   FD weights are now built from a Vandermonde solve, never literals.
5. Source grids follow n(eps) = ceil(margin * (4/pi) * sqrt(ln(10/eps) ln(1/eps)))
   per dimension (midpoint-rule aliasing <= eps); reference margin 1.4.

## Data and figures

- `data/sweep_slab.csv`, `data/slab_errors.csv`, `data/raw/slab_*.jls`
  (+ archived no-edge-correction sweep in `data/raw/noedges_run1/`)
- `data/sweep_fig1.csv`, `data/fig1_errors.csv`, `data/raw/fig1_*.jls`,
  `data/raw/fig1_checks.jls`
- `figs/fig61_slab_convergence.pdf` (with the no-edge-correction contrast),
  `figs/fig61_fig1_convergence.pdf`,
  `figs/fig61_system_slab_geometry.pdf`, `figs/fig61_system_fig1_geometry.pdf`
  (RHS-adaptive meshes, source-refinement insets)

Compute time (wall, 96 threads): slab sweep 924 s + no-edge archive 1217 s +
Fig.-1 sweep 2779 s + checks ~360 s + diagnostics ~600 s ~= **1.6 h total**.

## Implications going forward

- All later experiments run edge-corrected by default; 6.2 gains variant
  "Bp" (upsampling on / edges off) to separate the two correction mechanisms.
- 6.4 Part A should report the adaptive edge-block assembly cost explicitly
  (uniform refinement multiplies edge-touching pairs ~4x per level).
- Reference protocol for the paper: edge-corrected, p one step above tests,
  eps 100x tighter, r two levels deeper, source grid margin 1.4 — and verify
  by extrapolation bracketing as done here.

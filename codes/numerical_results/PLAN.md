# Revised plan: Section 6 numerical experiments

Revision of `prompt.md` after reading the production monolayer/graphene pipeline
(`codes/graphene/monolayer/*`, `codes/graphene/bilayer_slab/src/ScreenedOrbitalSolve.jl`)
and the solver internals (`~/codes/BoundaryIntegral.jl`). The experiment *definitions*
(geometries, sweeps, metrics) are kept exactly as specified in `prompt.md`; the revision
fixes **how** they are run so every experiment exercises the same code path as the
physics application, and documents parameter mappings, solver-capability facts, and
feasibility constraints discovered in the code.

## 1. What changed vs prompt.md, and why

1. **All experiments run the production pipeline** used by the graphene application
   (`ScreenedOrbitalSolve.solve_screened_mode`):
   `screened source -> *_dielectric_box3d_rhs_adaptive -> Rhs_..._fmm3d ->
   Lhs_..._fmm3d_corrected -> Krylov.gmres -> (TKM incident + corrected-FMM scattered)
   target evaluation`. The shared harness (`common/Harness.jl`) inlines this body
   (same pattern as `monolayer/scripts/timing_density_solve.jl`) so each phase is
   timed and the near-correction structures are introspectable.
2. **Per-point dielectric screening of the source density.** With s = 0.05 Gaussians,
   the truncated support radius is 3.7s–7.1s (depends on truncation tol), which
   *crosses interfaces* in Systems I/II and 6.2 (standoff 0.1 = 2s, material-box
   half-height 0.1 = 2s). Scalar `eps_src` would be wrong; we scale the density
   per quadrature point by 1/eps(x) (multi-box-aware `SharpScreening` equivalent)
   and pass `eps_src = 1.0` to the RHS, exactly as the monolayer slab driver does.
3. **Instrumentation requires inlining the corrected LHS.**
   `Lhs_dielectric_box3d_fmm3d_corrected` returns an opaque `LinearMap`; near-pair
   count, p_up distribution and sparse-correction memory are obtained by calling
   `BI.build_neighbor_list` + `BI.laplace3d_DT_corrections` directly and assembling
   the operator in the harness.
4. **p_up cap fact (answers 6.4.3):** the upsampling order is clamped to `max_order`
   (`src/kernel/laplace3d_near.jl:360–367`); there is **no adaptive-quadrature
   fallback** for non-edge-touching pairs when the cap binds. Production uses
   `max_order = 64`; we adopt 64 for all runs and report when the cap binds.
5. **h0 is hard-coded** in `laplace3d_pottrg_fmm3d_corrected_hcubature`
   (`panel_size_limit = min panel length`). For 6.4B the harness calls
   `BI._refine_interface_for_targets` with an explicit `panel_size_limit = h0`
   (and `Inf` for the "no post-refinement" baseline), then assembles
   FMM + hcubature-near evaluation manually.
6. **Multi-box shared faces are emitted exactly once** with the correct eps pair
   (`src/shape/box3d_multi.jl:139–353`); 6.1.7's sanity check verifies this at run
   time (count panels on the x=0 plane, check the gamma value, plus the transmission
   check of 6.1.6 as the independent detector).
7. **FMM tolerance floor.** FMM3D cannot run at eps = 1e-13. Tied tolerance maps
   `fmm_tol = max(eps, 1e-12)`; for the eps = 1e-13 reference runs the FMM/near/GMRES
   tolerances clamp at 1e-12 while RHS-refinement/TKM run at 1e-13. This is the only
   deviation from the literal "all tolerances = eps" rule and is reported wherever it
   binds. (Verified during smoke test; if 1e-12 also fails the floor is raised and
   documented.)
8. **gamma -> 1 limit (6.3):** the code takes finite eps_in only; the
   "grounded-conductor" run uses eps_1 = 1e12 (gamma - 1 ~ 1e-12) and is labeled as such.
9. **Unspecified parameters in prompt.md, fixed here** (flagged, not silent):
   - 6.2 material box permittivity: eps_m = 2 (same as System II's material box).
   - 6.2/6.4 target density for V: same convention as System I — Gaussian, same
     width s = 0.05, displaced laterally by 0.2 from the source (inside the material box).
   - Edge-refinement depth at the tied-eps settings: r = 4 for all test runs (the
     production value) except where the sweep itself varies r (6.1.2) or the
     reference rule requires r+2.
   - GMRES: `rtol = eps`, `atol = 1e-14` (effectively rtol-bound), `itmax = 500`.
10. **6.4 feasibility staging:** Part A targets N from ~1e4 to ~3e6 guaranteed,
    1e7 attempted only if the wall clock allows (each 1e7 GMRES solve is
    O(hours)); Part B runs M in {1, 10, 100} fully, M = 1000 by strategy-(ii) full run
    + strategy-(i) subset extrapolation (as the prompt already allows); if (ii) at
    M = 1000 does not fit the night it is extrapolated and marked.
11. **Graphene application (future 6.5):** the realistic monolayer results
    (calibrated eps_in = 3.5, Lz = 3.35 A vs CoQui cRPA) already exist under
    `codes/graphene/monolayer/`; they are *not* re-run tonight. The paper's
    application subsection can cite them; out of scope here.

## 2. Parameter mapping (tied tolerance eps)

| prompt.md name | code parameter | value |
|---|---|---|
| eps_RHS | `rhs_atol` of `*_rhs_adaptive`, source truncation `tol` | eps |
| eps_near | `up_tol` (arg 2 of corrected LHS), `hcubature_atol` in evaluation | eps |
| (FMM) | `fmm_tol` (arg 1 of corrected LHS, RHS assembly, evaluation FMM) | max(eps, 1e-12) |
| eps_TKM | `volume_tol` of TKM incident evaluations | eps |
| eps_GMRES | `gmres_rtol` | eps |
| p | `n_quad` | {4, 6, 8}; ref 12 |
| edge depth r | `l_ec = (min box dim)/2^r * 1.01` | 4 unless swept; ref: deepest test +2 |
| p_up cap | `max_order` | 64 |
| h0 | `panel_size_limit` in target post-refinement | swept in 6.4B |
| rho* | `eps^(-1/2p)` (Bernstein criterion, diagnostics only) | derived |
| max_depth | `max_depth` | 128 |

DOF N = `BI.num_points(interface)`. N_iter = `stats.niter` from `Krylov.gmres`
(residual history from `stats.residuals` for the 6.3 inset).

**Per-run log (CSV row):** all input parameters + N, N_iter, near-pair count
(upsample + adaptive/edge), p_up max & mean, wall-clock split {source/TKM prep,
interface build, RHS assembly, near-block assembly, GMRES total, evaluation},
sparse-correction memory (`Base.summarysize` + nnz), V, phi at target sets,
GMRES residual, achieved errors vs reference. Raw rows append incrementally
(crash-safe); references stored as `.jls` snapshots with metadata.

**Reproducibility:** RNG seed 20260609 (target sets), 20260610 (6.4B centers);
target sets generated once per system and saved. Hardware (worker7018: 96-core
Genoa, 1.5 TB) and `JULIA_NUM_THREADS`/`OMP_NUM_THREADS` recorded once per CSV.

## 3. Geometry definitions (explicit coordinates)

**System I (6.1, 6.3):** unit cube center (0,0,0), eps_1 = 10 (default; swept in 6.3),
vacuum outside. Source Gaussian s = 0.05 at (0, 0, 0.6); target density Gaussian
s = 0.05 at (0.2, 0, 0.6). V = sum(w .* rho_tgt .* phi) over the target grid.
phi-error zones: (i) z = 0.51 plane, x,y ~ U(-0.5, 0.5); (ii) ball r <= 0.1 around
(0,0,0.6); (iii) sphere |x| = 5. 200 points each, fixed seed.

**System II (6.1.4–7, 6.3.3):** Omega_1 center (-0.5, 0, 0.25) size 1x1x0.5 eps 4;
Omega_2 center (0.5, 0, 0.25) size 1x1x0.5 eps 12 (shared face x = 0);
Omega_m center (0, 0, 0.85) size 0.6x0.6x0.2 eps 2 (gap g = 0.1 above substrate
top z = 0.5). Source Gaussian s = 0.05 at (0, 0, 0.85); target at (0.2, 0, 0.85).
Zones: (i) z = 0.51 plane over x,y ~ U(-1,1)x(-0.5,0.5); (ii) ball r <= 0.1 at
(0,0,0.85); (iii) sphere |x| = 5.

**6.2/6.4 slab:** slab AxAx0.5 center (0,0,0) eps 10 (A = 10 default); material box
0.6x0.6x0.2 center (0, 0, 0.4) eps 2 (gap 0.05 above slab top z = 0.25). Source
s = 0.05 at (0,0,0.4); target at (0.2, 0, 0.4).

## 4. Experiments (unchanged science, concrete recipes)

- **6.1.1/6.1.4** eps in {1e-3,…,1e-11} x p in {4,6,8}, r = 4. Reference per system:
  p = 12, eps = 1e-13 (FMM clamp note), r = 6, source grid n doubled (finer TKM).
  Outputs: V-err & phi-err vs DOF (3 p-curves, 3 zones), N_iter companion panel;
  Systems I & II side by side (1x2), N_iter table.
- **6.1.2** p = 8, eps = 1e-11 settings; sweep r = 0…8 alone. V-err vs r.
- **6.1.3 (optional, time-permitting)** d in {0.5, 0.1, 0.02}, eps = 1e-9: error/DOF/N_iter,
  DOF ~ log(1/d) check.
- **6.1.5** panelization figure colored by refinement level (top faces of System II)
  + sigma along a line approaching the triple junction (log-log).
- **6.1.6** transmission checks across Omega_1–Omega_2 and partial Omega_m:
  [phi] and [eps dn phi] via near-surface paired evaluations; sum(sigma dS) per body;
  Gauss law on an enclosing box. One table at eps = 1e-9, p = 8.
- **6.1.7** geometry sanity first: shared-face panel census + gamma_12 check.
- **6.2.1** variants A (full) / B (`Lhs_dielectric_box3d_fmm3d`, no correction) /
  C (uniform refinement DOF-matched ±10% to A, same l_ec) at p = 8,
  eps in {1e-3,…,1e-11}. Annotate B's predicted stagnation rho_min^(-2p) (min Bernstein
  rho over the slab-top/box-bottom pair set, from neighbor-list geometry).
- **6.2.2** side table at eps = 1e-9.
- **6.2.3** A in {1, 5, 10, 20, 50}, eps = 1e-9: V-err (A,B) + near-pairs & max p_up.
- **6.3.1** eps_1 in {2, 5, 10, 20, 80, 1000} (+1e12 as gamma = 1), eps = 1e-9, p = 8,
  per-contrast references; N_iter vs gamma; residual-history inset for {2, 80, 1000};
  V-err vs eps_1.
- **6.3.2** V/V_vacuum vs eps_1 (V_vacuum = direct TKM between the two Gaussians,
  no interface), saturation vs the gamma = 1 run.
- **6.3.3** System II with (eps_1, eps_2) scaled by {1, 2, 8}: N_iter table.
- **6.4A** slab problem, eps = 1e-6, p = 8; N grown by global uniform refinement
  (+ k x k tiling of the slab+box unit if needed, k <= 3): per-iteration FMM matvec,
  per-iteration near-block apply, setup (pair search + block assembly), GMRES total;
  fitted slopes; correction-memory vs N; p_up distribution; cap report.
- **6.4B** fixed sigma from the 6.2 run at eps = 1e-6: M in {1, 10, 100, 1000} Gaussian
  target densities (s = 0.05, centers ~ U inside Omega_m, seed fixed), grid 100^3 each.
  (i) no post-refinement (subset >= 10 + linear extrapolation for large M, marked);
  (ii) post-refinement at h0 + FMM far + hcubature near. Validation vs
  hcubature rtol 1e-12 on 10 targets. Total & per-target time vs M; crossover;
  error-vs-cost h0 sweep (5 values, M = 100).

## 5. Layout & deliverables

```
numerical_results/
  PLAN.md  prompt.md  Project.toml
  common/Harness.jl          # instrumented pipeline, screening, targets, logging, theme
  exp61_convergence/{scripts,data,figs,README.md}
  exp62_ablation/{...}  exp63_contrast/{...}  exp64_scaling/{...}
```

Runs execute on worker7018 (ssh, absolute juliaup binary, JULIA_NUM_THREADS=96),
launched as background jobs with per-script logs under each experiment's `data/raw/`.
Figures: CairoMakie, colorblind-safe (Okabe–Ito), consistent markers,
single-column (~3.4 in) or 1x2 (~7 in) PDF. Each README states geometry, parameters,
reference settings, hardware, and total compute time. Any infeasible run is reported,
not silently altered.

# Task: Numerical experiments for Section 6 of the paper

You are working with an existing adaptive BIE solver for screened Coulomb four-index integrals in piecewise-constant dielectric environments. The solver components (already implemented and unit-tested) are: TKM volume solver for the incident potential, Nyström discretization with panelwise p×p Gauss–Legendre quadrature, RHS-driven dyadic panel refinement (tolerance ε_RHS), geometry-driven edge refinement (threshold δ), Bernstein-radius near-field correction with upsampled blocks (tolerance ε_near, threshold ρ* = ε^(−1/2p)), FMM-accelerated GMRES, and post-refinement for repeated target evaluations (threshold h0). Your job is to set up, run, and produce publication-quality data and figures for experiments 6.1–6.4 below.

## Global conventions

- **Tied tolerance.** All component tolerances are controlled by a single parameter ε: set ε_RHS = ε_near = ε_TKM = ε_GMRES = ε. Every sweep over "tolerance" means sweeping this single ε.
- **Reference solutions** are self-convergence references. The reference run must differ from all test runs in every discretization parameter simultaneously: p = 12, ε = 10⁻¹³, edge refinement two levels deeper than the deepest test run, and a finer TKM volume grid. Never use a test run as a reference for another test run.
- **Error metrics.** Relative error of the final integral V; relative error of φ at fixed sets of ~200 random off-grid target points (sets defined per experiment); both against the reference run.
- **Always log per run:** DOF N, GMRES iteration count N_iter, number of near-correction pairs, max and mean upsampling order p_up, wall-clock time split into {RHS/TKM, refinement+setup, near-block assembly, GMRES total, evaluation}, peak memory of the sparse correction matrix.
- **Reproducibility.** Fixed RNG seeds for target-point sets and source placements. Save all raw data (JSON/CSV) separately from plotting scripts. Record hardware and thread count once.
- All densities are normalized isotropic Gaussians unless stated otherwise.

## 6.1 Convergence study (two systems)

**System I — single dielectric box.** Unit cube, ε₁ = 10, in vacuum. Source: Gaussian, width s = 0.05, at standoff d = 0.1 above the top face center. Target density: Gaussian of the same width, displaced laterally by 0.2 at the same height.

1. Sweep ε ∈ {10⁻³, 10⁻⁵, 10⁻⁷, 10⁻⁹, 10⁻¹¹} × p ∈ {4, 6, 8}. Report error of V and of φ in three target zones: (i) near-interface (dist 0.01 from the top face), (ii) inside the source support, (iii) far field (dist 5). Output: error vs DOF plot, three p-curves; N_iter vs DOF in a companion panel.
2. Edge-refinement isolation: fix p = 8 and all parameters at the ε = 10⁻¹¹ setting; sweep edge refinement depth r = 0…8 alone. Output: V-error vs r, showing geometric decay then flooring.
3. (Optional) Standoff sweep d ∈ {0.5, 0.1, 0.02} at fixed ε = 10⁻⁹: report achieved error, DOF, N_iter; verify DOF grows ~ log(1/d).

**System II — three-material configuration (paper Fig. 1).** Two coplanar substrate boxes Ω₁ (ε₁ = 4) and Ω₂ (ε₂ = 12), each 1×1×0.5, sharing one internal face; material box Ω_m (ε_m = 2), 0.6×0.6×0.2, centered over the junction line at gap g = 0.1; Gaussian source (s = 0.05) inside Ω_m.

4. Repeat sweep (1) with identical metrics. Output: a 1×2 figure with System I and System II convergence side by side (same axes), demonstrating no degradation of convergence order; N_iter for both systems in one table.
5. Qualitative figure: final adaptive panelization colored by refinement level; plus σ sampled along a line approaching the triple-junction edge, log-log.
6. Correctness checks: sample φ and ε∂ₙφ along lines crossing the Ω₁–Ω₂ interface and ∂Ω_m; verify both transmission conditions hold to ε. Also verify ∮σ dS = 0 on each closed body and Gauss's law on a box enclosing the source. Report as one small table.
7. **Geometry sanity check (must run first):** confirm the shared Ω₁–Ω₂ face enters Γ exactly once with γ₁₂ = (ε₁−ε₂)/(ε₁+ε₂). If geometry is built per-box and unioned, check explicitly for double-counted panels; the transmission check in (6) is the detector.

## 6.2 Component comparison (ablation on a high-aspect-ratio slab)

**Geometry (deliberately near-correction-stressing):** dielectric slab 10×10×0.5, ε = 10; material box 0.6×0.6×0.2 with Gaussian source at gap g = 0.05 above the slab center.

Three solver variants, identical in everything else:
- (A) full method;
- (B) no near correction: set all D^{T,corr}_{PQ} = 0;
- (C) no RHS adaptivity: replace the RHS criterion with uniform dyadic refinement matched to within 10% of variant A's DOF; keep edge refinement identical.

1. Sweep ε ∈ {10⁻³ … 10⁻¹¹} at p = 8 for all three variants. Output: V-error vs ε (or DOF), three curves. Annotate the predicted stagnation level of variant B from the Bernstein estimate ρ_min^(−2p) of the dominant slab-face pair as a dashed line.
2. Side table at ε = 10⁻⁹: DOF, N_iter, near-pair count, setup time, total matvec time, achieved V-error for A/B/C.
3. Aspect-ratio sweep at fixed ε = 10⁻⁹: slab A×A×0.5 with A ∈ {1, 5, 10, 20, 50}, gap fixed. Output: 1×2 figure — left: achieved V-error for variants A and B vs aspect ratio; right: near-pair count and max p_up vs aspect ratio.

## 6.3 Dielectric contrast

**Geometry:** System I of 6.1 (unit cube + external Gaussian source), fixed ε = 10⁻⁹, p = 8.

1. Sweep ε₁ ∈ {2, 5, 10, 20, 80, 1000}; the last value probes the γ→1 conductor limit. Report: N_iter vs γ; GMRES residual histories for ε₁ ∈ {2, 80, 1000} as an inset; achieved V-error vs ε₁ against per-contrast references (accuracy must not degrade with contrast).
2. Physics curve: V/V_vacuum vs ε₁; verify monotone screening that saturates as γ→1, and compare the saturation value against a γ = 1 run (grounded-conductor limit).
3. Repeat (1) for the three-material System II varying (ε₁, ε₂) jointly, e.g. scaling both by {1, 2, 8} relative to baseline, to show iteration behavior with multiple coexisting γ values. Compact: one N_iter table.

## 6.4 Scaling and post-refinement amortization

**Part A — solver scaling.** Take the 6.2 slab problem at ε = 10⁻⁶, p = 8. Grow N from ~10⁴ to ~10⁶–10⁷ by uniform refinement beyond what adaptivity requires (and/or by tiling k×k copies of the slab+box unit, k = 1, 2, 3, 4).

1. Output: log-log wall-clock vs N with separate curves for FMM matvec (per iteration), near-correction block apply (per iteration), setup (KD-tree pair search + upsampled block assembly), and GMRES total. Annotate fitted slopes. FMM matvec should be O(N) or O(N log N); explicitly characterize where (if anywhere) the correction cost ceases to be subdominant.
2. Report sparse-correction memory vs N, and the distribution (mean/max, or histogram) of p_up.
3. Document the implementation's p_up cap and the fallback behavior when the cap is hit; if no cap exists, flag this.

**Part B — post-refinement amortization.** Fix one solved σ from the 6.2 slab problem. Generate M ∈ {1, 10, 100, 1000} target densities (Gaussians, random centers inside Ω_m, fixed seed), each on an n³ ≈ 10⁶ point grid.

4. Compare: (i) no post-refinement — all near interactions via adaptive cubature per target (for large M, time on a subset of ≥10 targets and extrapolate linearly; mark extrapolated points); (ii) post-refinement (h0 fixed) + FMM far field + residual near set via adaptive cubature.
5. Accuracy: on a 10-target subset, validate both strategies against high-tolerance adaptive cubature (rtol 10⁻¹²).
6. Output: total time vs M (log-log; strategy (ii) should be setup + c·M with small c; mark the crossover M); amortized per-target time vs M; and at fixed M = 100, an error-vs-cost curve sweeping h0 over 4–5 values (this justifies the choice of h0 in the paper).

## Deliverables

- One directory per experiment with: run scripts, raw data, plotting scripts, generated PDFs.
- Figures sized for a single JCP column unless 1×2; consistent fonts/markers across all figures; colorblind-safe palette.
- A short README per experiment stating geometry, parameters, reference-run settings, and total compute time.
- Do not silently change any specified parameter; if a run is infeasible (memory/time), stop and report which one and why.


position of the article: /mnt/home/xgao1/Articles/four_indices_bie
code: /mnt/home/xgao1/codes/BoundaryIntegral.jl
skills: use-cluster, use julia
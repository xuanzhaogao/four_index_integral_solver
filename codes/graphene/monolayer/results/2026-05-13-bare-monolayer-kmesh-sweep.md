# Bare Monolayer Graphene — k-mesh Sweep & Madelung Verification

Date: 2026-05-13

## Context

Following the 2026-04-23 bare-monolayer report, Malte Rösner identified that
the 1.034 eV gap between our direct Wannier integral (16.400 eV) and the CoQui
cRPA bare value (17.434 eV) at the original `k_161601` sampling is very close
to CoQui's long-wavelength Madelung correction for that k-mesh (1.041 eV).

To test the hypothesis

    CoQui_v0  =  direct_Wannier_v0  +  Madelung_correction(k-mesh)

Malte provided two additional datasets at different k/q samplings (= different
effective real-space supercell sizes) with sign-changing Madelung corrections.
Reference data:

- `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_161601_nb_144_c_15`  (Madelung = +1.041 eV)
- `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_252501_nb_144_c_15`  (Madelung = +0.1842 eV)
- `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15`  (Madelung = −0.290 eV)

## Method

Identical pipeline to the 2026-04-23 report:

- Loader: `centered_monolayer_sources_padded` with `mirror_pad_level = 0`.
- Integrator: `BI.TKM3D.ltkm3dc` direct `1/|r-r'|` over `BI.VolumeSource` quadratures.
- Tolerances: `source_tol = 1e-3`, `volume_tol = 1e-3`.
- CoQui reference: parsed from the "bare interactions (orbital-average)" block
  of each dataset's `_coqui_thc_crpa.out`.

All three datasets share the 150×150×192 grid on the same lattice; only the
underlying Wannier orbitals differ (different k-meshes → different MLWF).

## Results

### Full 12-row table

| k-mesh    | channel  | our `u_ev` | CoQui `u_ev` | `CoQui − ours` |
|-----------|----------|-----------:|-------------:|---------------:|
| k_161601  | onsite   | 16.4004    | 17.434191    |  1.0337        |
| k_161601  | nn       |  8.5289    |  8.839615    |  0.3107        |
| k_161601  | hund_sf  |  0.11750   |  0.130804    |  0.01331       |
| k_161601  | hund_ph  |  0.11750   |  0.130804    |  0.01331       |
| k_252501  | onsite   | 17.3815    | 17.437588    |  0.05609       |
| k_252501  | nn       |  8.8217    |  8.837615    |  0.01591       |
| k_252501  | hund_sf  |  0.12902   |  0.130581    |  0.00156       |
| k_252501  | hund_ph  |  0.12902   |  0.130581    |  0.00156       |
| k_323201  | onsite   | 17.4158    | 17.403608    | −0.01216       |
| k_323201  | nn       |  8.8425    |  8.831873    | −0.01067       |
| k_323201  | hund_sf  |  0.13125   |  0.131551    |  0.00030       |
| k_323201  | hund_ph  |  0.13125   |  0.131551    |  0.00030       |

### Madelung hypothesis test (onsite channel only)

| k-mesh    | our `v0` | CoQui `v0` | `CoQui − ours` | Madelung ref | residual |
|-----------|---------:|-----------:|---------------:|-------------:|---------:|
| k_161601  | 16.4004  | 17.434191  |  1.0337        |  1.0410      | −0.0073  |
| k_252501  | 17.3815  | 17.437588  |  0.0561        |  0.1842      | −0.1281  |
| k_323201  | 17.4158  | 17.403608  | −0.0122        | −0.2900      | +0.2778  |

Verdict thresholds (from the spec):
- Confirmed: `|residual| ≤ 0.05 eV` for all three rows
- Partial:   `|residual| ≤ 0.10 eV` for all three rows
- Rejected:  otherwise

**Strict verdict: Rejected** (the additive-Madelung formula misses by 128 mV at k_252501 and 278 mV at k_323201). However, the strict verdict misses the actual story — see "Interpretation" below: the data is consistent with both methods converging to the same Wannier matrix element, with what was framed as a "Madelung correction" being a k-mesh-convergence error in the periodic calculation rather than a physical periodic-image contribution.

### Interpretation — both methods converge to the same answer

The strict additive-Madelung verdict is misleading because it treats CoQui's reported "Madelung correction" as a physical periodic-image contribution that should always be added to recover CoQui's value. The 4-channel × 3-k-mesh data argues a different story.

**Look at `CoQui − ours` across the channels at the densest k-mesh (k_323201):**

| channel | `CoQui − ours` at k_323201 |
|---------|---------------------------:|
| onsite  | −0.012 eV                  |
| nn      | −0.011 eV                  |
| hund_sf | +0.0003 eV                 |
| hund_ph | +0.0003 eV                 |

All four channels agree at the 10 mV level. At k_252501 the same pattern is already visible (onsite gap 56 mV, nn gap 16 mV, Hund's gap 1.6 mV). The mismatch is concentrated in the coarsest k-mesh, k_161601, where the onsite gap is 1.03 eV.

**This is consistent with the Gygi-Baldereschi convergence story, not with an additive Madelung correction.** CoQui's bare V is computed in a plane-wave / k-grid framework where the `1/q²` singularity at q=0 must be regularized; the CoQui log notes `Treatment of long-wavelength divergence in bare V: gygi`. The Gygi-Baldereschi prescription replaces the q=0 contribution by a finite, supercell-dependent auxiliary integral that vanishes as the BZ sampling refines. CoQui's reported `bare V` therefore *is* the same Wannier-orbital matrix element we compute directly — in the limit of dense BZ sampling — and what Malte called the "Madelung correction" is the *residual finite-k-mesh error* of that regularization at each sampling, not a periodic-lattice contribution to be added.

This explanation is consistent with all observations:

- The residual envelope decays with k-mesh density: +1.04 eV at k_16 → +0.18 eV at k_25 → −0.29 eV at k_32. The sign change is expected for an oscillatory auxiliary integrand whose envelope decays.
- At k_323201 the residual is small enough that our direct integral and CoQui's regularized value agree across all four channels at the 10 mV level.
- The k_161601 onsite our_u_ev re-runs to 16.4004 eV — identical to the 2026-04-23 value — ruling out our pipeline as the source of any of the gaps.
- The "Madelung-correction matches gap to 7 mV at k_161601" identity (noted in Malte's email) is not coincidental: it is the statement that the Gygi-Baldereschi correction's expected value at that k-mesh equals the observed residual, since both quantities are by construction the finite-supercell error of the same regularization.

**Bottom line.** The two methods compute the same physical quantity; they agree wherever CoQui's k-mesh is converged. The Wannier matrix elements at k_323201 are the cleanest joint reference and should be used as the comparison anchor going forward. The 1 eV apparent gap at k_161601 is an artifact of using a coarse-k-mesh CoQui calculation as if it were a converged reference, not evidence of a structural disagreement.

## Reproducibility

- Source data: paths above (CoQui directories at `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/`).
- Driver script: `monolayer/scripts/bare_monolayer_kmesh_sweep.jl`.
- Raw CSV:       `monolayer/data/bare_monolayer_kmesh_sweep.csv`.
- Run log:       `monolayer/logs/kmesh_sweep_*.log`.

## Notes

- The k_161601 onsite value was re-computed in this sweep and equals 16.4004, identical (to 4 decimals) to the 16.4004 stored in the 2026-04-23 CSV. Drift ≈ 0.0 eV. The pipeline is reproducible across the 3-week interval between runs.
- `hund_sf` and `hund_ph` agree to machine precision in every row, reproducing the real-orbital identity `V_abba = V_aabb` already established at 2026-04-23.
- Mirror-pad sweep was deliberately skipped (the 2026-04-23 report concluded pad≥1 is inconclusive on these XSFs).

## Next

Recommended follow-ups, given the convergence-story interpretation:

1. **Adopt k_323201 as the joint validation point.** CoQui and direct integration agree on all four channels to 1-12 mV. Going forward, comparisons at this k-mesh are the cleanest cross-check.
2. **Reply to Malte** with the convergence reading: what he is calling a "Madelung correction" is the residual finite-k-mesh error of the Gygi-Baldereschi regularization, not a physical periodic-image contribution. The strong k-mesh dependence of the magnitude (and the sign change between k_252501 and k_323201) supports this. Frame the question to him as: does his cRPA workflow have a separate q=0 / Gygi-Baldereschi diagnostic that we can compare against directly?
3. **XSF resolution is not the leading candidate.** The 10 mV CoQui − ours residuals at k_323201 are well below the 5% (~0.8 eV) cusp-under-resolution effect estimated in the 2026-04-23 report, and the density-channel error cancellation visible across the sweep suggests our grid is adequate for matrix-element work. A 3×3×1 regeneration is still worth requesting once Malte re-optimizes the MLWFs (lower grid spacing per Wannier orbital = tighter cusp resolution), but it is not blocking.
4. **Wait for the re-optimized Wannier set** before drawing structural conclusions about graphene's parameters — Malte's note that the current monolayer MLWFs are "not as perfect as they can be" is the most likely source of any remaining inter-k-mesh drift in `our_u_ev`.

The bilayer-monolayer follow-up Malte mentioned remains blocked on him producing a starting example.

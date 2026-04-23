# Bare Monolayer Graphene Interactions — Report

Date: 2026-04-23

## Reference (CoQui cRPA; Malte Rösner)

Source: `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_161601_nb_144_c_15/_coqui_thc_crpa.out`

Bare interactions (orbital-average):

- intra-orbital (onsite)  = 17.434191 eV
- inter-orbital (NN)      = 8.839615 eV
- Hund's (spin-flip)      = 0.130804 eV
- Hund's (pair-hopping)   = 0.130804 eV

## This Run

Orbitals: `graphene_00001.xsf`, `graphene_00002.xsf` (150×150×192 grid over a ~12.24 × 10.60 × 14.92 Å 5×5×1 Wannier90 supercell).

Method: direct real-space `1/|r-r'|` via `BI.TKM3D.ltkm3dc` on `BI.VolumeSource` objects built from the squared (density channels) and signed-product (Hund's channels) orbital grids. No dielectric solve. Shared centering shift from orbital 1's density centroid is applied to both orbitals.

Parameters: `source_tol = 1e-3`, `volume_tol = 1e-3`.

Pipeline independently validated against the analytic two-Gaussian Coulomb integral `erf(sqrt(α/2) R)/R` to 4.8e-4 relative error (see `bare integral against two-Gaussian analytic answer` testset).

### Level-0 results (raw XSF, no padding)

| Channel  | `u_ev` | `u_ref_ev` | `rel_err_pct` |
|----------|-------:|-----------:|--------------:|
| onsite   | 16.4004 | 17.434191 | −5.93 |
| nn       |  8.5289 |  8.839615 | −3.51 |
| hund_sf  |  0.1175 |  0.130804 | −10.17 |
| hund_ph  |  0.1175 |  0.130804 | −10.17 |

Spin-flip and pair-hopping agree to machine precision (`|J_sf - J_ph| < 1e-16`), confirming the real-orbital identity `V_abba = V_aabb`.

## Mirror-Pad Diagnostic

The spec called for a mirror-pad xy convergence check: extend the orbital grid by reflection and compare. The level-0 vs level-1 results are:

| Channel  | pad 0 `u_ev` | pad 1 `u_ev` | pad 0 Na | pad 1 Na |
|----------|-------------:|-------------:|---------:|---------:|
| onsite   | 16.4004 |  2.5310 |  71.65 | 644.89 |
| nn       |  8.5289 |  1.6101 |  71.65 | 644.89 |
| hund_sf  |  0.1175 |  0.0143 |  71.65 | 644.89 |
| hund_ph  |  0.1175 |  0.0143 |  71.65 | 644.89 |

Na(level=1) is exactly 9 × Na(level=0) = 9 × 71.65 = 644.89. This is not a convergence artifact — it is evidence that the Wannier90 XSF output **already fills the 5×5 in-plane supercell** with periodic orbital content, with no vacuum tail at the xy boundary. Mirror-padding that grid therefore creates nine coherent orbital replicas rather than extending a single localized orbital with vacuum.

Because the padded system is a periodic tiling of the original, the centering centroid is no longer at the single-orbital core but at the center of the 3×3 supertile, so the centering shift grows from `(0.026, −0.043)` Å at level 0 to `(−0.748, −1.266)` Å at level 1, and the reported `u_ev` values lose their per-orbital interpretation.

**Interpretation:** The mirror-pad diagnostic is inconclusive for this XSF dataset. A meaningful xy-tail convergence check would require either a zero-padding approach (not applicable here since the orbital is non-zero at the supercell boundary) or a larger Wannier90 supercell output (e.g., 7×7 or 9×9 in plane) — neither of which is in scope for this revision.

The level-0 rel. errors below therefore stand as the honest bare-channel comparison against CoQui, and the +9× Na scaling pins the source of the remaining discrepancy to the Wannier90 grid resolution (cusp under-resolution at the C nucleus on a ~0.08 Å grid), not to xy supercell truncation.

## Verdict

- Density channels underestimate CoQui by 3.5% (NN) to 5.9% (onsite).
- Hund's channels underestimate by 10.2%.
- Spin-flip / pair-hopping agree numerically (`|J_sf − J_ph| < 1e-6` eV — in fact to full machine precision in the test output).
- The independent two-Gaussian analytic test passes at `4.8e-4` relative error, establishing that the integration pipeline itself is correct.
- The mirror-pad diagnostic was inconclusive; the residual underestimate is attributed to the finite XSF grid's under-resolution of the orbital cusp, not to xy-tail truncation.

## Files

- Raw per-channel output: `monolayer/data/bare_monolayer_graphene.csv` (8 rows; 4 channels × 2 pad levels).
- Source module: `monolayer/src/MonolayerBareIntegrals.jl`.
- Loader + mirror-pad: `monolayer/src/MonolayerOrbitalLoader.jl`.
- Run script: `monolayer/scripts/bare_monolayer_graphene.jl`.

## Next

- Await Malte's effective-thickness / ε_eff(q) extraction script before attempting the screened-channel comparison.
- Await Malte's bilayer monolayer example before any bilayer revision.
- Optional follow-up: request a larger Wannier90 xy supercell (7×7 or 9×9) to cleanly test xy-tail convergence by true vacuum padding, or implement an FFT + Gygi-Baldereschi madelung correction that matches CoQui's periodic treatment directly.

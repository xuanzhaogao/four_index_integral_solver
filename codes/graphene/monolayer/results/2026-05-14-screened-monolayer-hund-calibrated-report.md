# Screened Monolayer Hund's at Calibrated (eps_in=3.5, LZ=3.35) — Report

Date: 2026-05-14

## Context

The 2026-05-14 eps_in sweep calibrated `(eps_in = 3.5, LZ = 3.35 Å)` as the slab-dielectric parameters that drive the screened density channels (onsite, nn) to within ~1.2 % of CoQui at `k_323201`. The Hund's channel was deferred — its previous value at `eps_in = 2.4` was 0.0902 eV vs CoQui 0.0937 eV (−3.7 %, Strong). This report measures U_hund at the new calibrated parameters.

Reference: `_coqui_crpa_loc.out` at k_323201.

| channel  | U_ijkl (eV) |
|----------|------------:|
| hund_sf  | 0.0937      |
| hund_ph  | 0.0937      |

## Method

| parameter | value |
|-----------|------:|
| LZ        | 3.35 Å (fixed) |
| Z_CENTER  | 1.675 Å (LZ/2, symmetric) |
| eps_in    | 3.5 (calibrated from density-channel sweep) |
| eps_out   | 1.0 |
| L (in-plane) | 90 Å |
| mode | Sharp only |
| source | `phi_1 · phi_2` (signed product VolumeSource) |
| solver tolerances | identical to k_323201 / eps_in-sweep runs |

Compute: Rusty `ccm`/genoa, 48 cores, Slurm job 6406269, 1m6s wall, 11.9 GB peak RSS. GMRES residual 2.87 × 10⁻⁵ (well below the 10⁻³ convergence target; ~2× higher than the eps_in=2.4 run's 1.5 × 10⁻⁵, but no instability).

## Result

| eps_in | ours `U_hund` (eV) | u_int (eV) | u_scatter (eV) | CoQui U (eV) | rel_err % |
|-------:|-------------------:|-----------:|---------------:|-------------:|----------:|
| 2.4    |            0.09021 |    0.09297 |       −0.00276 | 0.0937       |    −3.7 % |
| **3.5** |     **0.07979** |    0.08437 |       −0.00459 | 0.0937       | **−14.85 %** |

`u_total = u_int + u_scatter`. `rel_err = 100 * (ours − CoQui) / CoQui`.

## Verdict

Spec thresholds:
- Strong: `|rel_err| ≤ 5 %`
- Acceptable: `|rel_err| ≤ 15 %`
- Disagreement: otherwise

**Verdict: Acceptable, on the boundary.** At eps_in=3.5 the Hund's residual is −14.85 %, just inside the Acceptable threshold but well above the Strong band that the density channels reach at this parameter set.

## Interpretation

The pre-run expectation — that the Hund's channel would be approximately eps_in-insensitive because the `phi_1·phi_2` source has near-zero net charge — was wrong by an order of magnitude.

Decomposition of the change between eps_in=2.4 → 3.5:

- `u_int` (direct piece, screened-source volume integral) drops by 9.3 % (0.0930 → 0.0844). The "screened source" is `phi_1·phi_2 / eps_in` inside the slab; increasing eps_in from 2.4 to 3.5 directly divides this contribution by `(3.5/2.4) ≈ 1.46`, and we observe `0.0930/0.0844 ≈ 1.10`. The difference (1.10 vs 1.46) reflects the fact that a meaningful fraction of the signed-product distribution sits in the orbital overlap region — close to the slab midplane — where the rescaling is fully active; the remainder sits in the tails and is only weakly rescaled.

- `u_scatter` (interface contribution) grows by 66 % in magnitude (−0.00276 → −0.00459 eV). The Hund's-source's higher multipoles couple to the dielectric jump, and the stronger contrast at eps_in=3.5 (η = 0.556 vs 0.412 at eps_in=2.4) produces more interface charge.

Both changes pull `u_total` further away from CoQui's 0.0937 — at eps_in=2.4 we already undershot by 3.7 %, and increasing eps_in compounds the undershoot.

The qualitative consequence: **the slab model with a single scalar `eps_in` cannot simultaneously achieve Strong agreement on all three channels at this k-mesh.** The density channels prefer eps_in ~ 3.5 (Strong at +1.2 %); Hund's prefers eps_in ~ 2.4 or lower (it would over-shoot if eps_in dropped much below 2.4). Estimating from the two data points, the rate is roughly `d(rel_err)/d(eps_in) ≈ −10 %/Å⁻¹` for Hund's — so Hund's zero-rel-err would land around eps_in ≈ 2.0–2.1, where the density channels are at +20 % (Disagreement).

If we evaluate a joint verdict — max-abs-rel-err across all three channels — the trade-off looks like:

| eps_in | onsite rel_err | nn rel_err | hund rel_err (measured/est.) | joint max |
|-------:|---------------:|-----------:|------------------------------:|----------:|
| 2.4    | +20.1 %        | +18.3 %    | −3.7 % (measured)             |   20.1 %  |
| 3.5    |  +1.05 %       |  −1.18 %   | −14.85 % (measured)           |   14.85 % |

So eps_in=3.5 is actually the *better joint operating point* even with the Hund's degradation (14.85 % vs 20.1 % joint), but neither reaches the Strong band on all three channels at once.

## Reproducibility

- Driver:        `monolayer/scripts/screened_monolayer_hund_calibrated.jl`
- Submit script: `monolayer/scripts/submit_screened_monolayer_hund_calibrated.slurm`
- CSV:           `monolayer/data/screened_monolayer_hund_calibrated.csv`
- Run logs:      `monolayer/logs/hund_calibrated_6406269.{out,err}`

## Next

In order of expected payoff:

1. **Find the 3-channel joint optimum eps_in.** Run a denser eps_in sweep (e.g., 2.7, 3.0, 3.2) for the Hund's channel and overlay against the existing density-channel sweep. The joint optimum likely lies near eps_in ≈ 3.0 with max-abs-rel-err around 8–10 % (Acceptable, not Strong). One short Slurm job (~3 min, 3 Hund's solves).

2. **Reply to Malte** with the calibration story so far. Three channels, one scalar `eps_in`: density and Hund's pull in opposite directions, and the slab model with a step-function `eps(z)` has a structural limit on simultaneous fit quality. The Strong-band density result and the consistent Hund's trend are both useful inputs for the bigger conversation about what dielectric model to use.

3. **Consider an `eps(z)` profile or 2D-Lindhard model**, as a more capable replacement for the scalar slab. Justified now that we have direct evidence the scalar slab cannot exactly fit all four U_ijkl values.

4. **Defer**: k_252501 / k_161601 replication. With the calibrated `eps_in` still under negotiation, replicating before settling is premature.

The bilayer-monolayer follow-up Malte mentioned remains blocked on him producing a starting example.

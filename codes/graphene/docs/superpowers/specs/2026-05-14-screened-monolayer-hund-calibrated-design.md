# Screened Monolayer Hund's at Calibrated (eps_in=3.5, LZ=3.35) — Design

Date: 2026-05-14

## Background

The 2026-05-14 eps_in sweep at fixed `LZ = 3.35 Å` calibrated `eps_in = 3.5` as the value that drives the screened density channels to within ~1.2 % of CoQui (onsite +1.05 %, nn −1.18 %). Density channels were the only ones in the sweep — the Hund's channel was deliberately deferred, and its previous best value was at `eps_in = 2.4` from the 2026-05-13 k_323201 run: ours 0.0902 eV vs CoQui 0.0937 eV (−3.7 %, already Strong).

With the calibrated `eps_in = 3.5`, the screening environment seen by the Hund's-source `phi_1 · phi_2` distribution changes, and the resulting U_hund will shift. This report measures that shift.

Reference: `_coqui_crpa_loc.out` at k_323201.

| channel  | U_ijkl (eV) |
|----------|------------:|
| hund_sf  | 0.0937      |
| hund_ph  | 0.0937      |

## Goal

Measure `U_hund` at `(eps_in = 3.5, LZ = 3.35)` and compare against (a) CoQui's `U_hund_sf = U_hund_ph = 0.0937 eV` and (b) the previous `eps_in = 2.4` Hund's value of `0.0902` eV from the 2026-05-13 k_323201 CSV.

## Scope

In scope:

- One `BI.SharpScreening()` solve with source = target = `vs_product = phi_1 · phi_2` at the calibrated parameters.
- Output: 1-row CSV, brief dated report with the 2-row comparison table (eps_in=2.4 vs eps_in=3.5).

Out of scope:

- Soft-mix modes for Hund's.
- eps_in sweep on Hund's.
- k-mesh replication.

## Architecture

Black-box reuse of `monolayer/src/MonolayerScreenedSolve.jl` and `bilayer_slab/src/ScreenedOrbitalSolve.jl`. A new short driver script mirrors the eps_in-sweep driver but:

- Single eps_in (3.5), single channel (hund), single mode (Sharp).
- Source is `src.vs_product` (the signed product VolumeSource).
- Target spec is from `hund_target_specs(src)` (single entry: target = `src.vs_product`, Na = Nphi1, Nb = Nphi2).

**New files:**

- `monolayer/scripts/screened_monolayer_hund_calibrated.jl` — short driver (~50 lines).
- `monolayer/scripts/submit_screened_monolayer_hund_calibrated.slurm` — Slurm submit, partition `ccm`/genoa, 48 cores, 30 min wall.
- `monolayer/data/screened_monolayer_hund_calibrated.csv` — 1 data row.
- `monolayer/results/2026-05-14-screened-monolayer-hund-calibrated-report.md` — brief report.

**Reused unchanged:** `monolayer/src/MonolayerScreenedSolve.jl`, `bilayer_slab/src/ScreenedOrbitalSolve.jl`, `monolayer/src/MonolayerOrbitalLoader.jl`.

## Slab parameters

| parameter | value |
|-----------|------:|
| `LZ`      | 3.35 Å |
| `eps_in`  | 3.5 (calibrated) |
| `eps_out` | 1.0 |
| `Z_CENTER` | 1.675 Å |
| `L`       | 90 Å |
| mode      | `BI.SharpScreening()` |
| solver tolerances | reused from k_323201 / eps_in sweep |

## Data flow

```
coqui = parse_coqui_loc(_coqui_crpa_loc.out)

src = monolayer_screened_sources(orbital_1, orbital_2, source_tol=1e-3, z_center=Z_CENTER)
hund_specs = hund_target_specs(src)

result = solve_screened_mode(src.vs_product, hund_specs,
                             L, L, LZ, EPS_IN, EPS_OUT, BI.SharpScreening(); …)

for pr in result.pair_results:   # exactly one entry, pair="hund"
    row = screened_run_record(pr, result; ...) ∪ (coqui_v_ijkl, coqui_U_ijkl, diff_ev, rel_err_pct)
    push to CSV
```

## CSV schema

Same as the eps_in-sweep CSV but with 1 row and `channel = "hund"`. The `coqui_U_ijkl` is 0.0937 (using `hund_sf` from the parser — pair-hopping is degenerate to spin-flip).

## Verdict criteria

`rel_err = (ours − CoQui)/CoQui`:

- **Strong**: `|rel_err| ≤ 5 %`
- **Acceptable**: `|rel_err| ≤ 15 %`
- **Disagreement**: otherwise

Plus a delta-vs-baseline check: report the shift from the 2026-05-13 value 0.0902 to the new value. Expect order ~ a few percent shift in magnitude.

## Risks

1. **GMRES convergence on near-zero-charge source at higher eps_in.** η = (eps_in − eps_out)/(eps_in + eps_out) jumps from 0.412 (eps_in=2.4) to 0.556 (eps_in=3.5), so the BEM operator's contrast grows. The 2026-05-13 k_323201 Hund's solve at eps_in=2.4 had residual 1.5×10⁻⁵; we expect similar order at eps_in=3.5 but it could be 1–2× larger. If residual > 1e-3, flag in the report.
2. **Hund's-source dipole moment couples weakly to the long-range part of the response**, so unlike the density channels we don't expect dramatic eps_in sensitivity. A shift of ~few percent in U_hund (say 0.080–0.095) is the natural range. A surprise outside that range is worth investigating.

## Self-review

- Placeholders: none.
- Consistency: 1 eps_in × 1 mode × 1 channel = 1 CSV row.
- Scope: single sub-project, no decomposition.
- Ambiguity: "calibrated" = `eps_in = 3.5, LZ = 3.35` from the 2026-05-14 eps_in-sweep verdict. Explicit.

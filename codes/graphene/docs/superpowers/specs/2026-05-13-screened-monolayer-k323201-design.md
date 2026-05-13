# Screened Monolayer Graphene at k_323201 — Design

Date: 2026-05-13

## Background

The 2026-05-13 bare-monolayer k-mesh sweep showed our direct Wannier integral and CoQui's bare V agree across all 4 channels to ~10 mV at the densest k-mesh (`k_323201`). With that joint reference established, the natural next step is the screened-interaction comparison: reproduce CoQui's `U_ijkl` (cRPA-screened) at the same dataset using a slab dielectric model.

Reference U values from `_coqui_crpa_loc.out` in the k_323201 directory:

| i j k l | physical channel | v_ijkl (bare, eV) | U_ijkl (screened, eV) |
|---------|------------------|------------------:|----------------------:|
| 0 0 0 0 | onsite           | 17.4029           | 9.7824                |
| 0 0 1 1 | nn               |  8.8319           | 5.1678                |
| 0 1 0 1 | Hund's pair-hop  |  0.1316           | 0.0937                |
| 0 1 1 0 | Hund's spin-flip |  0.1316           | 0.0937                |

(The other 12 entries in the file are permutation-symmetric duplicates.)

## Goal

Run the existing dielectric-slab GMRES solver on the k_323201 monolayer Wannier orbitals using user-specified slab parameters (`LZ = 3.35` Å, `eps_in = 2.4`, `eps_out = 1.0`), evaluate all 4 distinct screened channels, and produce a CSV + report comparing against the CoQui U values above.

## Scope

In scope:

- 2 GMRES dielectric solves (density-source for onsite+nn; Hund's-source for spin-flip+pair-hop).
- 4 screened channel values produced and tabulated.
- Sweep over `Sharp` + 4 soft-screening bandwidths `[0.05, 0.1, 0.2, 0.4]` for the density solve (reusing the bilayer's protocol). For the Hund's solve we run only `Sharp` first (cheaper; soft modes can be added later if needed).
- Parse U_ijkl values for the 4 channels from `_coqui_crpa_loc.out`.
- Dated markdown report with verdict.

Out of scope (deferred):

- k-mesh sweep on screened (only k_323201).
- Sweep over `eps_in` or `LZ`.
- Bilayer monolayer comparison.
- Soft modes for Hund's channel (deferred to a follow-up if Hund's `Sharp` result is off).
- Refactoring or forking `bilayer_slab/src/ScreenedOrbitalSolve.jl` (consumed as a black-box dependency).

## Architecture

Black-box reuse of `bilayer_slab/src/ScreenedOrbitalSolve.jl`. The new code is a thin monolayer wrapper that builds the 4-channel sources/targets and calls `solve_screened_mode` twice.

**New files:**

- `monolayer/src/MonolayerScreenedSolve.jl` — wrapper module:
  - `include`s `bilayer_slab/src/ScreenedOrbitalSolve.jl` and uses `solve_screened_mode` as the workhorse.
  - Re-uses `MonolayerOrbitalLoader.centered_monolayer_sources_padded` for VolumeSource construction (sets `mirror_pad_level = 0` and applies a z-center shift so the orbital sits at slab midplane).
  - Constructs the `phi_0 * phi_1` signed-product VolumeSource for the Hund's channel using the same `_product_datagrid` helper logic as `MonolayerBareIntegrals`.
  - Builds two `pair_specs`-shaped lists (density-channel list of 2 entries; Hund's-channel list of 1 entry) and passes each to `solve_screened_mode`.
  - Provides `parse_coqui_loc(path)` that reads the U_ijkl table and returns a Dict of 4 channel → (v_ijkl, U_ijkl) tuples.
- `monolayer/scripts/screened_monolayer_k323201.jl` — top-level driver: loads sources, runs the 2 solves over the bandwidth sweep (density-channel only) plus a single Sharp Hund's solve, parses CoQui, writes CSV + summary print.
- `monolayer/data/screened_monolayer_k323201.csv` — produced by the script.
- `monolayer/results/2026-05-13-screened-monolayer-k323201-report.md` — dated report.

**Reused unchanged:**

- `bilayer_slab/src/ScreenedOrbitalSolve.jl` (consumed via include).
- `monolayer/src/MonolayerOrbitalLoader.jl` (for orbital loading + product datagrid building).
- `monolayer/src/MonolayerBareIntegrals.jl` (the bare-channel values from k_323201 are pulled directly from `monolayer/data/bare_monolayer_kmesh_sweep.csv`, no recomputation).

## Slab parameters

| parameter | value | rationale |
|-----------|------:|-----------|
| `LZ`         | 3.35 Å    | user-specified monolayer thickness |
| `EPS_IN`     | 2.4       | user-specified |
| `EPS_OUT`    | 1.0       | vacuum |
| `Z_CENTER`   | LZ/2 = 1.675 Å | symmetric placement: orbital centered between top and bottom interfaces |
| `L`          | 90 Å      | reused from bilayer; well beyond orbital extent |
| `SOURCE_TOL` | 1e-3      | reused |
| `N_QUAD`     | 6         | reused |
| `EDGE_REFINE_LEVEL` | 4   | reused |
| `RHS_TOL`    | 1e-3      | reused |
| `LHS_TOL`    | 1e-5      | reused |
| `GMRES_ATOL`/`RTOL` | 1e-5 | reused |
| `MAX_ORDER`  | 64        | reused |
| `MAX_DEPTH`  | 12        | reused |
| `BANDWIDTHS` | `[0.05, 0.1, 0.2, 0.4]` | reused (density-channel sweep) |
| Hund's modes | `Sharp` only | first pass; soft modes deferred |

## Data flow

```
load k_323201 squared orbitals (vs1, vs2)  via centered_monolayer_sources_padded
load k_323201 signed orbitals (datagrid_1, datagrid_2)
build product datagrid phi1*phi2  →  vs_product

# Density solve: source = vs1 (|phi1|^2)
for case in {Sharp, soft(0.05), soft(0.1), soft(0.2), soft(0.4)}:
    solve_screened_mode(vs1, [target=vs1 (:onsite), target=vs2 (:nn)],
                        L, L, LZ, EPS_IN, EPS_OUT, case.mode; ...)
    → 2 channel rows × this case

# Hund's solve: source = vs_product, target = vs_product
for case in {Sharp}:
    solve_screened_mode(vs_product, [target=vs_product (:hund_sf)],
                        L, L, LZ, EPS_IN, EPS_OUT, case.mode; ...)
    → 1 channel row × this case

parse CoQui U_ijkl from _coqui_crpa_loc.out  →  4 reference values

assemble CSV: (channel, mode, bandwidth, our_u_ev, coqui_u_ev, diff_ev, ...solver diagnostics)
```

For density channels the CSV has 2 channels × 5 modes = 10 rows. For Hund's, 1 channel × 1 mode = 1 row. We label the spin-flip / pair-hopping degeneracy with one `:hund` row labeled as both `hund_sf` and `hund_ph` (or duplicate the row for clarity in the CSV — TBD at implementation time, default: single row labeled `hund`).

## Hund's source treatment

For the screened solver the source charge density is `phi_1(r) * phi_2(r)` (signed product, integrates to ~0 over space because phi_1 and phi_2 are orthogonal Wannier orbitals). This is conceptually a "dipole-like" source — total charge ~0 but with non-trivial multipole content. The Hund's matrix element is the screened potential generated by this source, evaluated at the same source distribution.

Normalization for the eV conversion follows the bare-channel convention from `MonolayerBareIntegrals.bare_channel_integral`: divide by `(Nphi1 * Nphi2)` where `Nphi_i = integral |phi_i|^2 dr`, NOT by the product integrals (~0). The `pair_specs`-shaped entry uses `Na = Nphi1`, `Nb = Nphi2`.

Risk: the bilayer solver was developed for density (positive) sources. Whether GMRES converges cleanly on a near-zero-net-charge source is an open question for first run. If GMRES residual or potential are noisy, fall back to Sharp-only and note in the report.

## CoQui parser

Lines of interest in `_coqui_crpa_loc.out`:

```
  0  0  0  0   17.4029    9.7824
  0  0  1  1    8.8319    5.1678
  0  1  0  1    0.1316    0.0937
  0  1  1  0    0.1316    0.0937
```

Regex: `^\s*(\d)\s+(\d)\s+(\d)\s+(\d)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)`. Map (i,j,k,l) → channel:

- (0,0,0,0) → `:onsite`
- (0,0,1,1) → `:nn`
- (0,1,0,1) → `:hund_ph`
- (0,1,1,0) → `:hund_sf`

Returns `Dict{Symbol, NamedTuple{(:v_ijkl, :U_ijkl), Tuple{Float64, Float64}}}` with 4 entries.

## Verdict criteria

For each channel, define `rel_err = (our_u_ev - coqui_U_ijkl) / coqui_U_ijkl`. The report's verdict:

- **Strong agreement**: `|rel_err| ≤ 5%` for onsite and nn.
- **Acceptable**: `|rel_err| ≤ 15%` for onsite and nn.
- **Disagreement**: otherwise.

Hund's is informational only (smallest values, largest relative numerical noise expected).

## Reproducibility / outputs

- Driver script: `monolayer/scripts/screened_monolayer_k323201.jl`.
- Raw CSV:       `monolayer/data/screened_monolayer_k323201.csv`.
- Report:        `monolayer/results/2026-05-13-screened-monolayer-k323201-report.md`.
- Run log:       `monolayer/logs/screened_k323201_*.log` (committed selectively if interesting; otherwise gitignored — bilayer pattern doesn't commit logs).

## Risks / loose ends

1. **Z-placement convention.** The bilayer used `Z_CENTER = LZ/4` (asymmetric — single layer of a bilayer at quarter-height). The monolayer placement `LZ/2` is the natural symmetric choice but should be verified by a quick sanity output: the source-VolumeSource z-coordinates should be centered around `1.675 Å`. Easy first-output check.

2. **Slab thickness 3.35 Å.** This is the standard graphite interlayer spacing — a sensible monolayer "effective thickness" but not a derived constant. If U values disagree systematically with CoQui, this is a candidate parameter sweep. Not in scope for this first pass.

3. **eps_in = 2.4.** Same as the bilayer value. The user has committed to it for this first run. Sweep deferred.

4. **Hund's GMRES convergence on near-zero-charge source.** First-pass risk. Mitigation: if Hund's `Sharp` mode produces nonsense or fails to converge within `GMRES_ATOL`, drop Hund's from the CSV with a clear `missing` and flag in the report.

5. **Compute time.** Each `solve_screened_mode` invocation built 530k+ source points on the 2026-04-23 / 2026-05-13 bare runs. The screened solve adds a GMRES iteration loop. Wall-clock estimate: each mode invocation could be 5-30 min. With 5 density modes + 1 Hund's mode = 6 solves total, the script could run 30 min - 3 hours. Plan tasks accordingly (background julia, log polling — as done for the bare sweep in Task 4).

## Self-review

- **Placeholders.** One intentional "TBD at implementation time" — whether the Hund's row in the CSV is single or duplicated for sf/ph labeling. Default chosen (single row labeled `hund`). Acceptable.
- **Internal consistency.** Two GMRES solves, 5 + 1 mode cases, 4 channels, ~11 CSV rows. The Hund's-source treatment and normalization convention match what `MonolayerBareIntegrals` already does for the bare-channel Hund's. The slab geometry (`Z_CENTER = LZ/2`) matches the symmetric-monolayer rationale.
- **Scope.** One sub-project: reproduce 4 distinct U values at one (kmesh, slab-geometry) point. No decomposition needed.
- **Ambiguity.** "Reproduce" means within the verdict thresholds (5% / 15%). No looser meaning intended.

# Bare Monolayer Graphene — k-mesh Sweep Design

Date: 2026-05-13

## Background

The 2026-04-23 bare-monolayer report compared our direct `1/|r-r'|` integral of the Wannier orbitals against Malte Rösner's CoQui cRPA bare reference at k-mesh `16×16×1`:

- our `v0` (onsite) = 16.400 eV
- CoQui `v0`       = 17.434 eV
- discrepancy      = 1.034 eV

In his 2026-05 reply, Malte identified that this discrepancy is suspiciously close to CoQui's "long-wavelength Madelung correction" for that k-mesh (1.041 eV). To test the hypothesis that

  `CoQui_v0 = direct_Wannier_v0 + Madelung_correction(k-mesh)`

he provided two new datasets at different k/q samplings, whose Madelung corrections differ in sign and magnitude:

| k-mesh     | reported Madelung correction |
|------------|-----------------------------:|
| k_161601   |  +1.041 eV  (existing)       |
| k_252501   |  +0.1842 eV (new)            |
| k_323201   |  −0.290 eV  (new)            |

If the hypothesis holds, then for each dataset `CoQui_v0 − our_v0` should equal the reported Madelung correction within a tight tolerance.

Data location:
- `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_161601_nb_144_c_15`
- `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_252501_nb_144_c_15`
- `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15`

Each directory contains `graphene_00001.xsf`, `graphene_00002.xsf`, and a CoQui output file (the new ones use `_coqui_crpa_loc.out`; the original used `_coqui_thc_crpa.out`).

## Goal

Run the existing bare 4-channel pipeline on the two new k-mesh datasets and produce a single artifact that lets us (and Malte) verify the Madelung-correction hypothesis directly.

## Scope

In scope:

- Run `compute_all_bare_channels` for all three k-meshes × four channels (`onsite`, `nn`, `hund_sf`, `hund_ph`) at `mirror_pad_level = 0`.
- Parse CoQui bare-interaction reference values from each dataset's `_coqui_*_crpa*.out` file.
- Write a combined CSV with our values, CoQui values, their difference, and the Madelung-correction reference.
- Write a markdown report that states whether the Madelung-correction hypothesis is confirmed by the residuals.

Out of scope (deferred or owned elsewhere):

- Mirror-pad sweep on the new datasets (the 2026-04-23 report already concluded that pad-1 is inconclusive because the XSF tiles the supercell).
- Implementing our own Madelung correction in this codebase.
- Reply to Malte about regenerating XSF at higher real-space resolution (3×3×1 supercell): decided after seeing the residuals — if they're small the cusp-resolution question is moot for `v0`.
- Re-running existing 2026-04-23 artifacts; that report stands as the original snapshot.

## Architecture

One new self-contained script that drives the existing pipeline. No changes to `MonolayerOrbitalLoader` or `MonolayerBareIntegrals` — they already accept arbitrary XSF input.

**New files:**

- `monolayer/scripts/bare_monolayer_kmesh_sweep.jl` — loops over three dataset directories, invokes `compute_all_bare_channels`, parses CoQui output, writes the combined CSV.
- `monolayer/data/bare_monolayer_kmesh_sweep.csv` — 12 rows (3 k-meshes × 4 channels).
- `monolayer/results/2026-05-13-bare-monolayer-kmesh-sweep.md` — narrative + Madelung-hypothesis verdict.

**Reused unchanged:**

- `MonolayerOrbitalLoader.centered_monolayer_sources_padded` with `mirror_pad_level = 0`.
- `MonolayerBareIntegrals.compute_all_bare_channels`.

## Data flow

For each `kmesh in ("k_161601", "k_252501", "k_323201")`:

1. Build paths `orbital_1`, `orbital_2`, and the CoQui output file.
2. Call `compute_all_bare_channels(orbital_1, orbital_2, source_tol=1e-3, volume_tol=1e-3, mirror_pad_level=0)`.
3. Parse the CoQui bare-interaction lines for `:onsite`, `:nn`, `:hund_sf`, `:hund_ph` from the output file.
4. Look up the Madelung correction for this k-mesh from a small hard-coded table populated from Malte's email.
5. Emit four rows: `(kmesh, channel, our_u_ev, coqui_u_ev, diff_ev = coqui − ours, madelung_ref_ev, madelung_residual_ev = diff_ev − madelung_ref_ev, ...metadata)`.
6. `madelung_ref_ev` and `madelung_residual_ev` are populated only for `:onsite` (the channel for which the Madelung correction is defined); they are `missing` for the other three channels.

The existing `k_161601` row is recomputed inside this script so all three rows go through identical code on the same day, rather than mixing in a stored value from the 2026-04-23 CSV.

## CSV schema

Columns of `monolayer/data/bare_monolayer_kmesh_sweep.csv`:

| column | meaning |
|---|---|
| `kmesh` | string tag, e.g. `"k_161601"` |
| `channel` | one of `"onsite"`, `"nn"`, `"hund_sf"`, `"hund_ph"` |
| `our_u_ev` | result of `bare_channel_integral` (this work, eV) |
| `coqui_u_ev` | parsed CoQui bare value (eV); `missing` if parsing fails |
| `diff_ev` | `coqui_u_ev − our_u_ev` |
| `madelung_ref_ev` | Malte-supplied Madelung correction; populated only for `onsite` |
| `madelung_residual_ev` | `diff_ev − madelung_ref_ev`; only for `onsite` |
| `our_u_raw` | raw integral value before eV conversion |
| `Na`, `Nb` | source and target orbital norms used by `bare_channel_integral` |
| `n_source_points`, `n_target_points` | quadrature sizes |
| `tkm_kmax` | TKM3D `kmax` used |
| `source_tol`, `volume_tol` | tolerances used |
| `mirror_pad_level` | always `0` |
| `shift_x`, `shift_y`, `shift_z` | shared centering shift applied to both orbitals |

## CoQui output parsing

The original dataset's CoQui file is `_coqui_thc_crpa.out`; the two new datasets' file is `_coqui_crpa_loc.out`. The plan accepts both names — try each in order and use the first that exists.

Lines of interest (form determined by inspecting one of the new files at implementation time):

- intra-orbital / onsite "v0"
- inter-orbital / nearest-neighbor
- Hund's spin-flip
- Hund's pair-hopping

If line formats differ between `_coqui_thc_crpa.out` and `_coqui_crpa_loc.out`, the parser handles both forms (regex over a small whitelist of label patterns). If a channel can't be parsed, emit `missing` and continue — the Madelung test only needs the onsite line.

## Hypothesis check (in the report)

The Madelung-hypothesis verdict is the central artifact. For each k-mesh the report shows:

| k-mesh   | our `v0` (eV) | CoQui `v0` (eV) | `CoQui − ours` (eV) | Madelung (Malte, eV) | residual (eV) |
|----------|--------------:|----------------:|--------------------:|---------------------:|--------------:|
| k_161601 |          16.40 |          17.434 |               1.034 |                1.041 |        −0.007 |
| k_252501 |              ? |               ? |                   ? |               0.184  |             ? |
| k_323201 |              ? |               ? |                   ? |              −0.290  |             ? |

Verdict criteria:

- **Confirmed** if `|residual| ≤ 0.05 eV` for all three rows.
- **Partial** if `|residual| ≤ 0.1 eV` for all three rows.
- **Rejected** otherwise; in which case the report flags this for follow-up (XSF resolution, deeper Madelung formulation, etc.) without choosing a remediation here.

The Hund's and NN channels are reported in the same CSV but the report frames them as "informational; the Madelung correction is a long-wavelength onsite artifact and we don't expect it to apply to these short-range channels."

## Testing

No new automated tests. The sweep script reuses `compute_all_bare_channels`, which is already covered by the 2026-04-23 test suite. The script itself is a thin orchestration layer + a parser; manual inspection of the generated CSV is the verification.

If the CoQui parser is non-trivial enough to merit one, we may add a small parser test against a fixture string. Decided at implementation time, not committed to here.

## Risks / loose ends

- **CoQui parser** — the new files use a different naming convention from the original, and their line layout may differ. The plan inspects one new file at implementation time and adapts the parser. If parsing fails on a non-onsite channel, that row's `coqui_u_ev` is `missing` and we move on.
- **Lattice / supercell** — the new XSFs may use a different supercell than the original 5×5×1. The loader is generic and should handle any lattice, but we sanity-check the centering shift for each dataset (orbital-1 centroid should land at one of the carbon atoms of the unit cell, and orbital-2 should be shifted by the A→B bond vector).
- **Re-running k_161601** — recomputing inside this script may yield a value slightly different from the 16.4004 stored in the 2026-04-23 CSV (different machine state, BLAS thread count, etc.). Any drift > 1e-3 eV is itself worth a note in the report.

## Self-review

- Placeholders: the verdict-criteria thresholds (`0.05`, `0.1` eV) and the table cells marked `?` are intentional — the residuals are the experiment output. No `TBD`s remain in the design.
- Internal consistency: the data-flow, CSV schema, and report match. The Madelung-test framing is restricted to `:onsite` consistently throughout.
- Scope: a single self-contained sweep + report. No decomposition needed.
- Ambiguity: only the CoQui parser pattern is left until file inspection — explicit in the "Risks" section.

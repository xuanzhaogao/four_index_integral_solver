# Screened Monolayer Graphene — LZ Sweep at fixed eps_in=2.4 (Design)

Date: 2026-05-13

## Background

The 2026-05-13 screened-monolayer report at `k_323201` with `LZ = 3.35 Å, eps_in = 2.4` produced:

- onsite: 11.7471 eV vs CoQui 9.7824 eV (+20.1%)
- nn: 6.1130 eV vs CoQui 5.1678 eV (+18.3%)
- hund: 0.0902 eV vs CoQui 0.0937 eV (−3.7%)

The Hund's channel agreement (~4%) on a near-zero-charge dipole source under the same slab parameters rules out pipeline-level errors. The ~20% under-screening on the density channels is a calibration question: either `eps_in = 2.4` is too low, or the slab is too thin (orbital tails leak out into vacuum where they see no screening), or both.

The user has chosen to fix `eps_in = 2.4` and sweep `LZ` first — physically, `LZ` parameterizes how much of the orbital sits inside the screened region vs the unscreened vacuum, and a wider sweep is the natural way to find a calibration point.

## Goal

Find the slab thickness `LZ` that drives the density channels (onsite + nn) within 5% of CoQui's `U_ijkl` reference at the joint validation point `k_323201`, at fixed `eps_in = 2.4`.

## Scope

In scope:

- Density-source `solve_screened_mode` at `BI.SharpScreening()` for each of `LZ ∈ {2.0, 3.0, 3.35, 5.0, 8.0}` Å.
- Each solve evaluates two channels (`:onsite`, `:nn`).
- CoQui `_coqui_crpa_loc.out` is parsed once for the comparison values.
- Output: 10-row CSV + dated markdown report with the LZ-vs-rel_err curve and a verdict on the best LZ.

Out of scope (deferred):

- SoftMix bandwidth sweep — the k_323201 sweep showed they only worsen agreement on density channels.
- Hund's channel — already in Strong band (−3.7%), not the calibration target. Including it would roughly double the compute.
- `eps_in` sweep — held fixed at 2.4 for this round per the user's directive.
- Asymmetric Z_CENTER (orbital not at slab midplane).
- k_252501 / k_161601 replication.
- Bilayer monolayer comparison.

## Architecture

Black-box reuse of `monolayer/src/MonolayerScreenedSolve.jl` (parser, source builder, target_specs builder, run_record flattener) and the bilayer's `solve_screened_mode`. A new short driver script loops over LZ values.

**New files:**

- `monolayer/scripts/screened_monolayer_lz_sweep.jl` — driver. ~80 lines. Loops over `LZ_VALUES`, rebuilds the centered sources at `Z_CENTER = LZ/2` per LZ, runs one Sharp density solve, parses CoQui once, writes CSV.
- `monolayer/scripts/submit_screened_monolayer_lz_sweep.slurm` — Slurm submit script for Rusty `ccm` / genoa, 48 cores, 2-hour wall (generous margin over the ~15-min estimate).
- `monolayer/data/screened_monolayer_lz_sweep.csv` — 10 rows (5 LZ × 2 channels).
- `monolayer/results/2026-05-13-screened-monolayer-lz-sweep-report.md` — dated report.

**Reused unchanged:**

- `monolayer/src/MonolayerScreenedSolve.jl` — all four exported helpers consumed as-is. `monolayer_screened_sources` already accepts `z_center` as a kw arg, so the per-LZ source reconstruction is a single call.
- `bilayer_slab/src/ScreenedOrbitalSolve.jl` — `solve_screened_mode` consumed as a black box.
- `monolayer/src/MonolayerOrbitalLoader.jl` — orbital + product datagrid loaders.

## Slab parameters

| parameter | value | rationale |
|-----------|------:|-----------|
| `eps_in`  | 2.4   | user-fixed for this round |
| `eps_out` | 1.0   | vacuum |
| `LZ`      | {2.0, 3.0, 3.35, 5.0, 8.0} Å | five anchors bracketing 3.35 on both sides |
| `Z_CENTER` | `LZ/2` (per LZ) | symmetric: orbital at slab midplane |
| `L` (in-plane box) | 90 Å | reused; well beyond orbital extent |
| `source_tol` | 1e-3 | reused |
| `mode` | `BI.SharpScreening()` only | soft modes shown to worsen agreement |
| `N_QUAD`, `EDGE_REFINE_LEVEL`, `RHS_TOL`, `LHS_TOL`, `GMRES_ATOL`, `GMRES_RTOL`, `MAX_ORDER`, `MAX_DEPTH` | reused from k_323201 | proven values |

## Data flow

```
coqui = parse_coqui_loc(_coqui_crpa_loc.out path for k_323201)

for LZ in (2.0, 3.0, 3.35, 5.0, 8.0):
    z_center = LZ / 2
    src = monolayer_screened_sources(orbital_1, orbital_2, source_tol=1e-3, z_center=z_center)
    density_specs = density_target_specs(src)
    result = solve_screened_mode(src.vs1, density_specs,
                                 L, L, LZ, EPS_IN, EPS_OUT, BI.SharpScreening();
                                 …solver tolerances…)
    for pr in result.pair_results:                       # pr.pair ∈ {"onsite","nn"}
        emit row (LZ, channel, u_int_ev, u_scatter_ev, u_total_ev,
                  coqui_U_ijkl, diff_ev, rel_err_pct,
                  sigma_residual, n_interface_points, Nphi1, Nphi2, ...)
```

`coqui_U_ijkl` is looked up from the parsed Dict keyed by `:onsite` or `:nn`. `rel_err_pct = 100 * (u_total_ev − coqui_U_ijkl) / coqui_U_ijkl`.

## CSV schema

Columns of `monolayer/data/screened_monolayer_lz_sweep.csv` (10 rows):

| column | meaning |
|---|---|
| `LZ` | slab thickness (Å) |
| `channel` | `"onsite"` or `"nn"` |
| `u_int_ev` | direct part (screened source volume integral, eV) |
| `u_scatter_ev` | interface-scattering part (eV) |
| `u_total_ev` | `u_int_ev + u_scatter_ev` (eV) |
| `coqui_U_ijkl` | reference (constant per channel) |
| `diff_ev` | `coqui_U_ijkl − u_total_ev` |
| `rel_err_pct` | `100 * (u_total_ev − coqui_U_ijkl) / coqui_U_ijkl` |
| `sigma_residual` | GMRES final residual |
| `n_interface_points` | adaptive interface discretization size |
| `Nphi1`, `Nphi2` | orbital-norm integrals (constant per LZ; vary slightly with z_center via the adaptive source tol but should be ~78) |
| `eps_in`, `eps_out`, `Z_CENTER`, `L`, `source_tol`, `rhs_tol`, `lhs_tol`, `gmres_atol`, `gmres_rtol`, `n_quad`, `edge_refine_level`, `max_order`, `max_depth` | run metadata for reproducibility |

## Report structure

- Context: link back to k_323201 report and the user's choice to sweep LZ first.
- Method: 5-row parameter table.
- Results: 10-row CSV table; a "rel_err vs LZ" comparison table:

  | LZ (Å) | onsite ours | onsite rel_err | nn ours | nn rel_err |
  |--------|------------:|---------------:|--------:|-----------:|
  | 2.0    | <fill>      | <fill>         | <fill>  | <fill>     |
  | 3.0    | <fill>      | <fill>         | <fill>  | <fill>     |
  | 3.35   | <fill>      | <fill>         | <fill>  | <fill>     |
  | 5.0    | <fill>      | <fill>         | <fill>  | <fill>     |
  | 8.0    | <fill>      | <fill>         | <fill>  | <fill>     |

- Verdict (Strong/Acceptable/Disagreement on the best LZ).
- Interpretation (does a single LZ work for both channels? are the trends monotonic?).
- Next: depending on the verdict — either declare LZ calibrated and replicate at other k-meshes, or open the eps_in sweep, or recommend a non-step `eps(z)` model.

## Verdict criteria

For each LZ, compute `max(|rel_err_onsite|, |rel_err_nn|)`. The "best" LZ is the one minimizing this maximum.

- **Strong**: best LZ has `max ≤ 5%`. Declare LZ calibrated.
- **Acceptable**: best LZ has `max ≤ 15%`. Useful step forward but not done; consider eps_in sweep.
- **Disagreement**: best LZ has `max > 15%`. Either LZ alone cannot close the gap, or the optimum is outside [2, 8] Å — recommend widening the LZ range or moving to a 2D-Lindhard model.

## Risks / loose ends

1. **Small LZ (2.0, 3.0 Å) with the orbital tail extending well outside.** The orbital z-extent (from the XSF) is multiple Å; at LZ=2.0 most of the density may sit outside the slab. Pipeline should handle this — `BI.screened_volume_source` divides density by `eps_in` only inside the slab, outside it's untouched. The result is just a weaker effective screening. If GMRES has trouble converging at small LZ (residual > 1e-3), drop those rows and note in the report.
2. **Adaptive interface discretization scaling with LZ.** Larger LZ means a larger physical box and potentially more adaptive panels (and so more GMRES cost). The k_323201 run had `n_interface_points = 960552`; if LZ=8 grows this by 2-3×, wall-clock per solve goes up correspondingly. The 2-hour Slurm allocation has plenty of margin.
3. **`Z_CENTER = LZ/2` for every LZ.** Symmetric. Confirmed with the user; out-of-scope alternative is asymmetric placement.
4. **Bare `u_int_ev` will look different at each LZ.** The screened volume source's effective scaling depends on what fraction of the orbital is inside the slab. Worth tabulating in the report — it's part of the diagnostic.

## Self-review

- **Placeholders.** `<fill>` markers only appear in the report template (Task 7 of the plan-to-come will populate). No `TBD` elsewhere.
- **Internal consistency.** Single mode × 5 LZ × 2 channels = 10 CSV rows; matches the schema and report's table.
- **Scope.** One sub-project: parameter sweep at fixed eps_in. No decomposition needed.
- **Ambiguity.** "Find the LZ" means within the verdict thresholds (5% Strong / 15% Acceptable); explicit.

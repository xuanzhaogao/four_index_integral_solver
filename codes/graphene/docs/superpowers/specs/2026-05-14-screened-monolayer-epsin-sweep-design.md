# Screened Monolayer Graphene — eps_in Sweep at fixed LZ=3.35 (Design)

Date: 2026-05-14

## Background

The 2026-05-13 LZ sweep at fixed `eps_in = 2.4` showed a strictly monotonic improvement of agreement with CoQui as `LZ` grows, with no interior optimum in [2.0, 8.0] Å. At the best swept value LZ=8 the residual was still onsite +15.7%, nn +10.2% — i.e., on the boundary of the Acceptable verdict band. The report concluded that the slab-thickness lever alone cannot reach the Strong (≤5%) band within physically reasonable LZ; the next physical lever is the in-slab permittivity `eps_in`.

The user has chosen to return to the originally-motivated LZ=3.35 Å (graphite interlayer spacing) and sweep `eps_in` over [2.0, 3.5] instead.

## Goal

Find the `eps_in ∈ [2.0, 3.5]` at fixed `LZ = 3.35 Å` that drives both density channels (onsite + nn) within 5% of CoQui's `U_ijkl` reference at `k_323201`.

## Scope

In scope:

- Density-source `solve_screened_mode` at `BI.SharpScreening()` for each of `eps_in ∈ {2.0, 2.4, 2.7, 3.0, 3.2, 3.5}`.
- Each solve evaluates two channels (`:onsite`, `:nn`).
- Cross-check: the `eps_in = 2.4` row must reproduce the 2026-05-13 k_323201 baseline (onsite 11.747, nn 6.113) to 4 decimals. Failure here means pipeline regression.
- Output: 12-row CSV + dated markdown report with the eps_in-vs-rel_err curve and a verdict on the best eps_in.

Out of scope (deferred):

- SoftMix bandwidth sweep — the k_323201 work showed soft modes only worsen density-channel agreement.
- Hund's channel — already in Strong band (−3.7%) at the baseline.
- `LZ` further tuning — held fixed at 3.35 for this round per the user's directive.
- k_252501 / k_161601 replication.
- Widening the eps_in range outside [2.0, 3.5] (decided after seeing this sweep's results).

## Architecture

Black-box reuse of `monolayer/src/MonolayerScreenedSolve.jl` (parser, source builder, target_specs, run_record flattener) and the bilayer's `solve_screened_mode`. A new short driver script loops over eps_in values. Geometry is fixed across all rows: `Z_CENTER = LZ/2 = 1.675 Å`, `L = 90 Å`.

**New files:**

- `monolayer/scripts/screened_monolayer_epsin_sweep.jl` — driver. ~80 lines. Loops over `EPS_IN_VALUES`, runs one Sharp density solve per value, parses CoQui once, writes CSV.
- `monolayer/scripts/submit_screened_monolayer_epsin_sweep.slurm` — Slurm submit script. Mirror of `submit_screened_monolayer_lz_sweep.slurm` with a 1-hour wall (plenty for ~8 min compute).
- `monolayer/data/screened_monolayer_epsin_sweep.csv` — 12 rows (6 eps_in × 2 channels). Produced by the Slurm job.
- `monolayer/results/2026-05-14-screened-monolayer-epsin-sweep-report.md` — dated report.

**Reused unchanged:**

- `monolayer/src/MonolayerScreenedSolve.jl` — all four exported helpers consumed as-is.
- `bilayer_slab/src/ScreenedOrbitalSolve.jl` — `solve_screened_mode` consumed as a black box.
- `monolayer/src/MonolayerOrbitalLoader.jl` — orbital + product datagrid loaders.

## Slab parameters

| parameter | value | rationale |
|-----------|------:|-----------|
| `LZ`      | 3.35 Å (fixed) | user-fixed; graphite interlayer spacing |
| `Z_CENTER` | 1.675 Å (LZ/2) | symmetric placement |
| `eps_out` | 1.0 | vacuum |
| `eps_in`  | {2.0, 2.4, 2.7, 3.0, 3.2, 3.5} | six anchors covering the requested range |
| `L` (in-plane box) | 90 Å | reused |
| `mode` | `BI.SharpScreening()` only | soft modes proven to worsen agreement |
| solver tolerances | reused from k_323201 / LZ sweep | proven values |

## Data flow

```
coqui = parse_coqui_loc(_coqui_crpa_loc.out path for k_323201)

# sources only need to be built once — geometry is fixed across the sweep
src = monolayer_screened_sources(orbital_1, orbital_2, source_tol=1e-3, z_center=LZ/2)

for eps_in in (2.0, 2.4, 2.7, 3.0, 3.2, 3.5):
    density_specs = density_target_specs(src)
    result = solve_screened_mode(src.vs1, density_specs,
                                 L, L, LZ, eps_in, EPS_OUT, BI.SharpScreening();
                                 …solver tolerances…)
    for pr in result.pair_results:
        emit row (eps_in, channel, u_int_ev, u_scatter_ev, u_total_ev,
                  coqui_U_ijkl, diff_ev, rel_err_pct, ...solver diagnostics)
```

The sources do NOT need to be re-built per eps_in (unlike the LZ sweep, where Z_CENTER varies with LZ). This shaves the loader cost off every iteration after the first.

`rel_err_pct = 100 * (u_total_ev − coqui_U_ijkl) / coqui_U_ijkl`.

## CSV schema

Columns of `monolayer/data/screened_monolayer_epsin_sweep.csv` (12 rows):

| column | meaning |
|---|---|
| `eps_in` | swept in-slab permittivity |
| `channel` | `"onsite"` or `"nn"` |
| `u_int_ev` | direct part (screened-source volume integral, eV) |
| `u_scatter_ev` | interface-scattering part (eV) |
| `u_total_ev` | `u_int_ev + u_scatter_ev` (eV) |
| `coqui_U_ijkl` | reference (constant per channel) |
| `diff_ev` | `coqui_U_ijkl − u_total_ev` |
| `rel_err_pct` | `100 * (u_total_ev − coqui_U_ijkl) / coqui_U_ijkl` |
| `sigma_residual` | GMRES final residual |
| `n_interface_points` | adaptive interface discretization size |
| `Nphi1`, `Nphi2` | orbital norms (constant across rows — same source for every row) |
| `Lz`, `Z_CENTER`, `L`, `eps_out`, `source_tol`, `rhs_tol`, `lhs_tol`, `gmres_atol`, `gmres_rtol`, `n_quad`, `edge_refine_level`, `max_order`, `max_depth` | run metadata for reproducibility |

## Report structure

- Context: links to the LZ-sweep report (the immediate predecessor) and to the k_323201 baseline.
- Method: parameter table.
- Results: 12-row CSV + the eps_in-vs-rel_err table.
- Cross-check: the `eps_in = 2.4` row reproduces the k_323201 baseline (must match to ~4 decimals).
- Verdict (Strong / Acceptable / Disagreement on the best eps_in).
- Interpretation: are the curves monotonic? does a single eps_in work for both channels?
- Next: depending on verdict — declare eps_in calibrated and move to k-mesh replication, widen the range, or open a 2D-Lindhard model.

## Verdict criteria

For each `eps_in`, compute `max(|rel_err_onsite|, |rel_err_nn|)`. The "best" eps_in is the one minimizing this maximum.

- **Strong**: best eps_in has `max ≤ 5%`. Declare eps_in calibrated.
- **Acceptable**: best eps_in has `max ≤ 15%`. Useful step forward but not done; consider widening the range or moving to LZ + eps_in joint fit.
- **Disagreement**: best eps_in has `max > 15%`. Either eps_in alone cannot close the gap, or the optimum is outside [2.0, 3.5] — recommend widening or moving to a 2D-Lindhard model.

If `onsite` and `nn` prefer different eps_in (more than ~0.3 apart in optimum), flag that the slab model is structurally limited and a single scalar `eps_in` cannot perfectly capture both channels.

## Risks / loose ends

1. **GMRES convergence at large eps_in.** As eps_in grows, the interface response is stronger and the BEM operator becomes more ill-conditioned (η = (eps_in − eps_out)/(eps_in + eps_out) → 1). If GMRES residual at any eps_in is > 1e-3, note in the report; the trend across rows tells us whether the conditioning degrades or stays under control.
2. **Source identity across rows.** Because the source is built once outside the loop, every row uses the exact same `screened_vs` input — but `BI.screened_volume_source(...)` is called inside `solve_screened_mode` and depends on `eps_in`, so the screened source is rebuilt per row. The orbital-norm `Nphi1`, `Nphi2` (computed from the *bare* `vs1`, `vs2` outside the loop) are constants across the CSV.
3. **Baseline cross-check at eps_in=2.4.** If the new sweep's `eps_in=2.4` row doesn't reproduce the 2026-05-13 k_323201 baseline to 4 decimals, the pipeline has regressed. Stop and investigate before trusting the other rows.

## Self-review

- **Placeholders.** None in the spec.
- **Internal consistency.** 6 eps_in × 2 channels = 12 CSV rows; matches schema and report table.
- **Scope.** Single sub-project at fixed LZ. No decomposition needed.
- **Ambiguity.** "Find the eps_in" means within the verdict thresholds (5% Strong / 15% Acceptable); explicit.

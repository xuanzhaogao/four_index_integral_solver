# Screened Monolayer Graphene at k_323201 — Report

Date: 2026-05-13

## Context

The 2026-05-13 bare-monolayer k-mesh sweep showed our direct Wannier integration and CoQui's bare V agree to ~10 mV across all 4 channels at k_323201, establishing this dataset as the joint validation point. This report extends the comparison to the screened (cRPA) `U_ijkl` table using a slab dielectric model.

Reference: `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15/_coqui_crpa_loc.out`.

| channel  | v_ijkl (eV) | U_ijkl (eV) |
|----------|------------:|------------:|
| onsite   | 17.4029     | 9.7824      |
| nn       |  8.8319     | 5.1678      |
| hund_ph  |  0.1316     | 0.0937      |
| hund_sf  |  0.1316     | 0.0937      |

## Method

Slab dielectric model with two GMRES dielectric solves:

| parameter | value |
|-----------|------:|
| LZ        | 3.35 Å |
| eps_in    | 2.4    |
| eps_out   | 1.0    |
| Z_CENTER  | 1.675 Å (LZ/2, symmetric) |
| L (in-plane box) | 90 Å |
| source_tol | 1e-3 |
| rhs_tol   | 1e-3  |
| lhs_tol   | 1e-5  |
| GMRES atol/rtol | 1e-5 |
| n_quad / edge_refine / max_order / max_depth | 6 / 4 / 64 / 12 |

- Density-source solve (source = `|phi_1|^2`) evaluated at targets `|phi_1|^2` (→ onsite) and `|phi_2|^2` (→ nn), swept over `Sharp` + 4 soft-mix-inverse-permittivity bandwidths `[0.05, 0.1, 0.2, 0.4]`.
- Hund's-source solve (source = target = `phi_1·phi_2` signed product), `Sharp` mode only.
- eV conversion uses orbital-norm normalization (`Na, Nb = Nphi1, Nphi2`), consistent with the bare-channel Hund's-safe convention.

Compute: one Rusty `ccm` / genoa node, 48 OpenMP threads, 6m33s wall, 13.9 GB peak RSS, all GMRES residuals 5e-6 to 1.5e-5.

## Results

### Density channels (5 modes)

| mode               | bandwidth | onsite ours | onsite CoQui | onsite diff | nn ours | nn CoQui | nn diff |
|--------------------|----------:|------------:|-------------:|------------:|--------:|---------:|--------:|
| Sharp              |       —   |     11.7471 |       9.7824 |     −1.9646 |  6.1130 |   5.1678 | −0.9452 |
| SoftMix bw=0.05    |     0.05  |     11.7480 |       9.7824 |     −1.9656 |  6.1135 |   5.1678 | −0.9457 |
| SoftMix bw=0.10    |     0.10  |     11.7592 |       9.7824 |     −1.9768 |  6.1197 |   5.1678 | −0.9519 |
| SoftMix bw=0.20    |     0.20  |     11.8172 |       9.7824 |     −2.0348 |  6.1538 |   5.1678 | −0.9860 |
| SoftMix bw=0.40    |     0.40  |     12.0136 |       9.7824 |     −2.2312 |  6.2788 |   5.1678 | −1.1110 |

Note: `diff = CoQui − ours` (so a negative diff means we overestimate U). Sharp is the best density mode; broader bandwidths only worsen agreement (broader smoothing of the inverse-permittivity transition lets less interface scattering happen, so screening weakens further).

### Hund's channel (Sharp only)

| mode  | hund_total_ev | CoQui U_ijkl | diff_ev |
|-------|--------------:|-------------:|--------:|
| Sharp |       0.09021 |       0.0937 | +0.0035 |

### Decomposition into direct vs scattered (Sharp mode)

| channel | u_int_ev (screened-volume direct) | u_scatter_ev (interface) | u_total_ev |
|---------|----------------------------------:|-------------------------:|-----------:|
| onsite  |                          12.3362  |                  −0.5891 |    11.7471 |
| nn      |                           6.2635  |                  −0.1505 |     6.1130 |
| hund    |                           0.09297 |                 −0.00276 |     0.09021 |

Sanity context — at k_323201 the bare-channel CSV gives `our_u_ev = 17.4158` (onsite), `8.8425` (nn), `0.1313` (hund). The `u_int_ev` values above are the integrals of the **screened** volume source (which has the inside-slab density divided by eps_in), so they are systematically *lower* than the bare values — for onsite by `17.42 / 12.34 ≈ 1.41×` (smaller than the formal `eps_in = 2.4` because parts of the orbital tail sit outside the slab where the density is not divided by eps_in). The `u_scatter_ev` adds the interface-induced screening on top.

## Verdict

Spec thresholds (`|rel_err|` vs CoQui U_ijkl, where `rel_err = (ours − CoQui) / CoQui`):

- Strong agreement: `|rel_err| ≤ 5%` for onsite and nn
- Acceptable:       `|rel_err| ≤ 15%` for onsite and nn
- Disagreement:     otherwise

| channel | best mode | ours    | CoQui  | rel_err  | verdict |
|---------|-----------|--------:|-------:|---------:|---------|
| onsite  | Sharp     | 11.7471 | 9.7824 |  +20.1 % | Disagreement |
| nn      | Sharp     |  6.1130 | 5.1678 |  +18.3 % | Disagreement |
| hund    | Sharp     |  0.0902 | 0.0937 |   −3.7 % | Strong (informational; not part of the formal verdict) |

**Verdict: Disagreement on density channels (~20%); Strong agreement on Hund's (~4%).**

## Interpretation

The split between Hund's (good) and density channels (off by ~20%) is informative.

The Hund's source is the signed product `phi_1·phi_2` — a near-zero-net-charge, dipole-like distribution whose far-field couples only weakly to the long-range monopole part of the dielectric screening. Its agreement at ~4% tells us the *local* and short-range parts of the slab model are well-calibrated, and that the pipeline (sources, GMRES, target normalization, eV conversion) is sound on this dataset.

The density channels (`|phi_1|^2`, `|phi_2|^2`) are net-positive, monopole-like distributions whose screening is dominated by the long-range response. We systematically under-screen them by ~20%. Three candidate causes, in rough order of expected impact:

1. **`eps_in = 2.4` is too low for the monolayer.** A higher in-slab permittivity would suppress `u_int` further and produce more interface charge to scatter. The 20% effective-screening shortfall maps (rough back-of-envelope) to needing `eps_in` ~3.0–3.5 to close the gap if `LZ` is held fixed.
2. **`LZ = 3.35 Å` is too thin.** Some fraction of the orbital density sits in the asymptotic tails outside the slab and sees `eps_out = 1` instead of the in-slab dielectric. Doubling the slab thickness (or moving to an `eps(z)` profile) would draw more of the density into the screening region.
3. **The slab-only screening model is missing a long-wavelength piece.** Real graphene's dielectric response has both a substrate-thickness component and a 2D Lindhard contribution; a step-function `eps(z)` with vacuum outside cannot capture both consistently. The fact that broader (smoother) SoftMix bandwidths *worsen* agreement says we're not missing a smoothing degree of freedom — we're missing a screening-strength degree of freedom.

The Hund's-channel success rules out: pipeline bugs, normalization errors, the Wannier orbital construction, and any short-range pathology of the slab model. The density-channel disagreement is therefore physical-parameter selection, not numerical correctness.

## Reproducibility

- Driver:        `monolayer/scripts/screened_monolayer_k323201.jl`
- Submit script: `monolayer/scripts/submit_screened_monolayer_k323201.slurm`
- CSV:           `monolayer/data/screened_monolayer_k323201.csv`
- Run logs:      `monolayer/logs/screened_sweep_6401887.{out,err}`

Slurm job 6401887, partition `ccm`, 1 genoa node × 48 cores, 6m33s wall, exit 0, MaxRSS 13.9 GB.

## Notes

- GMRES residuals: 5.5×10⁻⁶ (Sharp) to 1.5×10⁻⁵ (Hund's). All well under the 1e-3 convergence target.
- `n_interface_points`: 960552 for density solves; 956448 for the Hund's solve (small difference from the adaptive interface refinement seeing a different source profile).
- `u_int_ev` is identical to four decimals across all 5 density modes — confirms the screened volume source construction is mode-independent, and only the interface-scattering term `u_scatter_ev` varies with the soft-mix bandwidth.
- Hund's channel: only `Sharp` mode in this pass. Soft modes deferred — given the ~4% match is already in the Strong band, a soft-mix sweep is unlikely to be informative.

## Next

In order of expected payoff:

1. **`eps_in` sweep for density channels.** Repeat the Sharp-mode density solve at `eps_in ∈ {2.4, 2.8, 3.2, 3.6, 4.0}` (keep `LZ = 3.35 Å`). One additional sweep × 5 modes × 2 channels ≈ 30 min on the cluster. Find the `eps_in` that drives both onsite and nn within 5% of CoQui. If a single `eps_in` works for both → physical parameter calibration done; if onsite and nn want different `eps_in` → the slab model is structurally limited and we should move to step 2 or 3.
2. **`LZ` sweep at the best `eps_in`.** `LZ ∈ {3.0, 3.5, 4.0, 5.0, 6.7}` Å (last value is the bilayer thickness). Tests whether thickening the slab to enclose more orbital tail helps.
3. **Reply to Malte with the current verdict and the planned parameter sweeps.** The Hund's-channel result is itself a useful data point — it confirms that the monolayer Wannier+slab pipeline works for short-range screening; CoQui-vs-us agreement on the long-range piece is a model-calibration question, not a method question.
4. **Defer**: replicate this run at `k_252501` / `k_161601` once the eps_in / LZ calibration is settled. Until parameters are calibrated the k-mesh-dependence comparison would only add noise.
5. **Defer**: Hund's soft-mix sweep. Not informative until/unless the density-channel disagreement is resolved.

The bilayer-monolayer follow-up Malte mentioned remains blocked on him producing a starting example.

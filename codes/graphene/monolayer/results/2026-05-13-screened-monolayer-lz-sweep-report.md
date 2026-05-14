# Screened Monolayer Graphene — LZ Sweep at eps_in=2.4 — Report

Date: 2026-05-13

## Context

The 2026-05-13 screened-monolayer report at `k_323201` (`LZ = 3.35 Å, eps_in = 2.4`) showed +20% over-estimation on density channels (`U_onsite`, `U_nn`) while the Hund's channel agreed within 4%. The user's directive: keep `eps_in = 2.4` fixed and sweep `LZ` instead, on the hypothesis that the slab thickness controls how much of the orbital tail sits inside the screened region. This report runs that sweep.

Reference: `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15/_coqui_crpa_loc.out`.

| channel | U_ijkl (eV) |
|---------|------------:|
| onsite  | 9.7824      |
| nn      | 5.1678      |

## Method

| parameter | value |
|-----------|------:|
| eps_in    | 2.4 (fixed) |
| eps_out   | 1.0 (fixed) |
| LZ values | 2.0, 3.0, 3.35, 5.0, 8.0 Å |
| Z_CENTER  | LZ/2 (symmetric) |
| L (in-plane) | 90 Å |
| mode | Sharp only |
| solver tolerances | identical to the 2026-05-13 k_323201 run |

Per-LZ procedure: rebuild source `VolumeSource`s with `z_center = LZ/2`, run one Sharp density solve, evaluate onsite (target=vs1) and nn (target=vs2). Hund's channel skipped — already in Strong agreement at the baseline LZ=3.35.

Compute: Rusty `ccm`/genoa, 48 cores, Slurm job 6402133, 6m40s wall, 22.4 GB peak RSS.

## Results

| LZ (Å) | onsite ours | onsite rel_err% | nn ours | nn rel_err% | u_scatter onsite | u_scatter nn |
|-------:|------------:|----------------:|--------:|------------:|-----------------:|-------------:|
|  2.0   |     12.2487 |          +25.21 |  6.5618 |     +26.97  |          −0.0928 |       +0.2937 |
|  3.0   |     11.8339 |          +20.97 |  6.1942 |     +19.86  |          −0.5023 |       −0.0693 |
|  3.35  |     11.7471 |          +20.08 |  6.1130 |     +18.29  |          −0.5891 |       −0.1505 |
|  5.0   |     11.5007 |          +17.57 |  5.8771 |     +13.72  |          −0.8354 |       −0.3864 |
|  8.0   |     11.3138 |          +15.66 |  5.6937 |     +10.18  |          −1.0223 |       −0.5698 |

`u_total = u_int + u_scatter`. `rel_err = 100 * (ours − CoQui) / CoQui`. All GMRES residuals 3–10 ppm, all converged.

### Direct part vs interface part — Sharp mode

`u_int_ev` is essentially flat across LZ (12.336 ± 0.005 for onsite, 6.263 ± 0.0001 for nn for LZ ≥ 3) — the screened volume source's effective integral is independent of slab thickness once the orbital is enclosed. All the LZ-dependence sits in `u_scatter_ev` (interface contribution), which grows monotonically in magnitude as LZ increases. The (small, positive) nn `u_scatter` at LZ = 2 is a boundary artifact: the orbital extends well outside the 2 Å slab, and the interface charge distribution there does not have its standard "screening" topology.

## Verdict

Best LZ by `max(|rel_err_onsite|, |rel_err_nn|)`:

| LZ (Å) | max abs rel_err |
|-------:|----------------:|
| 2.0    |          27.0 % |
| 3.0    |          21.0 % |
| 3.35   |          20.1 % |
| 5.0    |          17.6 % |
| 8.0    |    **15.7 %**   |

**Best LZ in the swept range: 8.0 Å, with max abs rel_err = 15.7%.**

Spec thresholds:
- Strong: best max ≤ 5%
- Acceptable: best max ≤ 15%
- Disagreement: otherwise

**Verdict: Disagreement (just barely — 15.7% vs the 15.0% Acceptable threshold).** The nn channel at LZ=8 (10.2%) is well inside Acceptable; the onsite channel (15.7%) sits ~0.7 percentage points above. Without ambiguity-aware thresholds the formal verdict is Disagreement, but the data is squarely on the boundary.

## Interpretation

The four key observations:

1. **The rel_err curves are strictly monotonic in LZ across [2.0, 8.0] Å.** No interior optimum exists in this range. Going to *smaller* LZ — which I included as a sanity probe at the user's prompt — gives strictly worse agreement, confirming that the direction of improvement is toward larger LZ and ruling out any sub-graphite slab-thickness sweet spot.

2. **`u_int` saturates quickly; all LZ-dependence is in `u_scatter`.** From LZ = 3 onward `u_int_ev` is constant to 4 decimals. The screened volume source `BI.screened_volume_source` integrates the inside-slab density (divided by `eps_in`) plus the outside-slab density (unscaled); for LZ ≥ 3 essentially all of the relevant orbital density is inside the slab, so further widening the slab just adds vacuum to the box without changing the source's contribution. What actually changes is the interface-induced scattering term: the dielectric step is pushed further from the orbital core, allowing more interface area to participate in the screening response.

3. **nn converges faster than onsite.** At LZ = 8 the nn channel is already at 10.2 % (Acceptable) while onsite is at 15.7 %. Onsite is dominated by short-range density-density coupling, where the slab-model q→0 limit is most sensitive to the chosen `eps_in`; nn averages over a larger spatial range where finite-q effects partially compensate. The fact that the two channels have systematically different `eps_in`-sensitivity is itself a signal that a single scalar `eps_in` cannot perfectly capture the screening for both.

4. **Extrapolating linearly in LZ, neither channel reaches the 5% Strong band within physically reasonable slab thicknesses.** The onsite rel_err drops from 25.2 % (LZ=2) → 15.7 % (LZ=8) — a ~0.16 %/Å slope. To reach 5 % would require LZ ≈ 75 Å, at which point the "slab" is no longer a physical monolayer description. The slab model with `eps_in = 2.4` has an asymptotic floor on its agreement; closing the last few percent requires either a higher `eps_in` or a different screening model (e.g., 2D Lindhard, `eps(z)` profile).

## Reproducibility

- Driver:        `monolayer/scripts/screened_monolayer_lz_sweep.jl`
- Submit script: `monolayer/scripts/submit_screened_monolayer_lz_sweep.slurm`
- CSV:           `monolayer/data/screened_monolayer_lz_sweep.csv`
- Run logs:      `monolayer/logs/lz_sweep_6402133.{out,err}`

Slurm job 6402133, partition `ccm`, 1 genoa node × 48 cores, 6m40s wall, exit 0, MaxRSS 22.4 GB.

## Next

In order of expected payoff:

1. **Open the `eps_in` sweep at the best LZ (= 8 Å).** Fix `LZ = 8.0` and sweep `eps_in ∈ {2.4, 2.8, 3.2, 3.6, 4.0}`. The LZ sweep has done the work it can; the remaining +15.7 % / +10.2 % shortfall is too large to be closed by widening the slab further (would need LZ ~ 75 Å for onsite), so the next physical lever is the in-slab permittivity. Estimated compute is the same as this sweep (~7 min).
2. **Reply to Malte** with this LZ-sweep table and the planned eps_in follow-up. The monotonic-in-LZ pattern and the saturation of `u_int` are independently useful for him to see — they characterize what a step-function slab model can and cannot do for the monolayer.
3. **Defer**: replicate this sweep at k_252501 / k_161601. The k-mesh-dependence question is downstream of getting the model parameters right at one k-mesh.
4. **Defer**: try an asymmetric `Z_CENTER` (e.g., orbital pinned at a fixed z while LZ varies). Symmetric placement is the natural physical setup; an asymmetric variant is only interesting if a substrate model is in play.
5. **Defer**: `eps(z)` profile or 2D-Lindhard correction. These are larger code changes and only justified if the eps_in sweep also stalls.

The bilayer-monolayer follow-up Malte mentioned remains blocked on him producing a starting example.

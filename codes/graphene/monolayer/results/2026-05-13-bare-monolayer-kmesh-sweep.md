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

**Verdict: Rejected.** k_161601 confirms the hypothesis to 7 mV, but k_252501 and k_323201 exceed the 0.10 eV "Partial" threshold (|−0.128| and |+0.278| eV respectively). The Madelung correction alone does not reconcile our direct Wannier integral with CoQui's bare v₀ across all three k-meshes.

### Interpretation

The pattern of the three residuals is informative:

- **k_161601** (Madelung +1.041 eV) — residual −7 mV. The original 2026-04-23 hypothesis match was not coincidence; on the original Wannier dataset the Madelung-correction model is essentially exact.
- **k_252501** (Madelung +0.184 eV) — residual −128 mV. Our integral landed at 17.382 eV (vs CoQui 17.438). The Madelung correction overshoots the observed gap by ~70%.
- **k_323201** (Madelung −0.290 eV) — residual +278 mV. Our integral landed at 17.416 eV, slightly above CoQui 17.404. The Madelung correction would predict our value should be ~0.290 eV higher than CoQui's, but it is barely above it.

The k_161601 onsite result re-runs to 16.4004 eV — identical (to 4 decimals) to the 2026-04-23 value, confirming run reproducibility and ruling out pipeline drift as a cause.

Three possible explanations for the failure on the new k-meshes:

1. **Wannier construction differences.** Malte flagged in his email that "the Wannier constructions for the monolayer case were not as perfect as they can be," and is currently re-optimizing them. The MLWFs at k_252501 and k_323201 may differ in shape from the well-tested k_161601 set in ways that affect both v₀ and the Madelung formula's applicability.
2. **Madelung convention mismatch.** CoQui's long-wavelength correction in the bare V uses the Gygi-Baldereschi prescription (visible in the original `_coqui_thc_crpa.out`: `Treatment of long-wavelength divergence in bare V: gygi`). The numerical Madelung values Malte provided may be reported in a convention or normalization that doesn't map directly onto `CoQui_v0 = Wannier_v0 + Madelung` for these k-meshes.
3. **XSF resolution.** Cusp under-resolution on the 150×150×192 grid could bias our integral by amounts that grow with the size of v₀; this is the same mechanism noted in the 2026-04-23 report as the residual ~5% underestimate on k_161601. A 3×3×1 regeneration at the same nominal grid resolution would increase the effective resolution per unit cell by ~2.78×.

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

Given the **Rejected** verdict, the recommended follow-ups (in order):

1. **Confirm the Madelung convention with Malte.** Ask whether the three values he supplied (1.041, 0.1842, −0.290 eV) are intended as `Madelung = CoQui_v0 − Wannier_v0` directly, or whether they involve a normalization (e.g., per unit-cell volume) that we are misapplying. Forward this report as the concrete observation.
2. **Wait for his re-optimized Wannier set.** If the Wannier construction is being improved anyway, repeat the k_252501 / k_323201 sweep on the new XSFs before drawing structural conclusions.
3. **Defer the XSF-resolution decision (3×3×1 regeneration).** The largest residual (k_323201 = +0.278 eV) is comparable to the ~5% raw cusp-under-resolution effect on k_161601, but the *opposite-sign* k_252501 residual rules out a pure-resolution explanation. Resolution is therefore not the leading candidate, and a 3×3×1 regeneration should wait until the Wannier set is finalized.

The bilayer-monolayer follow-up Malte mentioned remains blocked on him producing a starting example.

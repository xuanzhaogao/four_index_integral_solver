# Bare Bilayer Graphene Interactions — Report

Date: 2026-05-15 (revised 2026-05-16 with per-pair CoQui reference)

## Reference (CoQui cRPA; Malte Rösner)

Source directory: `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/bilayer/k_161601_nb_288_c_15/`

Two reference levels were extracted:

1. **Orbital averages** from `_coqui_thc_crpa.out` (`Downfolded Coulomb Elements Summary`):
   - intra-orbital              = 17.669390 eV
   - inter-orbital              =  5.785749 eV
   - Hund's (spin-flip)         =  0.041867 eV
   - Hund's (pair-hopping)      =  0.041867 eV

2. **Full bare 4-index tensor** from `crpa.mbpt_clean.h5`, dataset
   `/scf/iter0/downfolded_model/Vloc_abcd` (Fortran shape `{4,4,4,4,2}` with last
   axis real/imag, stored in Hartree). Index convention `Vloc[a,b,c,d]` =
   `(ab|cd)` is verified by reproducing the orbital averages from the text
   output to 6 decimal places (see `bilayer/scripts/extract_coqui_bare_tensor.jl`).
   The imaginary part of `Vloc_abcd` is exactly zero, confirming the real-orbital
   convention.

   Mapping to our channels:

   | channel    | Vloc index pattern |
   |------------|--------------------|
   | onsite (i)            | `Vloc[i,i,i,i]` |
   | density–density (i,j) | `Vloc[i,i,j,j]` |
   | Hund spin-flip (i,j)  | `Vloc[i,j,j,i]` |
   | Hund pair-hop  (i,j)  | `Vloc[i,j,i,j]` |

## This Run

Orbitals: `graphene_00001.xsf` … `graphene_00004.xsf` (150×150×225 grid over a
12.24 × 10.60 × 18.35 Å Wannier90 supercell).

Wannier centres (Å, from `graphene_centres.xyz`):

| orbital | x | y | z | site |
|--------:|---:|---:|---:|---|
| 1 | 0.0000 | 0.0000 |  7.5016 | layer 1, sublattice A |
| 2 | 0.0000 | 1.4232 |  7.4968 | layer 1, sublattice B |
| 3 | 0.0000 | 1.4232 | 10.8532 | layer 2, dimer site (above B1) |
| 4 | 1.2325 | 0.7116 | 10.8484 | layer 2, non-dimer (hollow over layer 1) |

Method: identical to the monolayer pipeline — direct real-space `1/|r−r'|` via
`BI.TKM3D.ltkm3dc` on `BI.VolumeSource` objects built from squared (density) and
signed-product (Hund's) orbital grids. The monolayer driver
`compute_all_bare_channels(; orbital_1, orbital_2, …)` is invoked on every
distinct orbital pair `(i, j), i ≤ j`; centering uses orbital `i`'s density
centroid for both members of the pair, so the in-plane / inter-layer offsets
between the two centres are preserved. No dielectric solve.

Parameters: `source_tol = 1e-3`, `volume_tol = 1e-3`, `mirror_pad_level = 0`.

Pipeline validation: the underlying TKM3D path is the same one independently
verified against the analytic two-Gaussian Coulomb integral to 4.8 × 10⁻⁴
relative error in the monolayer testset.

## 1:1 Per-Pair Comparison Against `Vloc_abcd`

All 22 channels we compute, each lined up with the corresponding CoQui entry
from `Vloc_abcd`:

| (i, j) | class       | channel | ours (eV) | CoQui (eV) | Δ (eV) | rel.err |
|--------|-------------|---------|----------:|-----------:|-------:|--------:|
| (1, 1) | onsite      | onsite  | 17.6867 | 17.6526 | +0.0341 | +0.193 % |
| (2, 2) | onsite      | onsite  | 17.7189 | 17.6859 | +0.0330 | +0.186 % |
| (3, 3) | onsite      | onsite  | 17.7187 | 17.6860 | +0.0327 | +0.185 % |
| (4, 4) | onsite      | onsite  | 17.6845 | 17.6531 | +0.0315 | +0.178 % |
| (1, 2) | in-plane    | nn      |  8.8710 |  8.8586 | +0.0124 | +0.141 % |
| (3, 4) | in-plane    | nn      |  8.8705 |  8.8586 | +0.0119 | +0.134 % |
| (1, 3) | inter-layer | nn      |  4.0754 |  4.1393 | −0.0640 | **−1.545 %** |
| (1, 4) | inter-layer | nn      |  4.0805 |  4.1441 | −0.0636 | **−1.534 %** |
| (2, 3) | inter,dimer | nn      |  4.5082 |  4.5746 | −0.0664 | **−1.451 %** |
| (2, 4) | inter-layer | nn      |  4.0753 |  4.1393 | −0.0640 | **−1.547 %** |
| (1, 2) | in-plane    | hund_sf | 0.122452 | 0.122142 | +3.10e−4 | +0.254 % |
| (3, 4) | in-plane    | hund_sf | 0.122445 | 0.122140 | +3.04e−4 | +0.249 % |
| (1, 3) | inter-layer | hund_sf | 0.000566 | 0.000565 | +1.05e−6 | +0.186 % |
| (1, 4) | inter-layer | hund_sf | 0.000935 | 0.000935 | −4.63e−7 | −0.049 % |
| (2, 3) | inter,dimer | hund_sf | 0.004862 | 0.004854 | +8.04e−6 | +0.166 % |
| (2, 4) | inter-layer | hund_sf | 0.000566 | 0.000565 | +9.66e−7 | +0.171 % |

`hund_ph = hund_sf` to machine precision for every pair (real-orbital identity
`V_abba = V_abab`); the six `hund_ph` rows in `bare_bilayer_vs_coqui.csv` are
numerically identical to the `hund_sf` rows and are omitted from the table
above for brevity.

### Orbital-average summary (derived from the same rows)

| channel               | this work (eV) | CoQui (eV) | rel. err |
|-----------------------|---------------:|-----------:|---------:|
| intra-orbital         | 17.7022 | 17.6694 | +0.186 % |
| inter-orbital         |  5.7468 |  5.7857 | −0.673 % |
| Hund's (spin-flip)    |  0.04199 |  0.04187 | +0.248 % |
| Hund's (pair-hopping) |  0.04199 |  0.04187 | +0.248 % |

## Observations

1. **All 22 channels agree with CoQui to within ±1.55 %**, and 16/22 within
   ±0.3 %. The only systematic deviation is the inter-layer density–density
   channel.

2. **Uniform inter-layer NN deficit of −1.45 % to −1.55 %.** The four
   inter-layer density–density values — including the dimer pair (2,3) at
   3.35 Å and the three non-dimer pairs at ≈ 3.64 Å — all sit below CoQui by a
   nearly identical fraction. The flatness across geometrically different pairs
   suggests a small systematic in how the orbital z-extent is captured by our
   finite-grid `|φ|²` density rather than a per-distance kernel issue. A
   z-direction convergence study (`mirror_pad_z` or denser z grid) is the
   natural follow-up; a 1.5 % effect at 4 eV is ≈ 60 meV, which is within the
   bandwidth resolution we use for the screened-channel comparison and so does
   not block downstream work.

3. **Onsite and in-plane NN are +0.13 to +0.19 % *above* CoQui** — the opposite
   sign of the inter-layer deficit. The orbital averages partially cancel these
   contributions (intra +0.19 %, inter −0.67 %), making the orbital-average
   numbers look tighter than the per-pair structure actually is.

4. **Hund's channels per-pair match within ±0.3 %**, including the dimer pair
   to 8 µeV out of 4.86 meV (one part in 600). The exchange-like Hund's integral
   is dominated by the orbital cores where our grid is well-resolved, and there
   is no inter-layer-specific deficit; this is consistent with the density-only
   nature of the −1.5 % inter-layer artifact.

5. **Per-pair physics is clean.** Dimer pair (2,3), |Δr| = 3.35 Å, has the
   largest inter-layer density-density V (4.51 eV) and the largest inter-layer
   Hund's J (4.86 meV); the three non-dimer pairs (|Δr| ≈ 3.64 Å) cluster at
   V ≈ 4.08 eV and J ≈ 0.69 meV. Inter-layer J is suppressed by ≈ 70× relative
   to in-plane J (122 meV), as expected from the exponential decay of the `pz`
   tail at 3.35–3.64 Å with `ξ ~ 0.7 Å`.

## Verdict

- 1:1 per-pair agreement with CoQui's `Vloc_abcd`: ±1.55 % worst case, ±0.3 %
  on 16/22 channels.
- A single small systematic remains: a uniform ≈ −1.5 % deficit on the four
  inter-layer density–density channels, of order 60 meV. Hund's, onsite, and
  in-plane channels are all sub-0.3 %.
- The orbital-average numbers from the text output (intra 17.67, inter 5.79,
  J 0.0419) are reproduced to ±0.7 % and to ±0.3 % respectively, with a
  cancellation between the +0.19 % intra and the −1.5 % inter-layer entering
  the inter-orbital average.

## Files

- Our per-pair output: `bilayer/data/bare_bilayer_graphene.csv` (22 rows).
- CoQui reference per-pair (extracted from h5): `bilayer/data/coqui_bare_vloc_abcd.csv` (22 rows).
- Joined 1:1 comparison: `bilayer/data/bare_bilayer_vs_coqui.csv` (22 rows with Δ and rel.err).
- Bilayer driver: `bilayer/src/BilayerBareIntegrals.jl`.
- Run script: `bilayer/scripts/bare_bilayer_graphene.jl`.
- Reference extractor: `bilayer/scripts/extract_coqui_bare_tensor.jl` (also doubles as the index-convention verifier).
- Slurm script: `bilayer/scripts/submit_bare_bilayer_graphene.slurm` (this run completed locally in ~5 min on `ccmlin078` with `JULIA_NUM_THREADS=16`).
- Run log: `bilayer/logs/bare_bilayer_20260515-152456.log`.

## Next

- (Optional) z-convergence check on the inter-layer NN channel: re-run a single
  inter-layer pair at a tighter `source_tol`, and/or with a `mirror_pad_z`
  diagnostic analogous to the monolayer `mirror_pad_xy`, to confirm the
  −1.5 % deficit is a finite-grid z-resolution effect.
- Bilayer screened-channel comparison can now hook into this per-pair table
  once the dielectric model is wired up — `Vloc_abcd` has its screened sibling
  `Uloc_wabcd(ω)` in the same h5, so the same per-pair join works for `W`.

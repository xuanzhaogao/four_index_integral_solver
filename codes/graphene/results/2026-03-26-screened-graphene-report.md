# Screened Graphene Results Report

Date: 2026-03-26

## Scope

This report summarizes the verified results currently available in this repository for the screened graphene density study and the new full screened-orbital solve workflow.

The current state splits into two parts:

1. The density-analysis workflow is complete and has generated figures.
2. The full screened Hubbard/orbital workflow is implemented, tested, and now has a validated sharp-screening production run for the bilayer-as-slab geometry, yielding `U_00` through `U_05`.

## Problem Setup

The current bilayer-slab production workflow is configured with:

- `Lx = Ly = 90.0`
- `Lz = 6.7`
- `eps_in = 2.4`
- `eps_out = 1.0`
- orbital centers shifted to `z = d/4 = 1.675`

The shared centering shift applied to the orbital densities is:

- `(-4.9383594533135144e-5, -1.1431743832672156e-5, -6.325157062349036)`

At `SOURCE_TOL = 1e-3`, the centered orbital sources resolve to:

- orbital 1 source points: `314860`
- orbital 2 source points: `315606`
- orbital 1 normalization: `84.17965021252988`
- orbital 2 normalization: `84.18004419096994`

## Density Analysis

### Generated Figures

- `figs/screened_density_xslice.png`
- `figs/screened_density_kz_decay.png`
- `figs/screened_density_kz_decay_refined_qz.png`
- `figs/screened_density_kz_decay_global_nufft_qz.png`

### Centering Check

After the centering fix, the screened-density centroid is numerically at the box center:

- centroid = `(1.0542272718855331e-16, 1.9503204529882363e-16, 1.0529460035413301e-15)`

This confirms that the `x = 0` slice plots and the `\hat{\rho}(0,0,k_z)` curves are using the same centered density.

### Slice Results

The `x = 0` `yz` slice is plotted using `log10(value)` rather than the raw density. That change was needed because the original linear heatmap hid the small-density structure away from the orbital core. The log-scale slice makes the near-interface tails visible while preserving the centered placement of the density in the box.

### Fourier-Decay Results

The `k_z` study fixes `k_x = k_y = 0` and compares `|\hat{\rho}(0,0,k_z)|` for:

- `SharpScreening()`
- `SoftMixInversePermittivity(b)` with `b = 0.05, 0.1, 0.2, 0.4`

The qualitative trend is consistent across the generated figures:

- smaller `b` gives a slower-decaying spectral tail
- `b = 0.05` is the most difficult case
- the coarse-grid behavior indicates that the narrow surface layer is underresolved when `b` becomes comparable to or smaller than the `z` mesh spacing

### Refined `Q(z)` Experiments

Two separate refined `Q(z)` experiments were carried out for the `k_z` tail:

1. periodic type-1/type-2 NUFFT upsampling
2. padded/global type-1/type-2 NUFFT upsampling

The padded/global construction reduces wrap-around ringing relative to the purely periodic reconstruction for `b = 0.1, 0.2, 0.4`, but it does not materially improve the smallest-bandwidth case.

Normalized tail values at the largest sampled `k_z`:

| `b` | periodic refined `Q(z)` | padded/global NUFFT `Q(z)` |
| --- | ---: | ---: |
| `0.05` | `2.116851517633355e-3` | `2.139159318511159e-3` |
| `0.1` | `3.335779482031743e-4` | `3.156327295004313e-4` |
| `0.2` | `1.7911491158156037e-4` | `1.2563963029966617e-4` |
| `0.4` | `1.74063737938745e-4` | `1.2497173321172067e-4` |

Interpretation:

- the periodic Fourier interpolant is not a good reconstruction for `Q(z)` because `Q(z)` is not smoothly periodic across the box boundary
- padding reduces the boundary-continuation artifact, but it cannot recover surface structure that the original mesh never resolved
- the slow tail for very small `b` is therefore partly physical and partly mesh-limited

## Full Screened Orbital Solve Workflow

### Implemented Pipeline

The full workflow has been implemented in:

- `src/ScreenedOrbitalSolve.jl`
- `scripts/screened_hubbard_graphene.jl`

The solve path now does the following:

1. load and square the orbital densities from XSF
2. apply one shared centering shift from orbital 1 to both orbitals
3. construct the dielectric box interface
4. screen the source density with either sharp screening or softmix
5. assemble the RHS and corrected LHS
6. solve for the surface density `sigma` with GMRES
7. evaluate the direct volume potential `u_int`
8. evaluate the scattered potential `u_scatter`
9. integrate `u_int + u_scatter` against the six target orbitals `U_00`, `U_01`, `U_02`, `U_03`, `U_04`, and `U_05`

The script is configured to write:

- `data/screened_hubbard_graphene.csv`

and supports:

- `julia --project=. scripts/screened_hubbard_graphene.jl sharp`
- `julia --project=. scripts/screened_hubbard_graphene.jl 0.05 0.1 0.2 0.4`

### Current Status

The workflow is implemented, covered by tests, and now has a completed local `sharp` production run for the bilayer-slab case. The output file is present:

- `data/screened_hubbard_graphene.csv`

The sharp-screening interaction values are:

| Pair | `u_int_ev` | `u_scatter_ev` | `u_total_ev` |
| --- | ---: | ---: | ---: |
| `U_00` | `7.18640` | `1.23061` | `8.41701` |
| `U_01` | `3.69732` | `1.15564` | `4.85296` |
| `U_02` | `2.33554` | `1.05490` | `3.39044` |
| `U_03` | `2.03792` | `1.01758` | `3.05550` |
| `U_04` | `1.57792` | `0.93555` | `2.51347` |
| `U_05` | `1.39494` | `0.89469` | `2.28963` |

This means the repo now contains validated analysis figures, validated solve infrastructure, and a completed sharp-screening interaction table for the bilayer-slab geometry.

### Labeling Note

For the current slab workflow, `U_00` through `U_05` are labeled to follow the shell ordering used in Rosner et al. for bilayer graphene:

- `U_00`: on-site
- `U_01`: nearest neighbor
- `U_02`: next-nearest neighbor
- `U_03`: third shell
- `U_04`: fourth shell representative
- `U_05`: fifth shell representative

This is a labeling alignment only. The present workflow is still a dielectric-slab approximation with scalar box parameters `eps_in = 2.4` and `eps_out = 1.0`. The paper instead uses AB-stacked bilayer graphene together with a momentum-dependent effective dielectric function `ε_eff^2D(q)`, so the current values should not be interpreted as a full WFCE/cRPA reproduction of the paper’s bilayer interaction channels.

## Verification

The current test suite passes with:

```bash
julia --project=. test/runtests.jl
```

Relevant tested items include:

- centering and shifted slicing
- log10 slice transform behavior
- refined `Q(z)` spectral paths
- global NUFFT boundary-ringing reduction
- pair-target construction for `U_00/U_01/U_02/U_03/U_04/U_05`
- target-potential integration
- synthetic validation of the volume potential evaluation

## Bottom Line

The density-analysis results are complete enough to support two conclusions:

1. The centered density and `x = 0` log-scale slices behave as expected in the dielectric box.
2. Smaller softmix bandwidths produce slower `k_z` decay, and simple NUFFT upsampling of the reduced `Q(z)` only partially mitigates numerical artifacts because the main limitation is the original mesh resolution near the interface.

The full screened-orbital solver is now in place and has produced a verified sharp-screening bilayer-slab interaction table through `U_05`. Additional softmix production sweeps remain optional follow-on runs rather than a missing core result.

# Screened Graphene Results Report

Date: 2026-03-26

## Scope

This report summarizes the verified results currently available in this repository for the screened graphene density study and the new full screened-orbital solve workflow.

The current state splits into two parts:

1. The density-analysis workflow is complete and has generated figures.
2. The full screened Hubbard/orbital workflow is implemented and tested, but a real-data production run has not completed locally yet, so there is no validated CSV of `U_00`, `U_01`, `U_02`, and `U_03` values in this repo at the time of this report.

## Problem Setup

The full-solve workflow is configured with:

- `Lx = Ly = 90.0`
- `Lz = 2.8`
- `eps_in = 2.4`
- `eps_out = 1.0`

The shared centering shift applied to the orbital densities is:

- `(-4.9383594533135144e-5, -1.1431743832672156e-5, -8.000157062349036)`

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
9. integrate `u_int + u_scatter` against the four target orbitals `U_00`, `U_01`, `U_02`, and `U_03`

The script is configured to write:

- `data/screened_hubbard_graphene.csv`

and supports:

- `julia --project=. scripts/screened_hubbard_graphene.jl sharp`
- `julia --project=. scripts/screened_hubbard_graphene.jl 0.05 0.1 0.2 0.4`

### Current Status

The workflow is implemented and covered by tests, but there is no completed real-data output file in this repository yet:

- `data/screened_hubbard_graphene.csv` is absent
- a local `sharp` production run did not finish within the interactive session, so no validated orbital interaction table is available yet

This means that the repo currently contains validated analysis figures and validated solve infrastructure, but not final Hubbard numbers from a completed production solve.

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
- pair-target construction for `U_00/U_01/U_02/U_03`
- target-potential integration
- synthetic validation of the volume potential evaluation

## Bottom Line

The density-analysis results are complete enough to support two conclusions:

1. The centered density and `x = 0` log-scale slices behave as expected in the dielectric box.
2. Smaller softmix bandwidths produce slower `k_z` decay, and simple NUFFT upsampling of the reduced `Q(z)` only partially mitigates numerical artifacts because the main limitation is the original mesh resolution near the interface.

The full screened-orbital solver is now in place, but the actual production interaction values still need a completed real-data run.

# Screened Density kz-Slab Spectrum Design

Date: 2026-03-26

## Goal

Replace the old `screened_density_kz_decay.png` observable based on `\hat{\rho}(0,0,k_z)` with a 3D type-1 NUFFT spectrum reduced by

`max_{k_x,k_y} |\hat{\rho}(k_x,k_y,k_z)|`.

## Requirements

- Use the same TKM3D-style mode construction as the solver path:
  - `kmax = π / h` from the estimated source spacing
  - centered mode axes built from TKM3D-style `Δk_i`
- Compare the same screening cases:
  - `Sharp`
  - `SoftMixInversePermittivity(b)` for `b = 0.05, 0.1, 0.2, 0.4`

## Design

1. Build the screened density vector on the centered `VolumeSource`.
2. Form source charges as `weights .* rho`.
3. Use a 3D type-1 NUFFT on a TKM3D-style mode box.
4. For each nonnegative `k_z` slab inside the spherical `kmax` cutoff, take the maximum coefficient magnitude over `(k_x, k_y)`.
5. Normalize by the `k_z = 0` slab maximum and overwrite `figs/screened_density_kz_decay.png`.

# Screened Density z-Upsampling Design

Date: 2026-03-26

## Goal

Add a separate decay figure that tests whether globally upsampling the original density only in `z` changes the TKM3D-style screened-density spectral decay.

## Requirements

- Keep the current `figs/screened_density_kz_decay.png` unchanged.
- Start from the original centered density, before screening.
- For each fixed `(x, y)` column:
  - use a type-1 NUFFT in `z` to obtain Fourier modes
  - use a type-2 NUFFT to evaluate on a finer uniform `z` grid
- Try upsampling factors `2` and `3`.
- After reconstruction, apply the same screening cases as before:
  - `Sharp`
  - `SoftMixInversePermittivity(b)` for `b = 0.05, 0.1, 0.2, 0.4`
- Evaluate the same observable as `screened_density_kz_decay.png`:
  - `max_{k_x,k_y} |\hat{\rho}(k_x,k_y,k_z)|`
  - with the same TKM3D-style `k_i` construction and `kmax`

## Design

1. Recover the tensor-product `x`, `y`, `z` axes and `weights`/`density` arrays from the centered `VolumeSource`.
2. Reconstruct each `(x, y)` density column globally in `z` using a type-1/type-2 NUFFT round-trip at a finer uniform `z` spacing.
3. Scale the `z` weights by the new `Δz` so the refined source preserves volume measure.
4. Rebuild a refined `VolumeSource` with unchanged `x/y` axes and refined `z`.
5. For each screening case and each upsampling factor in `{1, 2, 3}`, compute the screened density and the same TKM3D slab-max decay used by the current figure.
6. Save a new multi-panel figure so each screening case shows three curves: original, `z x 2`, and `z x 3`.

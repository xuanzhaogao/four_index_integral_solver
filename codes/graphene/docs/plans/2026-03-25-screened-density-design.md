# Screened Density Slice And kz Decay Design

**Goal:** Compare sharp screening against soft inverse-permittivity mixing on the graphene density by plotting the `x = 0` `yz` slice and the decay of `\hat{\rho}(0, 0, k_z)` for several softmix bandwidths.

## Scope

- Reuse the existing XSF input and dielectric-box parameters from `scripts/screened_density_tkm.jl`.
- Compare `BoundaryIntegral.SharpScreening()` with `BoundaryIntegral.SoftMixInversePermittivity(b)` for a small bandwidth sweep.
- Save one figure for real-space slices and one figure for spectral decay.

## Real-Space Slice

The underlying orbital density lives on an oblique XSF grid, so a Cartesian `x = 0` plane is not aligned with a single stored grid index. The slice will therefore be evaluated by:

1. Building the affine map of the XSF datagrid.
2. Sampling a Cartesian `yz` grid on the plane `x = 0`.
3. Interpolating the unscreened density on that plane with BoundaryIntegral's trilinear datagrid evaluator.
4. Applying the screening factor pointwise on the sampled plane for each mode.

This keeps the slice physically aligned with Cartesian `x`.

## Spectral Decay

For `\hat{\rho}(0, 0, k_z)`, only the `z` dependence matters. The full 3D screened density will be computed on the existing `VolumeSource`, then reduced to plane charges along `z`:

`\hat{\rho}(0, 0, k_z) = \sum_{z_j} Q(z_j) e^{-i k_z z_j}`

where `Q(z_j)` is the quadrature-weighted charge on the `z_j` plane. This avoids an unnecessary 3D FFT and works directly with the stored non-orthogonal volume grid.

## Outputs

- `figs/screened_density_xslice.png`
- `figs/screened_density_kz_decay.png`

## Testing

- Unit test the Cartesian slice evaluation on a synthetic linear datagrid, where trilinear interpolation is exact.
- Unit test the `\hat{\rho}(0, 0, k_z)` reduction against direct quadrature on a small synthetic volume source.

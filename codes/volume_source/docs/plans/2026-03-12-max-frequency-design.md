# Max Frequency Design

## Goal

Add a standalone analysis script that estimates the smallest radial cutoff `k_cut`
for the squared graphene orbital density such that

`max_{|k| > k_cut} |F(k)| / max_k |F(k)| < tol`.

The orbital to analyze is
`density_data/graphene_00001_5x5x1_shifted.xsf`, loaded the same way as the
existing orbital RHS scripts.

## Chosen Approach

Use the raw XSF grid directly instead of a pruned `VolumeSource`.

1. Load the XSF datagrid with `BoundaryIntegral.read_xsf`.
2. Square `datagrid.values` in place, matching the existing orbital RHS script.
3. Convert the full regular grid into physical-space source coordinates and
   preweighted cell masses.
4. Use a single 3D type-1 NUFFT to evaluate the centered Fourier coefficients on
   the full grid Nyquist box implied by the original grid spacings.
5. Compute `|F(k)|` and scan the radial shells to find the smallest `k_cut` whose
   remaining tail stays below the requested relative pointwise tolerance.

This is the most faithful answer to the user’s question because it measures the
actual orbital file’s spectral tail instead of the tail of a sparsified source
representation.

## Scope

The script will:

- define a user-editable `tol`
- print grid shape, box vectors, spacings, and axis-wise Nyquist limits
- print the global spectral maximum
- print the estimated `k_cut`
- print a few tail checkpoints so the spectral decay can be inspected manually

The script will not add plotting or CLI parsing in this change.

## Implementation Notes

- The script should be structured as reusable helper functions plus a
  `main()` entrypoint so the cutoff logic can be tested directly.
- The relative tail metric is
  `tail_ratio(k_cut) = max_{|k| > k_cut} |F(k)| / max_k |F(k)|`.
- The key helper will operate only on arrays of radii and magnitudes, making it
  easy to verify with a small synthetic test before wiring in the orbital data.

## Verification

- Add a focused test for the cutoff helper on a synthetic spectrum where the
  correct cutoff is obvious.
- Run that test in red/green order.
- Run the full script on the graphene orbital and report the measured `k_cut`.

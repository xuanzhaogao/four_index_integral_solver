# Screened Gaussian Refined-NUFFT Design

## Goal

Resolve the Fourier-space decay of the `erf`-screened Gaussian for narrow
bandwidths by locally refining the source sampling near `x = 0` and evaluating
the spectrum with a weighted type-1 NUFFT.

## Chosen Approach

Keep the work in `scripts/screened_gauss.jl`, but replace the current
uniform-grid decay study with a nonuniform quadrature rule. The source interval
will remain finite and centered, but the cell widths will be made much smaller
near `x = 0`, where the `erf` screen differs from the step function.

The sampled transform will be approximated by

`fhat(k_m) ≈ Σ_j w_j f(x_j) exp(-i k_m x_j)`

where `x_j` are nonuniform midpoint nodes and `w_j` are the corresponding cell
widths. A type-1 NUFFT will then evaluate these modal sums efficiently on a
uniform `k` grid.

This is the right fix for the current issue because the screening-layer
information is localized near `x = 0`, while the Gaussian tail away from `0`
does not need the same resolution.

## Scope

The update will:

- add a locally refined midpoint quadrature helper on a symmetric finite interval
- add a weighted type-1 NUFFT helper for nonuniform nodes
- update the decay study to use the refined quadrature and weighted NUFFT
- keep the unscreened Gaussian as an exact reference on the same `k` grid
- keep the step-screened and `erf`-screened curves in the same figure
- preserve the separate convolution study already added

The update will not remove the finite-domain truncation; it only changes how the
source is sampled inside that interval.

## Numerical Design

- Outer interval: the same truncated Gaussian domain already used for the decay
  and convolution studies
- Quadrature nodes: midpoint samples on piecewise-uniform cells
- Inner refined zone: centered at `x = 0` with a user-editable half-width
- Inner spacing: coarse spacing divided by a user-editable refinement factor
- Modal grid: choose the NUFFT mode extent from the smallest cell width unless
  the user overrides it

## Plot Design

- Horizontal axis: `k * sigma` instead of `k / (pi / dx)`, because the refined
  nonuniform sampling has no single global `dx`
- Vertical axis: `|fhat(k)|` on a log scale
- Reference curves: `unscreened` and `step`
- Smooth-screen curves: one line per `erf` bandwidth

## Testing

Add a focused test for the refined quadrature and weighted NUFFT on the plain
Gaussian, whose exact Fourier transform is already available. Also keep the
study-shape test so the step and `erf` curves continue to share the same
positive-frequency axis and contain nonnegative magnitudes.

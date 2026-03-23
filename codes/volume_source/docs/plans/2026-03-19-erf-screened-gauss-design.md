# Erf-Screened Gaussian Decay Design

## Goal

Extend the screened Gaussian analysis so it visualizes how the Fourier-space
decay changes across several `erf` screening bandwidths while also keeping both
the hard-step screened case and the unscreened Gaussian as reference curves on
the same axes.

## Chosen Approach

Keep the work in `scripts/screened_gauss.jl` and add a second study path for
smooth screening. The new path will define an `erf`-based screening profile,
sample the screened Gaussian on a centered uniform grid, compute the continuous
FFT on that grid, and overlay `|fhat(k)|` for a user-editable list of
bandwidths on one set of axes. The figure will also include the existing
hard-step screened Gaussian as a separate reference curve and the unscreened
Gaussian as a clean baseline for the native Gaussian decay.

This is the most direct answer to the request because it makes the decay-rate
comparison visible without forcing the user to inspect separate figures or infer
the tail behavior from real and imaginary parts.

## Scope

The update will:

- keep the existing hard-step helpers available
- add an `erf` screening function with a bandwidth parameter
- add a reusable study helper that returns the sampled spectrum magnitude for a
  chosen bandwidth
- generate a single figure with one decay curve per bandwidth plus the hard-step
  and unscreened reference curves
- use a logarithmic `y` axis so the tail decay is readable

The update will not derive or plot a closed-form Fourier transform for the
`erf`-screened case. It will use the existing exact Gaussian transform for the
unscreened reference.

## Plot Design

- Horizontal axis: normalized frequency `k / (pi / dx)` to stay consistent with
  the current script
- Vertical axis: `|fhat(k)|` on a log scale
- One line per bandwidth with a legend entry showing the bandwidth
- One additional line labeled `step` for the hard-step screened Gaussian
- One additional line labeled `unscreened` for the exact Gaussian transform
- Restrict the plot to nonnegative `k` because the magnitude decay is even and
  the positive half is sufficient for comparison

## Testing

Add a focused unit test for the new screening helpers and decay-study helper so
the script logic is exercised without relying on figure generation. The test
should verify basic invariants: smaller `erf` bandwidth produces a sharper
transition in physical space, the study helper returns nonnegative spectrum
magnitudes with consistent array lengths, the step and unscreened reference
curves share the same positive-frequency axis as the `erf` curves, and the
script can be included without running `main()` as a side effect.

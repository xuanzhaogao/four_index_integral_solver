# Screened Gaussian Convolution Design

## Goal

Compare the same-line convolution of screened Gaussians with the `1 / r` kernel,
using `HCubature`, for the hard-step screening and several `erf` screening
bandwidths.

## Chosen Approach

Keep the work in `scripts/screened_gauss.jl` and add a second numerical study
for the regularized line potential

`u_delta(x0) = ∫ f(x) / |x - x0| dx`

on a finite interval, excluding the singular neighborhood
`(x0 - delta, x0 + delta)` from the integral. Each remaining subinterval will
be evaluated with `hcubature`.

This is the most practical comparison because it makes the singular treatment
explicit and applies uniformly to the step-screened and `erf`-screened cases.

## Scope

The update will:

- add a reusable `HCubature` helper for the exclusion-radius regularized
  same-line convolution
- support optional internal breakpoints so the step screen can be integrated
  piecewise across `x = 0`
- evaluate one step-screened curve and one curve per `erf` bandwidth on a shared
  target grid
- generate a figure with the regularized potentials and the difference
  `u_delta^erf - u_delta^step`
- keep the existing decay figure available

The update will not attempt a finite-part or principal-value formulation.

## Plot Design

- Top panel: `u_delta(x0)` versus `x0 / sigma`
- Bottom panel: `u_delta^erf(x0) - u_delta^step(x0)` versus `x0 / sigma`
- One dashed black line for the step-screened case
- One colored line per `erf` bandwidth
- The exclusion radius `delta` will be fixed across all curves and reported by
  the script output

## Numerical Design

- Source interval: the same finite truncation already used for the decay study
- Target grid: a smaller user-editable centered interval to keep the `HCubature`
  cost reasonable
- Step-screened case: integrate with an extra breakpoint at `x = 0`
- `erf`-screened case: integrate without extra breakpoints

## Testing

Add a focused unit test for the regularized convolution helper using a constant
source on a finite interval. For `f(x) = 1`,

`u_delta(x0) = log(((x0 - xmin) * (xmax - x0)) / delta^2)`

when `x0` is inside the interval and `delta` stays away from the endpoints.

Also add a small study-shape test that verifies the returned step and `erf`
curves share the same target axis and contain finite values.

# Design: distinct markers per line — B&W / colorblind-safe standard

Date: 2026-07-07

## Goal

Today the shared figure standard (`codes/fig_gen/fig_style.jl`) distinguishes
line/scatter series **by color only**. Series become indistinguishable on a
grayscale (B&W) printer and for colorblind readers. Add a distinct **marker**
per series (and a distinct **linestyle** where markers would clutter), locked to
the existing color palette, so every series is separable without color.

Scope: the standard file plus the figures that consume it and have color-only
line sweeps — `fig2`, `fig3`, `fig5`. `fig4` already pairs distinct color+marker
per category (compliant). `fig8` is a 3D schematic (no line sweep).

Out of scope (off-standard, separate follow-up): `fig6` and `fig7` do not
`include` `fig_style.jl` (hardcoded palettes/canvas). `fig7` already assigns a
distinct marker per `p` (already meets the goal). `fig6` panel (b) has 6 series,
which exceeds the 5-color/5-marker palette and needs a palette extension.

## Style mechanism (additions to `fig_style.jl`)

```julia
# Distinct marker per palette index, locked 1:1 to LINE_COLORS / sweep_colors.
const LINE_MARKERS = [:circle, :rect, :utriangle, :diamond, :dtriangle]
sweep_markers(n::Integer) = LINE_MARKERS[1:n]

# named marker handles mirroring QUAL (category comparisons)
const QUAL_MK = (blue = :circle, red = :rect, green = :utriangle,
                 orange = :diamond, purple = :dtriangle)

# Distinct linestyle per index — companion to sweep_colors for *dense* line
# series where a marker at every point would be noise.
const LINE_STYLES = [:solid, :dash, :dashdot, :dot]
sweep_linestyles(n::Integer) = LINE_STYLES[1:n]
```

Index-locked pairing: `#4477AA→:circle`, `#D62728→:rect`, `#228833→:utriangle`,
`#EE7733→:diamond`, `#AA3377→:dtriangle`. `:xcross` stays reserved for the
source annotation (fig2), not a data series.

## The rule (documented in the header comment)

- **Point / scatter / convergence series** → distinct color AND distinct marker,
  index-locked: `sweep_colors(i)` + `sweep_markers(i)`, or `QUAL.x` + `QUAL_MK.x`.
- **A second binary category on the same series** (adaptive/uniform, FMM/HCub,
  data/reference) → **solid vs dashed line**, not the marker. The marker keeps
  tracking the sweep.
- **Dense line-only profiles** (no per-point markers) → distinct color +
  distinct linestyle (`sweep_linestyles(i)`).

## Per-figure application

- **fig2** — marker follows `p` (`sweep_markers`) on both adaptive and uniform;
  keep solid=adaptive / dashed=uniform (already the linestyle idiom).
- **fig3(a)** — dense σ profiles: give the colored data lines distinct
  linestyles from `[:solid, :dashdot, :dot, :dashdotdot]` (skip `:dash` to avoid
  colliding with the per-contrast dashed theory guides, which are unchanged).
- **fig3(b)** — eps scatter: `marker = sweep_markers(n)[i]` per eps.
- **fig5(a,b)** — TKM curves: per-eps marker via `sweep_markers`. The black
  FMM/Direct-Sum reference keeps a distinct marker (`:star5`) + dashed line so it
  does not collide with an eps marker.
- **fig4** — already compliant; no change.

## Verification

Run each edited plotter under `--project=codes/fig_gen` with the juliaup binary;
confirm it writes its PDF without error. Eyeball each output (optionally desatured
to grayscale) to confirm series are separable by marker/linestyle alone. No
numerical values change.

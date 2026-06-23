#=
fig_style.jl — shared style for the numerical-methods figures
(fig2, fig3, fig4_near_correction, fig5, fig8).

Single source of truth for figure size, font size, and color style.
`include` this right after `using CairoMakie` in each plot script, then build
figures with the constants/helpers below. Tweak the look here once and every
figure follows.

Conventions
  * All figures share a common WIDTH (FIG_W); height is tuned per content.
  * One base fontsize (FS_BASE) set globally via set_theme!; legend labels and
    in-axis annotations use FS_LEGEND / FS_ANNOT.
  * Color style:
      - continuous fields (contourf, mesh, refinement levels) -> FIELD_CMAP (viridis)
      - every line / scattered series -> the bright qualitative LINE_COLORS
        palette: sweep_colors(n) for n-series sweeps (p, eps, contrast), and
        QUAL.* for the two-way method comparisons (fig4b standard/upsampled,
        fig8 focal/neighbor)
=#

using CairoMakie
using CairoMakie.Colors   # parse(RGBAf, "#...")

# ---- Canvas ---------------------------------------------------------------
const FIG_W       = 1000   # common width (px) for every figure
const FIG_H       = 420    # default height: two-panel line-plot figures
const FIG_H_3D    = 480    # taller height: fig8's single 3D panel
const PX_PER_UNIT = 4      # raster upscaling on save (rasterized elements)

# ---- Font sizes -----------------------------------------------------------
const FS_BASE   = 20       # base / axis-label size
const FS_LEGEND = 16       # legend labels
const FS_ANNOT  = 20       # in-axis text() annotations

set_theme!(fontsize = FS_BASE, Legend = (; labelsize = FS_LEGEND))

# ---- Colors ---------------------------------------------------------------
const FIELD_CMAP = :viridis   # continuous fields (contourf, mesh, levels)

# Bright qualitative palette for every line / scattered series: Paul Tol's
# "bright" scheme with the olive-yellow swapped for a brighter orange. Crisp
# on white and (mostly) colorblind-safe.
const LINE_COLORS = parse.(RGBAf,
    ["#4477AA",    # blue
     "#D62728",    # red    (darker, purer red so it doesn't read like orange)
     "#228833",    # green
     "#EE7733",    # orange  (replaces Tol-bright's #CCBB44 yellow)
     "#AA3377"])   # purple

# first `n` distinct colors, for an n-series sweep (p, eps, contrast, ...)
sweep_colors(n::Integer) = LINE_COLORS[1:n]

# named handles for qualitative method comparisons (fig4b, fig8)
const QUAL = (blue   = LINE_COLORS[1],
              red    = LINE_COLORS[2],
              green  = LINE_COLORS[3],
              orange = LINE_COLORS[4],
              purple = LINE_COLORS[5])

# ---- Secondary consistency knobs ------------------------------------------
const MS       = 11    # default marker size
const LW_DATA  = 2     # data-line width
const LW_GUIDE = 1.6   # guide / theory / reference dashed-line width

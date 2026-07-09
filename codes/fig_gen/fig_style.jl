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
  * Distinguishing series WITHOUT color (B&W printout / colorblind readers):
    color alone is not enough, so pair each series with a redundant marker or
    linestyle, locked 1:1 to the color index.
      - point / scatter / convergence series -> distinct color AND marker:
        sweep_colors(i) + sweep_markers(i), or QUAL.x + QUAL_MK.x.
      - a second binary category on the same series (adaptive/uniform, FMM/HCub,
        data/reference) -> solid vs dashed LINE, not the marker; the marker keeps
        tracking the sweep index.
      - dense line-only profiles (a marker at every point would be noise) ->
        distinct color + distinct linestyle: sweep_linestyles(i).
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

# ---- Markers & linestyles (color-free discriminators) ---------------------
# One distinct marker per palette index, locked 1:1 to LINE_COLORS: a series
# drawn with sweep_colors(i)[k] should use sweep_markers(i)[k], so the two
# encodings reinforce each other and the series stays separable in grayscale.
# These five shapes are the classic maximally-distinct-in-B&W set; :xcross is
# deliberately NOT here — it is reserved for the source annotation (fig2).
const LINE_MARKERS = [:circle, :rect, :utriangle, :diamond, :dtriangle]

# first `n` distinct markers, index-locked to sweep_colors(n)
sweep_markers(n::Integer) = LINE_MARKERS[1:n]

# named marker handles mirroring QUAL (category comparisons)
const QUAL_MK = (blue   = LINE_MARKERS[1],
                 red    = LINE_MARKERS[2],
                 green  = LINE_MARKERS[3],
                 orange = LINE_MARKERS[4],
                 purple = LINE_MARKERS[5])

# One distinct linestyle per index, for DENSE line-only sweeps where a marker at
# every sample would be noise (companion to sweep_colors for such series).
const LINE_STYLES = [:solid, :dash, :dashdot, :dot]

# first `n` distinct linestyles, index-locked to sweep_colors(n)
sweep_linestyles(n::Integer) = LINE_STYLES[1:n]

# ---- Secondary consistency knobs ------------------------------------------
const MS       = 11    # default marker size
const LW_DATA  = 2     # data-line width
const LW_GUIDE = 1.6   # guide / theory / reference dashed-line width

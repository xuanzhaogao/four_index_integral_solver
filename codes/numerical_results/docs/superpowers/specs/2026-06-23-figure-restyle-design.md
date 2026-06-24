# Design: restyle the Section 6 numerical-results figures to the `fig_gen` look

Date: 2026-06-23

## Goal

Regenerate the publication figures under `codes/numerical_results/exp*/` so they
adopt the **plotting style** defined by `codes/fig_gen/fig_style.jl` — the shared
Tol-bright palette, common canvas width, font sizes, line widths, marker size,
and raster upscaling used by the paper's other figures (fig2–fig8). This is a
restyle + light content edit of the *existing* figures; it is not a data re-run.

Out of scope: regenerating any underlying `.jls`/`.csv` data (all required data
already exists on disk), and experiments **6.2 (ablation)** and **6.4
(scaling/amortization)**, which have only `run_*.jl` and no data/plots yet.

## Style mechanism

Each plotter `include`s the canonical style file directly (single source of
truth, no vendored copy):

```julia
include(joinpath(@__DIR__, "..", "..", "..", "fig_gen", "fig_style.jl"))
```

All plot scripts live at `exp*/scripts/`, so the relative path is identical for
every one. `fig_style.jl` only needs `CairoMakie` (+ `CairoMakie.Colors`), which
the `numerical_results` project already provides, so it runs under
`--project=codes/numerical_results` unchanged. Include it immediately after the
script's own `using CairoMakie`, per the file's documented convention. It calls
`set_theme!` globally and exports: `FIG_W`, `FIG_H`, `PX_PER_UNIT`, `FS_BASE`,
`FS_LEGEND`, `FS_ANNOT`, `FIELD_CMAP`, `LINE_COLORS`, `sweep_colors(n)`, `QUAL`,
`MS`, `LW_DATA`, `LW_GUIDE`.

Common conversions applied to every line/scatter figure:
- canvas `Figure(size = (FIG_W, FIG_H))` (drop the ad-hoc `820×330` / `1400×400`);
- replace the hardcoded Okabe–Ito hex palette with `sweep_colors(n)` (n-series
  sweeps) or `QUAL.*` (category comparisons);
- data series at `LW_DATA` + `markersize = MS`; guide/reference/theory dashed
  lines at `LW_GUIDE`;
- drop per-call `fontsize=`/`labelsize=` overrides — inherit `FS_BASE`/`FS_LEGEND`
  from the theme; in-axis annotations use `FS_ANNOT`;
- save PDF with `px_per_unit = PX_PER_UNIT`; keep an accompanying PNG where the
  script already wrote one.

## Figure set (5 figures)

### 1. exp61 slab convergence — `exp61_convergence/scripts/plot_slab.jl`
- `sweep_colors(3)` for `p = 2,4,6` solid curves (circles, `LW_DATA`).
- Keep the archived **no-edge-correction** contrast: same color per `p`, dashed
  `:utriangle`, `LW_GUIDE` (matches `fig_gen`'s solid-data / dashed-guide idiom).
  Keep the gray legend proxy entry "no edge corr.".
- **Remove** the `ε = 10⁻⁴` `hlines!` + its `text!` label.
- Switch data loader `Harness → Lite` (records are scalar `NamedTuple`s; avoids
  loading `BoundaryIntegral`). `using .Lite` after `include("common/Lite.jl")`.

### 2. exp61 Fig.-1 (slab + 2 cubes) convergence — `exp61_convergence/scripts/plot_fig1.jl`
- `sweep_colors(3)` for `p = 2,4,6` (circles, `LW_DATA`); panels (a) V-error vs
  DOF, (b) GMRES iterations vs DOF.
- **Remove** the `ε = 10⁻⁴` `hlines!` + its `text!` label.
- Switch data loader `Harness → Lite`.

### 3. exp61 slab + 2 cubes geometry — `exp61_convergence/scripts/visualize_fig1_system.jl`
- 3D schematic; not a line plot, so only partial styling: inherit `FS_BASE`
  fonts; pull the legend's interface-category colors from `QUAL.*` and the panel
  refinement-level coloring from `FIELD_CMAP`.
- Keeps its `BoundaryIntegral` dependency (it builds/loads the actual interface
  to draw panels).

### 4. exp63 contrast — `exp63_contrast/scripts/plot_contrast.jl`
- `sweep_colors(4)` for `ε₂ ∈ {6, 20, 60, 200}` (keep the `γ₁₂` legend labels).
- **Remove** the `ε = 10⁻⁴` `hlines!` + its `text!` label.
- Already uses `Lite`; just apply canvas/palette/fonts.

### 5. exp66 multicube scaling — `exp66_multicube/scripts/plot_scaling.jl`
- Reduce from 3 panels to **2**: **drop the peak-RAM panel**; keep
  (a) runtime vs `K` and (b) GMRES-iterations vs `K`.
- Panel (a): map the categories to `QUAL.*` — precompute (dashed), block solve,
  eval, total (`:rect` marker). Add a **faint dashed reference line**
  `naive(K) = K · total(K=1)` (the cost of repeating a single-RHS solve K times),
  and **text-label the multi-RHS speedup** `naive(K)/total(K)` at the points:
  `K=4 → 1.9×`, `K=10 → 3.0×`, `K=13 → 3.2×`, `K=19 → 3.4×`. Speedup is computed
  from the gathered records / `data/multicube.csv` (no new data needed).
- Panel (b): GMRES iterations vs `K`, single series.
- This script already has no `BoundaryIntegral` dependency (plain `Serialization`).

### Dropped (not regenerated)
- exp61 slab geometry (`visualize_slab_system.jl` output) — redundant; the slab
  is fully described in one sentence and appears as the top plate of figure 3.
- exp65 orbital multi-RHS scaling (`fig65_multirhs_scaling`) — near-duplicate of
  figure 5; benchmark/campaign artifact, not a paper figure.
- exp66 system geometry (`fig66_system`) — same topology as figure 3.
- exp62 / exp64 — no data yet; deferred to a separate effort.

## Deferred (not in this pass, flagged for later)
- exp63 V/V_vacuum screening-saturation curve (prompt 6.3.2): `v_vacuum.jls` +
  `contrast_ratio.csv` exist; could become an added panel later.
- Content question (not styling): the convergence sweeps are `p ∈ {2,4,6}` at
  `ε = 10⁻⁴`, whereas `prompt.md`/`PLAN.md` specify `p ∈ {4,6,8}` down to
  `10⁻¹¹`. If the paper needs the prompt's grid, the data must be re-run — out of
  scope here.

## Verification

For each regenerated figure: run its plotter under
`--project=codes/numerical_results` (juliaup binary, no heavy BI for figs 1/2/4/5),
confirm it writes the PDF/PNG without error, and visually compare the output
against the pre-restyle version (same data, new look; the dropped `ε` lines and
RAM panel gone; exp66 speedup labels present). No numerical values change.

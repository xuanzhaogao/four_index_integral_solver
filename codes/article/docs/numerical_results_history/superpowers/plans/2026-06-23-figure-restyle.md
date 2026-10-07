# Numerical-Results Figure Restyle Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restyle the 5 kept Section-6 figures under `codes/numerical_results/exp*/` to the shared `fig_gen` plotting style, with the per-figure content edits agreed in the spec.

**Architecture:** Each plotter `include`s `codes/fig_gen/fig_style.jl` directly (single source of truth) and replaces its ad-hoc canvas/palette/fonts with the shared `FIG_W`/`FS_*`/`sweep_colors`/`QUAL`/`LW_*`/`MS`/`PX_PER_UNIT`. Convergence plotters also switch their data loader from `Harness` (loads `BoundaryIntegral`) to the BI-free `common/Lite.jl`. No underlying data is regenerated — all `.jls`/`.csv` already exist on disk.

**Tech Stack:** Julia 1.12 (juliaup), CairoMakie, the `numerical_results` Julia project.

**Spec:** `docs/superpowers/specs/2026-06-23-figure-restyle-design.md`

## Global Constraints

- Julia binary (juliaup-managed; never `module load julia`): `/mnt/home/xgao1/.juliaup/bin/julia`
- Project root (`NR`): `/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results`
- Run pattern (absolute paths, no `cd`): `/mnt/home/xgao1/.juliaup/bin/julia --project=$NR $NR/<script>`
- Shared style file, included from every `exp*/scripts/` plotter via:
  `include(joinpath(@__DIR__, "..", "..", "..", "fig_gen", "fig_style.jl"))` — it `using CairoMakie`s internally, so the line must come **after** the script's own `using CairoMakie`. It exports `FIG_W=1000`, `FIG_H=420`, `PX_PER_UNIT=4`, `FS_BASE=20`, `FS_LEGEND=16`, `FS_ANNOT=20`, `FIELD_CMAP`, `sweep_colors(n)`, `QUAL.{blue,red,green,orange,purple}`, `MS=11`, `LW_DATA=2`, `LW_GUIDE=1.6`, and calls `set_theme!(fontsize=FS_BASE, Legend=(;labelsize=FS_LEGEND))` globally.
- Do **not** regenerate data. Do **not** touch exp62, exp64, exp65, or the dropped geometry/system figures.
- Keep each plotter's panel **titles** verbatim (incl. their "ε = 10⁻⁴" wording); only the dotted ε *line* + its in-axis `text!` label are removed.
- Every plotter saves both a PDF (at `px_per_unit = PX_PER_UNIT`) and a review PNG (`px_per_unit = 2`).
- Verification is run-and-look (no unit tests for figures): the script must exit 0 and write its files; then visually compare the PNG against the pre-restyle look.
- Commits: the repo default branch is `main`; all work lands on the branch created in Task 0. Commit per task.

---

### Task 0: Create work branch

**Files:** none (git only)

- [ ] **Step 1: Create and switch to the feature branch**

```bash
git -C /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results checkout -b figure-restyle
```

- [ ] **Step 2: Confirm the style file resolves from a plotter dir**

Run:
```bash
ls /mnt/home/xgao1/work/four_index_integral_solver/codes/fig_gen/fig_style.jl
```
Expected: the path prints (file exists). This is the target of every `../../../fig_gen/fig_style.jl` include.

---

### Task 1: exp61 slab convergence — `exp61_convergence/scripts/plot_slab.jl`

**Files:**
- Modify: `exp61_convergence/scripts/plot_slab.jl`

**Interfaces:**
- Consumes: `Lite.load_ref(path)::NamedTuple` (fields `.N`, `.V`, `.niter`); `sweep_colors`, `QUAL`, `FIG_W`, `FIG_H`, `LW_DATA`, `LW_GUIDE`, `MS`, `PX_PER_UNIT` from `fig_style.jl`.
- Produces: `figs/fig61_slab_convergence.{pdf,png}` (no API consumed downstream).

- [ ] **Step 1: Swap the loader to Lite and add the style include**

Replace lines 5–7:
```julia
include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
using CairoMakie, LaTeXStrings, LinearAlgebra
```
with:
```julia
include(joinpath(@__DIR__, "..", "..", "common", "Lite.jl"))
using .Lite
using CairoMakie, LaTeXStrings, LinearAlgebra
include(joinpath(@__DIR__, "..", "..", "..", "fig_gen", "fig_style.jl"))
```

- [ ] **Step 2: Replace the palette and canvas**

Replace line 19:
```julia
const COL = Dict(2 => "#0072B2", 4 => "#D55E00", 6 => "#009E73")
```
with:
```julia
const _SC = sweep_colors(3)
const COL = Dict(2 => _SC[1], 4 => _SC[2], 6 => _SC[3])
```

Replace line 21:
```julia
fig = Figure(size = (820, 330))
```
with:
```julia
fig = Figure(size = (FIG_W, FIG_H))
```

- [ ] **Step 3: Apply data/guide line weights to the four `scatterlines!` calls**

Replace the solid-series block (lines 29–36):
```julia
for p in (2, 4, 6)
    s = load_series(joinpath(DATA, "raw"), p)
    Ns = [d.N for d in s]
    ev = [abs(d.V - ref.V) / abs(ref.V) for d in s]
    its = [d.niter for d in s]
    scatterlines!(ax1, Ns, ev; color = COL[p], marker = :circle, label = L"p = %$p")
    scatterlines!(ax2, Ns, its; color = COL[p], marker = :circle, label = L"p = %$p")
end
```
with:
```julia
for p in (2, 4, 6)
    s = load_series(joinpath(DATA, "raw"), p)
    Ns = [d.N for d in s]
    ev = [abs(d.V - ref.V) / abs(ref.V) for d in s]
    its = [d.niter for d in s]
    scatterlines!(ax1, Ns, ev; color = COL[p], marker = :circle,
        linewidth = LW_DATA, markersize = MS, label = L"p = %$p")
    scatterlines!(ax2, Ns, its; color = COL[p], marker = :circle,
        linewidth = LW_DATA, markersize = MS, label = L"p = %$p")
end
```

Replace the no-edge dashed block (lines 39–46):
```julia
for p in (2, 4, 6)
    s = load_series_uc(joinpath(DATA, "raw", "noedges_run1"), p)
    Ns = [d.N for d in s]
    ev = [abs(d.V - ref.V) / abs(ref.V) for d in s]
    its = [d.niter for d in s]
    scatterlines!(ax1, Ns, ev; color = (COL[p], 0.55), marker = :utriangle, linestyle = :dash)
    scatterlines!(ax2, Ns, its; color = (COL[p], 0.55), marker = :utriangle, linestyle = :dash)
end
```
with:
```julia
for p in (2, 4, 6)
    s = load_series_uc(joinpath(DATA, "raw", "noedges_run1"), p)
    Ns = [d.N for d in s]
    ev = [abs(d.V - ref.V) / abs(ref.V) for d in s]
    its = [d.niter for d in s]
    scatterlines!(ax1, Ns, ev; color = (COL[p], 0.55), marker = :utriangle,
        markersize = MS, linewidth = LW_GUIDE, linestyle = :dash)
    scatterlines!(ax2, Ns, its; color = (COL[p], 0.55), marker = :utriangle,
        markersize = MS, linewidth = LW_GUIDE, linestyle = :dash)
end
```

Replace the legend-proxy call (lines 48–49):
```julia
scatterlines!(ax1, [NaN], [NaN]; color = :gray40, marker = :utriangle,
    linestyle = :dash, label = "no edge corr.")
```
with:
```julia
scatterlines!(ax1, [NaN], [NaN]; color = :gray40, marker = :utriangle,
    linewidth = LW_GUIDE, linestyle = :dash, label = "no edge corr.")
```

- [ ] **Step 4: Remove the ε line + label and drop the legend `labelsize` overrides**

Replace lines 51–55:
```julia
hlines!(ax1, [1e-4]; color = :black, linestyle = :dot)
text!(ax1, 3.0e5, 1.22e-4; text = L"\varepsilon = 10^{-4}", fontsize = 12)
axislegend(ax1; position = :lb, framevisible = false, labelsize = 11)
ylims!(ax2, 0, 28)
axislegend(ax2; position = :lt, framevisible = false, labelsize = 11)
```
with:
```julia
axislegend(ax1; position = :lb, framevisible = false)
ylims!(ax2, 0, 28)
axislegend(ax2; position = :lt, framevisible = false)
```

- [ ] **Step 5: Save PDF at the shared raster + add a review PNG**

Replace line 57:
```julia
save(joinpath(FIGS, "fig61_slab_convergence.pdf"), fig)
```
with:
```julia
save(joinpath(FIGS, "fig61_slab_convergence.pdf"), fig; px_per_unit = PX_PER_UNIT)
save(joinpath(FIGS, "fig61_slab_convergence.png"), fig; px_per_unit = 2)
```

- [ ] **Step 6: Run and verify output is written**

Run:
```bash
/mnt/home/xgao1/.juliaup/bin/julia --project=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp61_convergence/scripts/plot_slab.jl
```
Expected: prints `wrote figs/fig61_slab_convergence.pdf` and exits 0; `figs/fig61_slab_convergence.png` exists and is newer than before.

- [ ] **Step 7: Visual check**

Open `exp61_convergence/figs/fig61_slab_convergence.png`. Confirm: Tol-bright blue/red/green solid `p` curves with circle markers; same-colored dashed `utriangle` "no edge corr." curves; no dotted ε line / no "ε = 10⁻⁴" floating label; larger fonts; both panel legends present; titles fit within their panels.

- [ ] **Step 8: Commit**

```bash
git -C /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results add exp61_convergence/scripts/plot_slab.jl exp61_convergence/figs/fig61_slab_convergence.pdf exp61_convergence/figs/fig61_slab_convergence.png
git -C /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results commit -m "fig61: restyle slab convergence to fig_gen look (Lite loader, drop eps line)"
```

---

### Task 2: exp61 Fig.-1 convergence — `exp61_convergence/scripts/plot_fig1.jl`

**Files:**
- Modify: `exp61_convergence/scripts/plot_fig1.jl`

**Interfaces:**
- Consumes: `Lite.load_ref`; `sweep_colors`, `FIG_W`, `FIG_H`, `LW_DATA`, `MS`, `PX_PER_UNIT`.
- Produces: `figs/fig61_fig1_convergence.{pdf,png}`.

- [ ] **Step 1: Swap the loader to Lite and add the style include**

Replace lines 4–6:
```julia
include(joinpath(@__DIR__, "..", "..", "common", "Harness.jl"))
using .Harness
using CairoMakie, LaTeXStrings, LinearAlgebra
```
with:
```julia
include(joinpath(@__DIR__, "..", "..", "common", "Lite.jl"))
using .Lite
using CairoMakie, LaTeXStrings, LinearAlgebra
include(joinpath(@__DIR__, "..", "..", "..", "fig_gen", "fig_style.jl"))
```

- [ ] **Step 2: Replace the palette and canvas**

Replace line 15:
```julia
const COL = Dict(2 => "#0072B2", 4 => "#D55E00", 6 => "#009E73")
```
with:
```julia
const _SC = sweep_colors(3)
const COL = Dict(2 => _SC[1], 4 => _SC[2], 6 => _SC[3])
```

Replace line 17:
```julia
fig = Figure(size = (820, 330))
```
with:
```julia
fig = Figure(size = (FIG_W, FIG_H))
```

- [ ] **Step 3: Apply line weights to the two `scatterlines!` calls**

Replace lines 25–31:
```julia
for p in (2, 4, 6)
    s = load_series(p)
    Ns = [d.N for d in s]
    ev = [abs(d.V - ref.V) / abs(ref.V) for d in s]
    scatterlines!(ax1, Ns, ev; color = COL[p], marker = :circle, label = L"p = %$p")
    scatterlines!(ax2, Ns, [d.niter for d in s]; color = COL[p], marker = :circle, label = L"p = %$p")
end
```
with:
```julia
for p in (2, 4, 6)
    s = load_series(p)
    Ns = [d.N for d in s]
    ev = [abs(d.V - ref.V) / abs(ref.V) for d in s]
    scatterlines!(ax1, Ns, ev; color = COL[p], marker = :circle,
        linewidth = LW_DATA, markersize = MS, label = L"p = %$p")
    scatterlines!(ax2, Ns, [d.niter for d in s]; color = COL[p], marker = :circle,
        linewidth = LW_DATA, markersize = MS, label = L"p = %$p")
end
```

- [ ] **Step 4: Remove the ε line + label and drop the legend `labelsize` overrides**

Replace lines 32–36:
```julia
hlines!(ax1, [1e-4]; color = :black, linestyle = :dot)
text!(ax1, 1.0e6, 1.15e-4; text = L"\varepsilon = 10^{-4}", fontsize = 12)
axislegend(ax1; position = :lb, framevisible = false, labelsize = 11)
ylims!(ax2, 0, 28)
axislegend(ax2; position = :lt, framevisible = false, labelsize = 11)
```
with:
```julia
axislegend(ax1; position = :lb, framevisible = false)
ylims!(ax2, 0, 28)
axislegend(ax2; position = :lt, framevisible = false)
```

- [ ] **Step 5: Save PDF at the shared raster + add a review PNG**

Replace line 38:
```julia
save(joinpath(FIGS, "fig61_fig1_convergence.pdf"), fig)
```
with:
```julia
save(joinpath(FIGS, "fig61_fig1_convergence.pdf"), fig; px_per_unit = PX_PER_UNIT)
save(joinpath(FIGS, "fig61_fig1_convergence.png"), fig; px_per_unit = 2)
```

- [ ] **Step 6: Run and verify output is written**

Run:
```bash
/mnt/home/xgao1/.juliaup/bin/julia --project=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp61_convergence/scripts/plot_fig1.jl
```
Expected: prints `wrote figs/fig61_fig1_convergence.pdf`, exits 0; the PNG exists.

- [ ] **Step 7: Visual check**

Open `exp61_convergence/figs/fig61_fig1_convergence.png`. Confirm Tol-bright `p=2,4,6` curves, no ε dotted line/label, shared fonts/canvas, both legends present.

- [ ] **Step 8: Commit**

```bash
git -C /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results add exp61_convergence/scripts/plot_fig1.jl exp61_convergence/figs/fig61_fig1_convergence.pdf exp61_convergence/figs/fig61_fig1_convergence.png
git -C /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results commit -m "fig61: restyle Fig.-1 convergence to fig_gen look (Lite loader, drop eps line)"
```

---

### Task 3: exp63 contrast — `exp63_contrast/scripts/plot_contrast.jl`

**Files:**
- Modify: `exp63_contrast/scripts/plot_contrast.jl`

**Interfaces:**
- Consumes: `Lite.load_ref` (already used); `sweep_colors`, `FIG_W`, `FIG_H`, `LW_DATA`, `MS`, `PX_PER_UNIT`.
- Produces: `figs/fig63_contrast.{pdf,png}`.

- [ ] **Step 1: Add the style include**

After line 7 (`using CairoMakie, LaTeXStrings`), insert:
```julia
include(joinpath(@__DIR__, "..", "..", "..", "fig_gen", "fig_style.jl"))
```

- [ ] **Step 2: Replace the hardcoded palette**

Replace lines 20–22:
```julia
# Okabe-Ito, one color per eps2
const COLS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"]
col(i) = COLS[mod1(i, length(COLS))]
```
with:
```julia
# Tol-bright sweep palette, one color per eps2
col(i) = sweep_colors(5)[mod1(i, 5)]
```

- [ ] **Step 3: Replace the canvas**

Replace line 24:
```julia
fig = Figure(size = (820, 330))
```
with:
```julia
fig = Figure(size = (FIG_W, FIG_H))
```

- [ ] **Step 4: Apply line weights to the two `scatterlines!` calls**

Replace lines 39–40:
```julia
    scatterlines!(ax1, Ns, ev; color = col(i), marker = :circle, label = lab)
    scatterlines!(ax2, Ns, [t.niter for t in ts]; color = col(i), marker = :circle)
```
with:
```julia
    scatterlines!(ax1, Ns, ev; color = col(i), marker = :circle,
        linewidth = LW_DATA, markersize = MS, label = lab)
    scatterlines!(ax2, Ns, [t.niter for t in ts]; color = col(i), marker = :circle,
        linewidth = LW_DATA, markersize = MS)
```

- [ ] **Step 5: Remove the ε line + label and drop the legend `labelsize`**

Replace lines 42–44:
```julia
hlines!(ax1, [1e-4]; color = :black, linestyle = :dot)
text!(ax1, 2.0e5, 1.2e-4; text = L"\varepsilon = 10^{-4}", fontsize = 12)
axislegend(ax1; position = :rt, framevisible = false, labelsize = 10)
```
with:
```julia
axislegend(ax1; position = :rt, framevisible = false)
```

- [ ] **Step 6: Save PDF at the shared raster**

Replace lines 48–49:
```julia
save(joinpath(FIGS, "fig63_contrast.pdf"), fig)
save(joinpath(FIGS, "fig63_contrast.png"), fig; px_per_unit = 2)
```
with:
```julia
save(joinpath(FIGS, "fig63_contrast.pdf"), fig; px_per_unit = PX_PER_UNIT)
save(joinpath(FIGS, "fig63_contrast.png"), fig; px_per_unit = 2)
```

- [ ] **Step 7: Run and verify output is written**

Run:
```bash
/mnt/home/xgao1/.juliaup/bin/julia --project=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp63_contrast/scripts/plot_contrast.jl
```
Expected: prints `wrote figs/fig63_contrast.{pdf,png}`, exits 0.

- [ ] **Step 8: Visual check**

Open `exp63_contrast/figs/fig63_contrast.png`. Confirm 4 Tol-bright `ε₂` curves with `γ₁₂` legend labels, no ε dotted line/label, shared fonts/canvas.

- [ ] **Step 9: Commit**

```bash
git -C /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results add exp63_contrast/scripts/plot_contrast.jl exp63_contrast/figs/fig63_contrast.pdf exp63_contrast/figs/fig63_contrast.png
git -C /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results commit -m "fig63: restyle contrast figure to fig_gen look (drop eps line)"
```

---

### Task 4: exp66 multicube scaling — `exp66_multicube/scripts/plot_scaling.jl`

**Files:**
- Modify: `exp66_multicube/scripts/plot_scaling.jl`

**Interfaces:**
- Consumes: gathered records' `t_precompute`, `t_pottrg`, `t_solve_block`, `t_eval`, `niter`, `K` (already read into `t_pre`,`t_pot`,`t_blk`,`t_evl`,`niter`,`K`,`t_tot`); `QUAL`, `FIG_W`, `FIG_H`, `FS_ANNOT`, `LW_DATA`, `LW_GUIDE`, `MS`, `PX_PER_UNIT`.
- Produces: `figs/fig66_multicube_scaling.{pdf,png}` (2-panel) + `data/multicube.csv` (unchanged gather).

- [ ] **Step 1: Add the style include**

After line 12 (`using CairoMakie, LaTeXStrings, Serialization`), insert:
```julia
include(joinpath(@__DIR__, "..", "..", "..", "fig_gen", "fig_style.jl"))
```

- [ ] **Step 2: Drop the Okabe-Ito palette constant**

Delete lines 54–55:
```julia
const C = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00"]  # Okabe-Ito
xlab = L"K"
```
and replace with:
```julia
xlab = L"K"
# multi-RHS speedup vs repeating a single-RHS solve K times (precompute is ~3% of total)
naive = K .* t_tot[1]
speedup = naive ./ t_tot
```

- [ ] **Step 3: Rebuild the figure block as 2 panels with the speedup annotation**

Replace the whole figure block (lines 57–83):
```julia
begin
    fig = Figure(size = (1400, 400), fontsize = 18)

    ax1 = Axis(fig[1, 1]; xlabel = xlab, ylabel = "runtime (s)", xticks = K, xscale = log10, yscale = log10)
    # reference slopes (faint guides)
    lines!(ax1, K, fill(sum(t_pre) / length(t_pre), length(K)); color = :gray80)              # slope 0
    lines!(ax1, K, t_evl[1] .* K;                                color = :gray80, linestyle = :dot)  # slope 1 (∝K)
    # fixed cost as its average level (pottrg omitted — negligible here)
    hlines!(ax1, [sum(t_pre) / length(t_pre)]; color = C[1], linestyle = :dash, linewidth = 2, label = "precompute")
    # scaling costs vs K
    scatterlines!(ax1, K, t_blk; color = C[2], marker = :circle, label = "block solve")
    scatterlines!(ax1, K, t_evl; color = C[3], marker = :circle, label = "eval")
    scatterlines!(ax1, K, t_tot; color = :black, marker = :rect, linewidth = 2, label = "total")
    axislegend(ax1; position = :rb, labelsize = 16)
    ylims!(ax1, 1.0, 10^(3.3))

    ax2 = Axis(fig[1, 2]; xlabel = xlab, ylabel = "peak RAM (GB)", xticks = K)
    scatterlines!(ax2, K, rss; color = C[4], marker = :circle)
    length(K) > 1 && hlines!(ax2, [rss[1]]; color = :gray70, linestyle = :dot)
    ylims!(ax2, 0, 350)

    ax3 = Axis(fig[1, 3]; xlabel = xlab, ylabel = "GMRES iterations", xticks = K)
    scatterlines!(ax3, K, niter; color = C[5], marker = :circle)
    ylims!(ax3, 0, 40)

    fig
end
```
with:
```julia
begin
    fig = Figure(size = (FIG_W, FIG_H))

    # ----- Panel (a): runtime decomposition + multi-RHS speedup -----------
    ax1 = Axis(fig[1, 1]; xlabel = xlab, ylabel = "runtime (s)", xticks = K,
        xscale = log10, yscale = log10, title = "(a) multi-RHS runtime")
    # naive baseline: repeat a single-RHS solve K times (the gap to total = speedup)
    lines!(ax1, K, naive; color = (:gray50, 0.9), linestyle = :dash,
        linewidth = LW_GUIDE, label = L"K \times \mathrm{single\text{-}RHS}")
    # fixed cost as its average level (pottrg omitted — negligible here)
    hlines!(ax1, [sum(t_pre) / length(t_pre)]; color = QUAL.blue, linestyle = :dash,
        linewidth = LW_GUIDE, label = "precompute")
    scatterlines!(ax1, K, t_blk; color = QUAL.orange, marker = :circle,
        linewidth = LW_DATA, markersize = MS, label = "block solve")
    scatterlines!(ax1, K, t_evl; color = QUAL.green, marker = :circle,
        linewidth = LW_DATA, markersize = MS, label = "eval")
    scatterlines!(ax1, K, t_tot; color = :black, marker = :rect,
        linewidth = LW_DATA, markersize = MS, label = "total")
    # speedup labels at each K >= 4 (skip the K=1 baseline)
    for i in eachindex(K)
        K[i] == 1 && continue
        text!(ax1, K[i], t_tot[i]; text = string(round(speedup[i]; digits = 1), "×"),
            align = (:center, :top), offset = (0, -6), fontsize = FS_ANNOT - 4,
            color = :black)
    end
    axislegend(ax1; position = :rb)
    ylims!(ax1, 1.0, 10^(3.7))

    # ----- Panel (b): GMRES iterations vs K -------------------------------
    ax2 = Axis(fig[1, 2]; xlabel = xlab, ylabel = "GMRES iterations", xticks = K,
        title = "(b) GMRES iterations")
    scatterlines!(ax2, K, niter; color = QUAL.purple, marker = :circle,
        linewidth = LW_DATA, markersize = MS)
    ylims!(ax2, 0, 40)

    fig
end
```

- [ ] **Step 4: Save PDF at the shared raster**

Replace lines 85–86:
```julia
save(joinpath(FIGS, FIGNAME * ".pdf"), fig)
save(joinpath(FIGS, FIGNAME * ".png"), fig)
```
with:
```julia
save(joinpath(FIGS, FIGNAME * ".pdf"), fig; px_per_unit = PX_PER_UNIT)
save(joinpath(FIGS, FIGNAME * ".png"), fig; px_per_unit = 2)
```

- [ ] **Step 5: Run and verify output is written**

Run:
```bash
/mnt/home/xgao1/.juliaup/bin/julia --project=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp66_multicube/scripts/plot_scaling.jl
```
Expected: prints `gathered 5 cutoffs (K = [1, 4, 10, 13, 19]) -> ...`, exits 0; 2-panel PDF/PNG written.

- [ ] **Step 6: Visual check**

Open `exp66_multicube/figs/fig66_multicube_scaling.png`. Confirm: exactly **2 panels** (no RAM panel); panel (a) has precompute (dashed)/block solve/eval/total + a faint dashed `K × single-RHS` line above `total`, with `1.9×`, `3.0×`, `3.2×`, `3.4×` labels at K=4,10,13,19; panel (b) is GMRES iterations vs K; Tol-bright/QUAL colors, shared fonts.

- [ ] **Step 7: Commit**

```bash
git -C /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results add exp66_multicube/scripts/plot_scaling.jl exp66_multicube/figs/fig66_multicube_scaling.pdf exp66_multicube/figs/fig66_multicube_scaling.png exp66_multicube/data/multicube.csv
git -C /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results commit -m "fig66: 2-panel multicube scaling, drop RAM, add multi-RHS speedup labels, fig_gen look"
```

---

### Task 5: exp61 slab+2cubes geometry — `exp61_convergence/scripts/visualize_fig1_system.jl`

**Files:**
- Modify: `exp61_convergence/scripts/visualize_fig1_system.jl`

**Interfaces:**
- Consumes: `Harness` + `BoundaryIntegral` (kept — builds the actual interface to draw panels); `QUAL`, `FIG_W` from `fig_style.jl`.
- Produces: `figs/fig61_system_fig1_geometry.{png,pdf}`.

Note: this is a 3D schematic, so styling is limited to fonts (inherited via `set_theme!`) and the interface-category legend colors. The `BoundaryIntegral` dependency stays; this run does a modest RHS-adaptive interface build (no GMRES solve) and precompiles BI, so allow a couple of minutes.

- [ ] **Step 1: Add the style include**

After line 10 (`using CairoMakie, Printf`), insert:
```julia
include(joinpath(@__DIR__, "..", "..", "..", "fig_gen", "fig_style.jl"))
```

- [ ] **Step 2: Replace the category color map with QUAL + neutrals**

Replace lines 44–51:
```julia
COLOR = Dict(
    (1.0, 4.0) => "#E8C547",
    (1.0, 12.0) => "#9aa0a6",
    (4.0, 12.0) => "#CC79A7",
    (1.0, 10.0) => "#56B4E9",
    (4.0, 10.0) => "#D55E00",
    (10.0, 12.0) => "#009E73",
)
```
with:
```julia
# bright Tol palette for the interesting interfaces; neutral grays for the two
# outer cube|vac faces; source stays red (so QUAL.red is left unused here).
COLOR = Dict(
    (1.0, 4.0)   => RGBAf(0.80, 0.80, 0.82, 1.0),   # Ω₁ | vac   (neutral)
    (1.0, 12.0)  => RGBAf(0.62, 0.62, 0.66, 1.0),   # Ω₂ | vac   (neutral)
    (4.0, 12.0)  => QUAL.purple,                     # Ω₁ | Ω₂ shared
    (1.0, 10.0)  => QUAL.blue,                        # slab | vac
    (4.0, 10.0)  => QUAL.orange,                      # slab | Ω₁ contact
    (10.0, 12.0) => QUAL.green,                       # slab | Ω₂ contact
)
```

- [ ] **Step 3: Use the shared canvas width and drop the legend `labelsize`**

Replace line 76:
```julia
fig = Figure(size = (1000, 620))
```
with:
```julia
fig = Figure(size = (FIG_W, 620))
```

Replace line 93:
```julia
axislegend(ax; position = :rt, framevisible = false, labelsize = 11)
```
with:
```julia
axislegend(ax; position = :rt, framevisible = false)
```

- [ ] **Step 4: Save PDF at the shared raster**

Replace lines 112–113:
```julia
save(joinpath(FIGS, "fig61_system_fig1_geometry.png"), fig; px_per_unit = 2)
save(joinpath(FIGS, "fig61_system_fig1_geometry.pdf"), fig)
```
with:
```julia
save(joinpath(FIGS, "fig61_system_fig1_geometry.png"), fig; px_per_unit = 2)
save(joinpath(FIGS, "fig61_system_fig1_geometry.pdf"), fig; px_per_unit = PX_PER_UNIT)
```

- [ ] **Step 5: Run and verify output is written**

Run:
```bash
/mnt/home/xgao1/.juliaup/bin/julia --project=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp61_convergence/scripts/visualize_fig1_system.jl
```
Expected: prints the panel/point counts and `wrote figs/fig61_system_fig1_geometry.{png,pdf}`, exits 0 (BI precompiles on first run — allow a few minutes).

- [ ] **Step 6: Visual check**

Open `exp61_convergence/figs/fig61_system_fig1_geometry.png`. Confirm: cubes' outer `|vac` faces are neutral gray; the slab-top, two contacts, and shared face are blue/orange/green/purple; source dot still red; legend readable at the larger inherited font; the slab-top RHS-adaptive panel inset still renders.

- [ ] **Step 7: Commit**

```bash
git -C /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results add exp61_convergence/scripts/visualize_fig1_system.jl exp61_convergence/figs/fig61_system_fig1_geometry.png exp61_convergence/figs/fig61_system_fig1_geometry.pdf
git -C /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results commit -m "fig61: restyle slab+2cubes geometry (QUAL palette, shared fonts)"
```

---

### Task 6: Cleanup — remove dropped figure outputs

**Files:**
- Delete: `exp61_convergence/figs/fig61_system_slab_geometry.{png,pdf}`, `exp61_convergence/figs/fig61_slab_convergence` stray previews are kept; `exp65_orbital_bench/figs/fig65_multirhs_scaling.pdf`; `exp66_multicube/figs/fig66_system.{png,pdf}`

Rationale: the spec drops these figures. Their generating scripts are left in place (data/analysis history), but the stale output files should not ship as paper figures.

- [ ] **Step 1: Confirm what will be removed**

Run:
```bash
ls -l /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp61_convergence/figs/fig61_system_slab_geometry.* \
      /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/figs/fig65_multirhs_scaling.pdf \
      /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp66_multicube/figs/fig66_system.*
```
Expected: the dropped output files list (and only these).

- [ ] **Step 2: Remove the dropped outputs**

```bash
git -C /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results rm \
  exp61_convergence/figs/fig61_system_slab_geometry.png exp61_convergence/figs/fig61_system_slab_geometry.pdf \
  exp65_orbital_bench/figs/fig65_multirhs_scaling.pdf \
  exp66_multicube/figs/fig66_system.png exp66_multicube/figs/fig66_system.pdf
```

- [ ] **Step 3: Commit**

```bash
git -C /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results commit -m "figs: drop slab-geometry, exp65 scaling, and exp66 system outputs per restyle spec"
```

---

## Self-Review

- **Spec coverage:** Fig 1 slab convergence → Task 1; Fig 2 Fig.-1 convergence → Task 2; Fig 4 contrast → Task 3; Fig 5 multicube scaling (2-panel, RAM dropped, GMRES kept, speedup labels) → Task 4; Fig 3 geometry → Task 5; dropped outputs → Task 6. Style mechanism (direct include) and `Harness→Lite` are in Tasks 1–2 and constraints. All spec items covered.
- **Placeholder scan:** every code step shows the exact old→new text; commands have expected output; no TBD/TODO.
- **Type consistency:** `Lite.load_ref` returns the same scalar `NamedTuple` (`.N`,`.V`,`.niter`) that `Harness.load_ref` did (both are `Serialization.deserialize`); `sweep_colors`, `QUAL.*`, `FIG_W`, `LW_DATA`, `LW_GUIDE`, `MS`, `FS_ANNOT`, `PX_PER_UNIT` names match `fig_style.jl` exactly; exp66's `naive`/`speedup` are defined in Task 4 Step 2 before use in Step 3.

# Screened Density z-Upsampling Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a new figure comparing the TKM3D slab-max decay before and after global z-only NUFFT upsampling of the original density.

**Architecture:** Convert the centered tensor-product `VolumeSource` back into structured arrays, apply a type-1/type-2 NUFFT reconstruction along `z` for each `(x, y)` column, rebuild refined sources for upsampling factors `2` and `3`, and reuse the existing TKM3D slab-spectrum helper on the screened refined sources.

**Tech Stack:** Julia, BoundaryIntegral.jl, FINUFFT.jl, CairoMakie.jl, Test

---

### Task 1: Add a failing regression for z-only global upsampling

**Files:**
- Modify: `test/test_screened_density_analysis.jl`

**Step 1: Write the failing test**

Add a small tensor-grid source with nontrivial `z` variation and assert that z-only global upsampling:
- preserves the original samples at every `factor`-th refined `z` index
- keeps the `x` and `y` axes unchanged
- multiplies the number of `z` planes by the requested factor

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/test_screened_density_analysis.jl`

Expected: FAIL because the helper does not exist yet.

### Task 2: Implement z-only NUFFT reconstruction helpers

**Files:**
- Modify: `src/ScreenedDensityAnalysis.jl`

**Step 1: Write minimal implementation**

Add helpers to:
- extract tensor-product axes and arrays from a `VolumeSource`
- upsample each `(x, y)` density column in `z` with a type-1/type-2 NUFFT round-trip
- rebuild a refined `VolumeSource` with refined `z` weights

**Step 2: Run the focused test**

Run: `julia --project=. test/test_screened_density_analysis.jl`

Expected: PASS

### Task 3: Add the new z-upsampled decay figure

**Files:**
- Modify: `scripts/screened_density_tkm.jl`

**Step 1: Wire the new figure**

Create a new figure that, for each screening case, plots the normalized slab-max decay for:
- original grid
- `z×2`
- `z×3`

Save it to `figs/screened_density_kz_decay_z_upsampled.png`.

**Step 2: Run broad verification**

Run:
- `julia --project=. test/runtests.jl`
- `julia --project=. scripts/screened_density_tkm.jl`

Expected: PASS, and the new figure is regenerated.

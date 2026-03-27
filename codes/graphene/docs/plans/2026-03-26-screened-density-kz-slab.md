# Screened Density kz-Slab Spectrum Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the old `\hat{\rho}(0,0,k_z)` decay figure with a TKM3D-style 3D NUFFT kz-slab maximum spectrum.

**Architecture:** Add a density-analysis helper that constructs the same spectral grid/cutoff style used by the TKM3D continuous solver path, then reduce each `k_z` slab by `max(abs(...))` and feed that observable into the existing plotting script.

**Tech Stack:** Julia, BoundaryIntegral.jl, TKM3D.jl internals, FINUFFT.jl, CairoMakie.jl, Test

---

### Task 1: Add a failing regression for the new slab observable

**Files:**
- Modify: `test/test_screened_density_analysis.jl`

**Step 1: Write the failing test**

Use a synthetic 2x2x2 source cloud with a single active source charge and assert that the TKM3D kz-slab maximum spectrum is constant across retained nonnegative `k_z` values.

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/test_screened_density_analysis.jl`

Expected: FAIL because the new helper does not exist yet.

### Task 2: Implement the TKM3D-style NUFFT helper

**Files:**
- Modify: `src/ScreenedDensityAnalysis.jl`

**Step 1: Write minimal implementation**

Add:
- TKM3D-style geometry construction for the mode axes
- 3D type-1 NUFFT of `weights .* rho`
- kz-slab reduction by `max(abs(coeff))` over `(k_x, k_y)` inside the spherical `kmax` cutoff

**Step 2: Run the focused test**

Run: `julia --project=. test/test_screened_density_analysis.jl`

Expected: PASS

### Task 3: Switch the main decay figure over to the new observable

**Files:**
- Modify: `scripts/screened_density_tkm.jl`

**Step 1: Update the decay plot**

Replace the old `normalized_kz_decay` call with the new TKM3D slab-spectrum helper and update the axis label/title accordingly.

**Step 2: Run broad verification**

Run:
- `julia --project=. test/runtests.jl`
- `julia --project=. scripts/screened_density_tkm.jl`

Expected: PASS, and `figs/screened_density_kz_decay.png` is regenerated.

# Screened Density Analysis Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Generate the requested `x = 0` screened-density slice plots and `\hat{\rho}(0, 0, k_z)` decay comparison for sharp screening and multiple softmix bandwidths.

**Architecture:** Add a small reusable Julia module for screened-density analysis, keep plotting in the script entrypoint, and verify the math with focused unit tests on synthetic data.

**Tech Stack:** Julia, BoundaryIntegral.jl, CairoMakie, FFT-free quadrature for `k_z` spectra, Julia `Test`

---

### Task 1: Add failing tests for the analysis helpers

**Files:**
- Create: `test/runtests.jl`
- Create: `test/test_screened_density_analysis.jl`

**Step 1: Write the failing tests**

- Add one test for screened slice evaluation on a synthetic linear datagrid.
- Add one test for `\hat{\rho}(0, 0, k_z)` against direct quadrature on a small synthetic volume source.

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/runtests.jl`

Expected: FAIL because `src/ScreenedDensityAnalysis.jl` does not exist yet.

### Task 2: Implement the analysis module

**Files:**
- Create: `src/ScreenedDensityAnalysis.jl`

**Step 1: Write minimal implementation**

- Add helpers to compute box bounds, screened density vectors, screened density at arbitrary Cartesian points, `x = 0` slice sampling, and `\hat{\rho}(0, 0, k_z)` from plane-summed charge.

**Step 2: Run tests to verify they pass**

Run: `julia --project=. test/runtests.jl`

Expected: PASS

### Task 3: Wire the plotting script

**Files:**
- Modify: `scripts/screened_density_tkm.jl`

**Step 1: Update the script**

- Load the new module.
- Build sharp and softmix cases.
- Save the slice and spectral-decay figures.

**Step 2: Run the script**

Run: `julia --project=. scripts/screened_density_tkm.jl`

Expected: two figures saved under `figs/`.

### Task 4: Verify outputs

**Files:**
- Output: `figs/screened_density_xslice.png`
- Output: `figs/screened_density_kz_decay.png`

**Step 1: Re-run tests**

Run: `julia --project=. test/runtests.jl`

Expected: PASS

**Step 2: Confirm figures exist**

Run: `find figs -maxdepth 1 -type f | sort`

Expected: both new figure files present.

### Task 5: Switch the x-slice heatmap to log10 values

**Files:**
- Modify: `src/ScreenedDensityAnalysis.jl`
- Modify: `test/test_screened_density_analysis.jl`
- Modify: `scripts/screened_density_tkm.jl`

**Step 1: Write the failing test**

- Add a test that verifies positive values are transformed with `log10`, nonpositive values are clamped to a specified floor before transforming, and `NaN` values stay `NaN`.

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/runtests.jl`

Expected: FAIL because the log-transform helper does not exist yet.

**Step 3: Write minimal implementation**

- Add a reusable helper for the log10 transform.
- Update the slice figure to use the transformed values and a matching colorbar label.

**Step 4: Run tests and regenerate the figure**

Run: `julia --project=. test/runtests.jl`
Run: `julia --project=. scripts/screened_density_tkm.jl`

Expected: PASS, and `figs/screened_density_xslice.png` is updated.

### Task 6: Make centering exact and consistent between slice and spectrum

**Files:**
- Modify: `src/ScreenedDensityAnalysis.jl`
- Modify: `test/test_screened_density_analysis.jl`
- Modify: `scripts/screened_density_tkm.jl`

**Step 1: Write the failing tests**

- Add one test that the computed centering shift moves the density centroid to `(0, 0, 0)`.
- Add one test that the Cartesian slice sampler respects a shifted datagrid origin.

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/runtests.jl`

Expected: FAIL because the centering/shift helpers do not exist yet.

**Step 3: Write minimal implementation**

- Add helpers for density centroid, exact centering shift, and datagrid shifting.
- Return the shifted datagrid from the default loader and remove the hard-coded approximate shift from the script.

**Step 4: Run tests and regenerate the figures**

Run: `julia --project=. test/runtests.jl`
Run: `julia --project=. scripts/screened_density_tkm.jl`

Expected: PASS, and both figures use the same centered density.

### Task 7: Add a separate refined 1D Q(z) decay figure

**Files:**
- Modify: `src/ScreenedDensityAnalysis.jl`
- Modify: `test/test_screened_density_analysis.jl`
- Modify: `scripts/screened_density_tkm.jl`

**Step 1: Write the failing tests**

- Add one test that periodic spectral upsampling reproduces a band-limited cosine exactly on a refined grid.
- Add one test that the refined `Q(z)` spectrum reduces to the existing coarse `k_z` quadrature when `upsample_factor = 1` and `eps_in = eps_out`.

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/runtests.jl`

Expected: FAIL because the 1D refinement helpers do not exist yet.

**Step 3: Write minimal implementation**

- Add helpers for plane-charge extraction, spectral upsampling in `z`, and refined `Q(z)`-based `\hat{\rho}(0,0,k_z)`.
- Add a new figure written separately from the current coarse plot.

**Step 4: Run tests and regenerate figures**

Run: `julia --project=. test/runtests.jl`
Run: `julia --project=. scripts/screened_density_tkm.jl`

Expected: PASS, and a third figure for the refined 1D method is saved.

### Task 8: Switch the refined Q(z) path to type-1/type-2 NUFFT upsampling

**Files:**
- Modify: `src/ScreenedDensityAnalysis.jl`
- Modify: `test/test_screened_density_analysis.jl`

**Step 1: Write the failing test**

- Add a round-trip test that type-1 followed by type-2 NUFFT reproduces a periodic cosine on a refined grid.

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/runtests.jl`

Expected: FAIL because the NUFFT upsampling helper does not exist yet.

**Step 3: Write minimal implementation**

- Add a `Q(z)` periodic NUFFT upsampling helper using `FINUFFT.nufft1d1` and `FINUFFT.nufft1d2`.
- Make the refined `Q(z)` figure use the NUFFT helper.

**Step 4: Run tests and regenerate figures**

Run: `julia --project=. test/runtests.jl`
Run: `julia --project=. scripts/screened_density_tkm.jl`

Expected: PASS, and the refined `Q(z)` figure is regenerated using NUFFT.

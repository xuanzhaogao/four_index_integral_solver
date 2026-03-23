# Screened Gaussian Convolution Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend `scripts/screened_gauss.jl` to compare the regularized same-line `1 / r` convolution of the step-screened and `erf`-screened Gaussians using `HCubature`.

**Architecture:** Add a small `HCubature`-based integration helper that evaluates the exclusion-radius regularized line potential over piecewise subintervals, then build a study helper that applies it to the step screen and each `erf` bandwidth on a shared target grid. Plot the resulting potentials and their difference to the step-screened baseline in a separate figure while preserving the existing decay figure.

**Tech Stack:** Julia, HCubature, CairoMakie, FFTW, FINUFFT, SpecialFunctions, Test

---

### Task 1: Add focused tests for the regularized convolution helper

**Files:**
- Modify: `/Users/xgao/Works/four_indices_integral_solver/codes/volume_source/test/screened_gauss.jl`
- Test: `/Users/xgao/Works/four_indices_integral_solver/codes/volume_source/test/screened_gauss.jl`

**Step 1: Write the failing test**

```julia
@testset "regularized line convolution" begin
    xmin = -2.0
    xmax = 2.0
    x0 = 0.25
    delta = 0.1

    result = regularized_line_convolution(
        x -> 1.0,
        x0;
        xmin,
        xmax,
        delta,
        rtol = 1.0e-10,
        atol = 1.0e-12,
    )

    expected = log(((x0 - xmin) * (xmax - x0)) / delta^2)

    @test result.value ≈ expected atol = 1.0e-8 rtol = 1.0e-8
    @test result.error >= 0.0
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/screened_gauss.jl`
Expected: FAIL because `regularized_line_convolution` does not exist yet

**Step 3: Write minimal implementation**

Add a reusable `regularized_line_convolution` helper in
`scripts/screened_gauss.jl` that splits the finite interval around the excluded
singular neighborhood and integrates each remaining piece with `hcubature`.

**Step 4: Run test to verify it passes**

Run: `julia --project=. test/screened_gauss.jl`
Expected: PASS

**Step 5: Commit**

```bash
git add test/screened_gauss.jl scripts/screened_gauss.jl
git commit -m "feat: add regularized line convolution helper"
```

### Task 2: Build the screened convolution study and comparison figure

**Files:**
- Modify: `/Users/xgao/Works/four_indices_integral_solver/codes/volume_source/scripts/screened_gauss.jl`
- Modify: `/Users/xgao/Works/four_indices_integral_solver/codes/volume_source/test/screened_gauss.jl`

**Step 1: Write the failing test**

Extend the same test file with a small-shape study check.

```julia
study = run_screened_line_convolution_study(;
    sigma = 1.0,
    bandwidths = [0.2],
    source_half_width_sigmas = 4.0,
    target_half_width_sigmas = 1.0,
    target_npoints = 9,
    exclusion_sigmas = 0.05,
    rtol = 1.0e-6,
    atol = 1.0e-8,
)

@test study.step_curve.label == "step"
@test length(study.targets) == 9
@test study.step_curve.targets == study.targets
@test length(study.curves) == 1
@test study.curves[1].targets == study.targets
@test all(isfinite, study.step_curve.potential)
@test all(isfinite, study.curves[1].potential)
@test all(isfinite, study.curves[1].difference_from_step)
```

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/screened_gauss.jl`
Expected: FAIL because `run_screened_line_convolution_study` does not exist yet

**Step 3: Write minimal implementation**

Update `scripts/screened_gauss.jl` to:

- import `HCubature`
- add the line-convolution study helper
- build the step and `erf` convolution curves on a shared target grid
- create a separate convolution comparison figure
- restore the `PROGRAM_FILE` guard so tests can include the script without
  running the plotting path

**Step 4: Run test to verify it passes**

Run: `julia --project=. test/screened_gauss.jl`
Expected: PASS

**Step 5: Commit**

```bash
git add test/screened_gauss.jl scripts/screened_gauss.jl
git commit -m "feat: compare screened Gaussian line convolutions"
```

### Task 3: Verify the script end to end

**Files:**
- Modify: `/Users/xgao/Works/four_indices_integral_solver/codes/volume_source/scripts/screened_gauss.jl`
- Test: `/Users/xgao/Works/four_indices_integral_solver/codes/volume_source/test/screened_gauss.jl`

**Step 1: Run the focused test file**

Run: `julia --project=. test/screened_gauss.jl`
Expected: PASS

**Step 2: Run the script**

Run: `julia --project=. scripts/screened_gauss.jl`
Expected: The script completes without error and writes the updated figure files

**Step 3: Confirm include-safe loading**

Run: `julia --project=. -e 'include("scripts/screened_gauss.jl")'`
Expected: Exit 0 with no plotting side effects

**Step 4: Inspect the figure paths**

Run: `ls -l figs/screened_gauss_erf_decay.svg figs/screened_gauss_line_convolution.svg`
Expected: Both configured output paths exist

**Step 5: Commit**

```bash
git add docs/plans/2026-03-19-screened-gauss-convolution-design.md docs/plans/2026-03-19-screened-gauss-convolution.md test/screened_gauss.jl scripts/screened_gauss.jl
git commit -m "docs: plan screened Gaussian convolution comparison"
```

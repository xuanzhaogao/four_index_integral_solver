# Erf-Screened Gaussian Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend `scripts/screened_gauss.jl` to plot the decay of `|fhat(k)|` for `erf`-screened Gaussians across multiple bandwidths with both the hard-step screened case and the unscreened Gaussian in the same figure.

**Architecture:** Keep the existing script as the home for one-dimensional screening studies, factor the numerical study code into reusable helpers, and add an `erf`-screening path that computes FFT-based decay curves on a common grid. The output figure becomes an overlay comparison of positive-frequency magnitude tails on a log-scale axis, with the hard-step screened Gaussian and the exact unscreened Gaussian added as shared reference curves.

**Tech Stack:** Julia, FFTW, CairoMakie, SpecialFunctions, Test

---

### Task 1: Add focused tests for the erf-screening study helpers

**Files:**
- Create: `/Users/xgao/Works/four_indices_integral_solver/codes/volume_source/test/screened_gauss.jl`
- Test: `/Users/xgao/Works/four_indices_integral_solver/codes/volume_source/test/screened_gauss.jl`

**Step 1: Write the failing test**

```julia
using Test

include(joinpath(@__DIR__, "..", "scripts", "screened_gauss.jl"))

@testset "erf-screened Gaussian decay helpers" begin
    sigma = 1.0
    narrow = 0.2
    wide = 1.0

    @test erf_screen_function(-5.0, narrow) < 1.0e-6
    @test erf_screen_function(5.0, narrow) > 1.0 - 1.0e-6
    @test erf_screen_function(0.0, narrow) ≈ 0.5
    @test erf_screen_function(0.5, narrow) > erf_screen_function(0.5, wide)

    study = run_erf_screen_decay_study(; sigma, bandwidths = [narrow, wide], npoints = 64)

    @test length(study.bandwidths) == 2
    @test length(study.curves) == 2
    @test all(length(curve.k_positive) == length(curve.magnitude_positive) for curve in study.curves)
    @test all(all(mag .>= 0.0) for mag in getfield.(study.curves, :magnitude_positive))
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/screened_gauss.jl`
Expected: FAIL because `erf_screen_function` and `run_erf_screen_decay_study` do not exist yet

**Step 3: Write minimal implementation**

Add the `erf` screening helper and a reusable decay-study helper in `scripts/screened_gauss.jl`, returning a structured result with one positive-frequency magnitude curve per bandwidth.

**Step 4: Run test to verify it passes**

Run: `julia --project=. test/screened_gauss.jl`
Expected: PASS

**Step 5: Commit**

```bash
git add test/screened_gauss.jl scripts/screened_gauss.jl
git commit -m "feat: add erf-screened Gaussian decay study"
```

### Task 2: Replace the current figure build with the bandwidth overlay decay plot

**Files:**
- Modify: `/Users/xgao/Works/four_indices_integral_solver/codes/volume_source/scripts/screened_gauss.jl`
- Test: `/Users/xgao/Works/four_indices_integral_solver/codes/volume_source/test/screened_gauss.jl`

**Step 1: Write the failing test**

Extend the same test file with assertions that the positive-frequency axis is sorted, starts at `0`, is shared across bandwidth curves when the common grid parameters are reused, and is also shared by the hard-step and unscreened reference curves.

```julia
@test first(study.curves[1].k_positive) ≈ 0.0
@test issorted(study.curves[1].k_positive)
@test study.curves[1].k_positive == study.curves[2].k_positive
@test study.step_curve.label == "step"
@test study.step_curve.k_positive == study.curves[1].k_positive
@test study.unscreened_curve.label == "unscreened"
@test study.unscreened_curve.k_positive == study.curves[1].k_positive
```

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/screened_gauss.jl`
Expected: FAIL until the helper returns the shared positive-frequency view expected by the plot path

**Step 3: Write minimal implementation**

Update the script’s plotting block to:

- define a user-editable bandwidth list
- call `run_erf_screen_decay_study`
- build a log-scale `|fhat(k)|` overlay for the positive-frequency tail
- label each curve by bandwidth
- add the hard-step screened case as a reference curve
- add the unscreened Gaussian as an exact reference curve
- restore the `PROGRAM_FILE` guard so tests can include the script without plot side effects
- save or return the resulting figure

**Step 4: Run test to verify it passes**

Run: `julia --project=. test/screened_gauss.jl`
Expected: PASS

**Step 5: Commit**

```bash
git add test/screened_gauss.jl scripts/screened_gauss.jl
git commit -m "feat: plot erf-screened Gaussian decay overlays"
```

### Task 3: Verify the script end to end

**Files:**
- Modify: `/Users/xgao/Works/four_indices_integral_solver/codes/volume_source/scripts/screened_gauss.jl`
- Test: `/Users/xgao/Works/four_indices_integral_solver/codes/volume_source/test/screened_gauss.jl`

**Step 1: Run the focused tests**

Run: `julia --project=. test/screened_gauss.jl`
Expected: PASS

**Step 2: Run the script**

Run: `julia --project=. scripts/screened_gauss.jl`
Expected: The script completes without error and produces the decay comparison figure

**Step 3: Inspect the figure path**

Run: `ls -l fig figs`
Expected: The configured output path exists or the script returns the figure without path errors

**Step 4: Commit**

```bash
git add docs/plans/2026-03-19-erf-screened-gauss-design.md docs/plans/2026-03-19-erf-screened-gauss.md test/screened_gauss.jl scripts/screened_gauss.jl
git commit -m "docs: plan erf-screened Gaussian decay visualization"
```

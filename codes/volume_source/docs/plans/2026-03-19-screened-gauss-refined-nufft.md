# Screened Gaussian Refined-NUFFT Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the current uniform-grid decay study in `scripts/screened_gauss.jl` with a locally refined, weighted type-1 NUFFT that resolves narrow `erf` screening bandwidths near `x = 0`.

**Architecture:** Build a small nonuniform midpoint quadrature helper and a weighted type-1 NUFFT helper, then route the decay study through those helpers while keeping the convolution study unchanged. The plotted `k` axis becomes `k * sigma` because there is no single uniform `dx` after local refinement.

**Tech Stack:** Julia, FINUFFT, CairoMakie, FFTW, HCubature, SpecialFunctions, Test

---

### Task 1: Add tests for the refined quadrature and weighted NUFFT helpers

**Files:**
- Modify: `/Users/xgao/Works/four_indices_integral_solver/codes/volume_source/test/screened_gauss.jl`
- Test: `/Users/xgao/Works/four_indices_integral_solver/codes/volume_source/test/screened_gauss.jl`

**Step 1: Write the failing test**

```julia
@testset "refined quadrature and weighted NUFFT" begin
    quad = screen_locally_refined_quadrature(;
        sigma = 1.0,
        half_width_sigmas = 6.0,
        npoints = 32,
        refinement_half_width_sigmas = 0.5,
        refinement_factor = 8,
    )

    @test length(quad.x) == length(quad.weights)
    @test all(quad.weights .> 0.0)
    @test sum(quad.weights) ≈ 12.0 atol = 1.0e-12
    @test minimum(quad.weights) < maximum(quad.weights)

    values = gaussian_screen_base.(quad.x, Ref(1.0))
    k, fhat = screen_continuous_nufft_type1_weighted(
        quad.x,
        values,
        quad.weights,
        quad.interval_length;
        mode_extent = 32,
    )

    positive = findall(k .>= 0.0)
    reference = gaussian_screen_base_ft.(k[positive][1:5], Ref(1.0))

    @test maximum(abs.(real.(fhat[positive][1:5]) .- reference)) < 1.0e-3
    @test maximum(abs.(imag.(fhat[positive][1:5]))) < 1.0e-3
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/screened_gauss.jl`
Expected: FAIL because the refined quadrature and weighted NUFFT helpers do not exist yet

**Step 3: Write minimal implementation**

Add the refined quadrature helper and the weighted NUFFT helper in
`scripts/screened_gauss.jl`.

**Step 4: Run test to verify it passes**

Run: `julia --project=. test/screened_gauss.jl`
Expected: PASS

**Step 5: Commit**

```bash
git add test/screened_gauss.jl scripts/screened_gauss.jl
git commit -m "feat: add refined weighted nufft helpers"
```

### Task 2: Route the decay study through the refined weighted NUFFT

**Files:**
- Modify: `/Users/xgao/Works/four_indices_integral_solver/codes/volume_source/scripts/screened_gauss.jl`
- Modify: `/Users/xgao/Works/four_indices_integral_solver/codes/volume_source/test/screened_gauss.jl`

**Step 1: Write the failing test**

Extend the existing decay-study test to require refined-grid metadata and a
shared positive-frequency axis.

```julia
study = run_erf_screen_decay_study(;
    sigma = 1.0,
    bandwidths = [0.05, 0.2],
    half_width_sigmas = 4.0,
    npoints = 32,
    refinement_half_width_sigmas = 0.5,
    refinement_factor = 8,
    mode_extent = 64,
)

@test study.min_weight < study.coarse_weight
@test study.step_curve.k_positive == study.curves[1].k_positive
@test study.unscreened_curve.k_positive == study.curves[1].k_positive
```

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/screened_gauss.jl`
Expected: FAIL until the decay study exposes the refined-grid path

**Step 3: Write minimal implementation**

Update `run_erf_screen_decay_study` and `build_erf_screen_decay_figure` to:

- use the refined quadrature and weighted NUFFT
- evaluate the step and `erf` curves on the nonuniform weighted rule
- keep the unscreened Gaussian as an exact reference on the same `k` grid
- change the horizontal axis label to `k * sigma`
- restore the `PROGRAM_FILE` guard so tests can include the script safely

**Step 4: Run test to verify it passes**

Run: `julia --project=. test/screened_gauss.jl`
Expected: PASS

**Step 5: Commit**

```bash
git add test/screened_gauss.jl scripts/screened_gauss.jl
git commit -m "feat: use refined weighted nufft for screened gaussian decay"
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
Expected: The script completes without error and writes both figure files

**Step 3: Confirm include-safe loading**

Run: `julia --project=. -e 'include("scripts/screened_gauss.jl")'`
Expected: Exit 0 with no plotting side effects

**Step 4: Inspect the figure paths**

Run: `ls -l figs/screened_gauss_erf_decay.svg figs/screened_gauss_line_convolution.svg`
Expected: Both configured output paths exist

**Step 5: Commit**

```bash
git add docs/plans/2026-03-19-screened-gauss-refined-nufft-design.md docs/plans/2026-03-19-screened-gauss-refined-nufft.md test/screened_gauss.jl scripts/screened_gauss.jl
git commit -m "docs: plan refined nufft screened gaussian decay study"
```

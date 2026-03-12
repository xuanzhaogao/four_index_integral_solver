# Max Frequency Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a script that computes the smallest radial Fourier cutoff for the squared graphene orbital density such that the relative pointwise tail stays below a chosen tolerance.

**Architecture:** Keep the analysis in a standalone script with reusable helpers. Test the cutoff-selection helper independently, then wire it to the XSF loader and one type-1 NUFFT over the full Nyquist box.

**Tech Stack:** Julia, BoundaryIntegral, TKM3D/FINUFFT, standalone script execution

---

### Task 1: Add a failing cutoff-helper test

**Files:**
- Create: `test/max_frequency.jl`
- Modify: `scripts/max_frequency.jl`

**Step 1: Write the failing test**

Create a focused test for a helper like `_smallest_cutoff_from_tail(radii, mags, tol)`
using a synthetic spectrum where the correct cutoff is known.

**Step 2: Run test to verify it fails**

Run: `julia --project=. --color=no test/max_frequency.jl`

Expected: FAIL because the helper is not implemented yet.

**Step 3: Write minimal implementation**

Add the helper and enough script structure for the test to load it.

**Step 4: Run test to verify it passes**

Run: `julia --project=. --color=no test/max_frequency.jl`

Expected: PASS.

### Task 2: Implement orbital grid to type-1 NUFFT analysis

**Files:**
- Modify: `scripts/max_frequency.jl`
- Test: `test/max_frequency.jl`

**Step 1: Add helpers**

Implement helpers to:

- load and square the orbital density
- convert the datagrid into `3 x N` source coordinates and weighted masses
- construct centered mode axes from the grid spacings
- build radial `|k|` values for the NUFFT coefficient box
- evaluate and summarize the tail checkpoints

**Step 2: Wire the script entrypoint**

Implement `main()` to run the analysis on
`density_data/graphene_00001_5x5x1_shifted.xsf` with a top-level `tol`.

**Step 3: Re-run the focused test**

Run: `julia --project=. --color=no test/max_frequency.jl`

Expected: PASS.

### Task 3: Verify on the real orbital

**Files:**
- Modify: `scripts/max_frequency.jl` if needed

**Step 1: Run the script**

Run: `julia --project=. --color=no scripts/max_frequency.jl`

Expected: prints grid metadata, Nyquist limits, spectral maximum, tail checkpoints, and the estimated `k_cut`.

**Step 2: If runtime or memory issues appear, trim only the reporting layer**

Do not change the raw-grid spectral analysis unless the result is technically invalid.

### Task 4: Final verification

**Files:**
- Verify: `test/max_frequency.jl`
- Verify: `scripts/max_frequency.jl`

**Step 1: Run the focused test fresh**

Run: `julia --project=. --color=no test/max_frequency.jl`

Expected: PASS.

**Step 2: Run the real script fresh**

Run: `julia --project=. --color=no scripts/max_frequency.jl`

Expected: completes and reports a concrete `k_cut` for the graphene orbital.

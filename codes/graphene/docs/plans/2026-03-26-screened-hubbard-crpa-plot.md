# Screened Hubbard cRPA Plot Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Plot the computed screened Hubbard interactions against softmix bandwidth and compare them to graphene cRPA reference values from the local PRL paper.

**Architecture:** Add a small helper module that owns the cRPA constants and CSV-to-curve transformation, then build a `CairoMakie` script on top of it. Cover the data transformation with focused tests so the plot has a verified input path.

**Tech Stack:** Julia, CSV.jl, DataFrames.jl, CairoMakie.jl, Test

---

### Task 1: Add failing tests for the comparison data path

**Files:**
- Create: `test/test_screened_hubbard_comparison.jl`
- Modify: `test/runtests.jl`

**Step 1: Write the failing tests**

Add tests that:
- assert the graphene cRPA constants are `9.3`, `5.5`, `4.1`, `3.6` for `U_00` through `U_03`
- build a small unsorted `DataFrame` and assert the helper returns pair curves sorted by bandwidth in the canonical pair order

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/test_screened_hubbard_comparison.jl`

Expected: FAIL because the helper module does not exist yet.

### Task 2: Implement the comparison helper

**Files:**
- Create: `src/ScreenedHubbardComparison.jl`

**Step 1: Write minimal implementation**

Implement:
- `DEFAULT_SCREENED_HUBBARD_CSV`
- `PAIR_ORDER`
- `GRAPHENE_CRPA_EV`
- `load_softmix_hubbard_results`
- `pair_curve_data`

**Step 2: Run the new tests**

Run: `julia --project=. test/test_screened_hubbard_comparison.jl`

Expected: PASS

### Task 3: Add the plot script and verify output

**Files:**
- Create: `scripts/plot_screened_hubbard_vs_crpa.jl`

**Step 1: Generate the combined plot**

Use `CairoMakie` to draw:
- four computed curves versus `b`
- four matching horizontal cRPA reference lines

Save to `figs/screened_hubbard_vs_crpa.png`.

**Step 2: Run broad verification**

Run:
- `julia --project=. test/runtests.jl`
- `julia --project=. scripts/plot_screened_hubbard_vs_crpa.jl`

Expected: PASS, and the figure file exists.

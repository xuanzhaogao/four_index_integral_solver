# Screened Orbital Solve Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a full screened dielectric-box orbital workflow that evaluates `u_int + u_scatter` for the graphene `U_00/U_01/U_02/U_03` channels over a bandwidth sweep.

**Architecture:** Add a new reusable module for centered orbital loading, pair construction, solver execution, and target integration. Keep the top-level sweep and CSV writing in a new script, and protect the core helpers with focused unit tests before using them in the end-to-end run.

**Tech Stack:** Julia, BoundaryIntegral.jl, BoundaryIntegral.TKM3D, Krylov.jl, CSV.jl, DataFrames.jl, Test stdlib

---

### Task 1: Add failing tests for orbital-solve helpers

**Files:**
- Create: `/mnt/home/xgao1/work/four_index_integral_solver/codes/graphene/test/test_screened_orbital_solve.jl`
- Modify: `/mnt/home/xgao1/work/four_index_integral_solver/codes/graphene/test/runtests.jl`
- Test: `/mnt/home/xgao1/work/four_index_integral_solver/codes/graphene/test/test_screened_orbital_solve.jl`

**Step 1: Write the failing test**

Add tests for:
- direct volume-potential evaluation against a small manual `1/(4πr)` sum,
- orbital-potential integration `sum(potential .* density .* weights)`,
- pair-source construction with the `a1` shift.

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/test_screened_orbital_solve.jl`
Expected: FAIL because the new module/functions do not exist yet.

**Step 3: Write minimal implementation**

Create only the helper signatures needed for the tests in the new source module.

**Step 4: Run test to verify it passes**

Run: `julia --project=. test/test_screened_orbital_solve.jl`
Expected: PASS

### Task 2: Implement reusable screened orbital solve helpers

**Files:**
- Create: `/mnt/home/xgao1/work/four_index_integral_solver/codes/graphene/src/ScreenedOrbitalSolve.jl`
- Modify: `/mnt/home/xgao1/work/four_index_integral_solver/codes/graphene/test/test_screened_orbital_solve.jl`

**Step 1: Write the failing test**

Add tests that exercise:
- source centering/pair construction on synthetic data,
- zero-scatter total field reduction when `sigma = 0`,
- mode sweep labeling/output row assembly.

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/test_screened_orbital_solve.jl`
Expected: FAIL with missing helper behavior.

**Step 3: Write minimal implementation**

Implement:
- centered graphene orbital loaders,
- `graphene_pair_sources`,
- `evaluate_volume_potential`,
- `integrate_target_potential`,
- `solve_screened_mode`,
- row assembly helpers for CSV output.

**Step 4: Run test to verify it passes**

Run: `julia --project=. test/test_screened_orbital_solve.jl`
Expected: PASS

### Task 3: Add the full-solve script

**Files:**
- Create: `/mnt/home/xgao1/work/four_index_integral_solver/codes/graphene/scripts/screened_hubbard_graphene.jl`
- Modify: `/mnt/home/xgao1/work/four_index_integral_solver/codes/graphene/src/ScreenedOrbitalSolve.jl`

**Step 1: Write the failing test**

Add a small smoke-level test for script-facing helpers, such as the default mode list and constant parameter wiring.

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/test_screened_orbital_solve.jl`
Expected: FAIL until the script-facing helpers are implemented.

**Step 3: Write minimal implementation**

Implement the script to:
- load centered orbitals,
- build the four target channels,
- loop over `Sharp` and SoftMix bandwidths,
- solve each mode,
- write `data/screened_hubbard_graphene.csv`,
- print the table.

**Step 4: Run test to verify it passes**

Run: `julia --project=. test/test_screened_orbital_solve.jl`
Expected: PASS

### Task 4: Verify the full workflow

**Files:**
- Modify: `/mnt/home/xgao1/work/four_index_integral_solver/codes/graphene/test/runtests.jl`
- Output: `/mnt/home/xgao1/work/four_index_integral_solver/codes/graphene/data/screened_hubbard_graphene.csv`

**Step 1: Run the test suite**

Run: `julia --project=. test/runtests.jl`
Expected: PASS

**Step 2: Run the full solver script**

Run: `julia --project=. scripts/screened_hubbard_graphene.jl`
Expected: PASS, and writes `data/screened_hubbard_graphene.csv`.

**Step 3: Spot-check output**

Confirm the CSV contains rows for:
- `Sharp`
- each SoftMix bandwidth
- `U_00`, `U_01`, `U_02`, `U_03`

**Step 4: Commit**

```bash
git add docs/plans/2026-03-25-screened-orbital-solve-design.md docs/plans/2026-03-25-screened-orbital-solve.md src/ScreenedOrbitalSolve.jl scripts/screened_hubbard_graphene.jl test/test_screened_orbital_solve.jl test/runtests.jl data/screened_hubbard_graphene.csv
git commit -m "feat: add screened graphene orbital solve workflow"
```

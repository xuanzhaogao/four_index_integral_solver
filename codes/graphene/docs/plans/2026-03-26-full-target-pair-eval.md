# Full-Target Pair Evaluation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Remove target batching from screened pair evaluation so each pair channel evaluates the full target orbital in one pass.

**Architecture:** Refactor `evaluate_screened_pair_interaction` to operate on the full target matrix and integrate once, then delete the now-unused `batch_size` plumbing from the solve and script layers. Cover the behavior with a regression that counts evaluator calls through injected test doubles.

**Tech Stack:** Julia, BoundaryIntegral.jl, Krylov, Test

---

### Task 1: Add the regression for one-shot target evaluation

**Files:**
- Modify: `test/test_screened_orbital_solve.jl`

**Step 1: Write the failing test**

Add a test that:
- constructs small synthetic `VolumeSource`s
- injects fake volume/scatter evaluators
- asserts each evaluator is called once on the full target set

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/test_screened_orbital_solve.jl`

Expected: FAIL because `evaluate_screened_pair_interaction` does not yet accept the injected evaluators and still uses the old chunked path.

### Task 2: Remove the target batching loop

**Files:**
- Modify: `src/ScreenedOrbitalSolve.jl`

**Step 1: Write minimal implementation**

Change `evaluate_screened_pair_interaction` to:
- evaluate `u_int` on `target_vs.positions`
- evaluate `u_scatter` on `target_vs.positions`
- integrate both vectors once against `target_vs.weights .* target_vs.density`

Remove the `batch_size` keyword from `evaluate_screened_pair_interaction` and `solve_screened_mode`.

**Step 2: Run test to verify it passes**

Run: `julia --project=. test/test_screened_orbital_solve.jl`

Expected: PASS

### Task 3: Remove stale script plumbing and verify broader behavior

**Files:**
- Modify: `scripts/screened_hubbard_graphene.jl`
- Test: `test/runtests.jl`

**Step 1: Delete the unused `BATCH_SIZE` constant and call-site forwarding**

Keep the top-level solve path aligned with the new full-target behavior.

**Step 2: Run broader verification**

Run: `julia --project=. test/runtests.jl`

Expected: PASS

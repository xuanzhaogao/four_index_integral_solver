# Paper Shell Labeling Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the slab-model `U_00` through `U_05` channels explicitly follow the paper-style shell ordering and document the dielectric-model mismatch.

**Architecture:** Keep the slab solver unchanged and only clarify the pair layout contract. Introduce one explicit paper-shell layout definition in the orbital helper module, verify it with a regression on shell distances, and update the report so the output is not mistaken for the paper’s full bilayer WFCE model.

**Tech Stack:** Julia, BoundaryIntegral.jl, Test, Markdown docs

---

### Task 1: Add a failing regression for the paper shell contract

**Files:**
- Modify: `test/test_screened_orbital_solve.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "paper shell layout labels are monotone in shell distance" begin
    xs = [0.0]
    ys = [0.0]
    zs = [6.7 / 4]
    weights = ones(1, 1, 1)
    density_1 = reshape([1.0], 1, 1, 1)
    density_2 = reshape([1.0], 1, 1, 1)

    vs1 = BI.VolumeSource((xs, ys, zs), weights, density_1)
    vs2 = BI.VolumeSource((xs, ys, zs), weights, density_2)
    pairs = pair_targets(vs1, vs2; layout = PAPER_SHELL_LAYOUT)

    distances = [
        norm(pairs.U_00.positions[:, 1] .- vs1.positions[:, 1]),
        norm(pairs.U_01.positions[:, 1] .- vs1.positions[:, 1]),
        norm(pairs.U_02.positions[:, 1] .- vs1.positions[:, 1]),
        norm(pairs.U_03.positions[:, 1] .- vs1.positions[:, 1]),
        norm(pairs.U_04.positions[:, 1] .- vs1.positions[:, 1]),
        norm(pairs.U_05.positions[:, 1] .- vs1.positions[:, 1]),
    ]

    @test distances == sort(distances)
end
```

- [ ] **Step 2: Run the focused test to verify it fails**

Run: `julia --project=. test/test_screened_orbital_solve.jl`
Expected: FAIL because `PAPER_SHELL_LAYOUT` does not exist yet.

- [ ] **Step 3: Commit**

```bash
git add test/test_screened_orbital_solve.jl
git commit -m "test: cover paper shell labeling contract"
```

### Task 2: Make the paper shell layout explicit in the solver helpers

**Files:**
- Modify: `src/ScreenedOrbitalSolve.jl`
- Test: `test/test_screened_orbital_solve.jl`

- [ ] **Step 1: Write the minimal implementation**

```julia
const PAPER_SHELL_LAYOUT = (
    (pair = :U_00, orbital = :vs1, shift = (0.0, 0.0, 0.0), description = "on-site"),
    (pair = :U_01, orbital = :vs2, shift = (0.0, 0.0, 0.0), description = "nearest neighbor"),
    (pair = :U_02, orbital = :vs1, shift = DEFAULT_A1, description = "next-nearest neighbor"),
    (pair = :U_03, orbital = :vs2, shift = DEFAULT_A1, description = "third shell"),
    (pair = :U_04, orbital = :vs2, shift = DEFAULT_A2, description = "fourth shell representative"),
    (pair = :U_05, orbital = :vs1, shift = (DEFAULT_A1[1] + DEFAULT_A2[1], DEFAULT_A1[2] + DEFAULT_A2[2], 0.0), description = "fifth shell representative"),
)
```

Use `PAPER_SHELL_LAYOUT` as the default layout for `pair_targets` and `pair_specs`.

- [ ] **Step 2: Run the focused test to verify it passes**

Run: `julia --project=. test/test_screened_orbital_solve.jl`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add src/ScreenedOrbitalSolve.jl test/test_screened_orbital_solve.jl
git commit -m "feat: make paper shell labeling explicit"
```

### Task 3: Document the approximation boundary

**Files:**
- Modify: `results/2026-03-26-screened-graphene-report.md`

- [ ] **Step 1: Update the report text**

```markdown
Add a note that:
- `U_00` through `U_05` follow the shell ordering used in Rosner et al.
- the current workflow remains a dielectric-slab approximation
- the paper’s true bilayer model instead uses AB-stacked bilayer graphene and a momentum-dependent effective dielectric function `ε_eff^2D(q)`
```

- [ ] **Step 2: Run the full test suite**

Run: `julia --project=. test/runtests.jl`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add results/2026-03-26-screened-graphene-report.md
git commit -m "docs: clarify paper shell labeling approximation"
```

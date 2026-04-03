# Bilayer Screened Graphene Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the screened graphene workflow to evaluate bilayer-as-slab interactions `U_00` through `U_05` at slab thickness `6.7 Å` with orbital centers fixed at `z = d/4`.

**Architecture:** Keep the dielectric-box solve path intact and modify only the orbital geometry helpers and the production script inputs. Represent the six interaction channels as an explicit table of pair labels and translation vectors so tests can verify the geometry directly and the script can emit the March 26 report style outputs for the new slab case.

**Tech Stack:** Julia, BoundaryIntegral.jl, Krylov, CSV.jl, DataFrames.jl, Test

---

### Task 1: Add geometry regression tests

**Files:**
- Modify: `test/test_screened_orbital_solve.jl`

- [ ] **Step 1: Write the failing tests**

```julia
@testset "pair targets support six graphene channels" begin
    xs = [0.0]
    ys = [0.0]
    zs = [1.675]
    weights = ones(1, 1, 1)
    density_1 = reshape([1.0], 1, 1, 1)
    density_2 = reshape([2.0], 1, 1, 1)

    vs1 = BI.VolumeSource((xs, ys, zs), weights, density_1)
    vs2 = BI.VolumeSource((xs, ys, zs), weights, density_2)
    pairs = pair_targets(vs1, vs2)

    @test keys(pairs) == (:U_00, :U_01, :U_02, :U_03, :U_04, :U_05)
    @test pairs.U_00.positions[:, 1] ≈ [0.0, 0.0, 1.675] atol = 1e-12
    @test pairs.U_01.positions[:, 1] ≈ [0.0, 0.0, 1.675] atol = 1e-12
end

@testset "centered graphene sources honor requested z center" begin
    sources = centered_graphene_sources(; tol = 1e-3, z_center = 1.675)
    z1 = sum(sources.vs1.positions[3, :] .* sources.vs1.weights .* sources.vs1.density) / sum(sources.vs1.weights .* sources.vs1.density)
    z2 = sum(sources.vs2.positions[3, :] .* sources.vs2.weights .* sources.vs2.density) / sum(sources.vs2.weights .* sources.vs2.density)

    @test z1 ≈ 1.675 atol = 1e-3
    @test z2 ≈ 1.675 atol = 1e-3
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `julia --project=. test/test_screened_orbital_solve.jl`
Expected: FAIL because `pair_targets` still returns four channels and `centered_graphene_sources` does not yet accept `z_center`.

- [ ] **Step 3: Commit**

```bash
git add test/test_screened_orbital_solve.jl
git commit -m "test: cover bilayer slab orbital geometry"
```

### Task 2: Implement configurable z centering and six-channel pairs

**Files:**
- Modify: `src/ScreenedOrbitalSolve.jl`
- Test: `test/test_screened_orbital_solve.jl`

- [ ] **Step 1: Write minimal implementation**

```julia
const DEFAULT_PAIR_TRANSLATIONS = (
    U_00 = ((0.0, 0.0, 0.0), :vs1),
    U_01 = ((0.0, 0.0, 0.0), :vs2),
    U_02 = ((2.465, 0.0, 0.0), :vs1),
    U_03 = ((2.465, 0.0, 0.0), :vs2),
    U_04 = ((...), :vs1),
    U_05 = ((...), :vs2),
)

function centered_graphene_sources(; orbital_1 = DEFAULT_ORBITAL_1, orbital_2 = DEFAULT_ORBITAL_2, tol = 1e-3, z_center = 0.0)
    ...
end

function pair_targets(vs1, vs2; pair_translations = DEFAULT_PAIR_TRANSLATIONS)
    ...
end

function pair_specs(vs1, vs2; pair_translations = DEFAULT_PAIR_TRANSLATIONS)
    ...
end
```

- [ ] **Step 2: Run tests to verify they pass**

Run: `julia --project=. test/test_screened_orbital_solve.jl`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add src/ScreenedOrbitalSolve.jl test/test_screened_orbital_solve.jl
git commit -m "feat: add bilayer slab orbital channels"
```

### Task 3: Update the production run entrypoint

**Files:**
- Modify: `scripts/screened_hubbard_graphene.jl`

- [ ] **Step 1: Update the run configuration**

```julia
const LZ = 6.7
const Z_CENTER = LZ / 4

function main(args = ARGS)
    sources = centered_graphene_sources(tol = SOURCE_TOL, z_center = Z_CENTER)
    specs = pair_specs(sources.vs1, sources.vs2)
    ...
end
```

- [ ] **Step 2: Ensure output metadata matches the new geometry**

```julia
push!(rows, (
    ...
    pair = pair_result.pair,
    shift_z = sources.shared_shift[3],
    orbital_z_center = Z_CENTER,
    slab_thickness = LZ,
))
```

- [ ] **Step 3: Run the focused tests**

Run: `julia --project=. test/test_screened_orbital_solve.jl`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add scripts/screened_hubbard_graphene.jl
git commit -m "feat: configure bilayer slab screened run"
```

### Task 4: Execute and verify the bilayer slab solve

**Files:**
- Modify: `data/screened_hubbard_graphene.csv`
- Modify: `results/2026-03-26-screened-graphene-report.md`

- [ ] **Step 1: Run the production solve**

Run: `julia --project=. scripts/screened_hubbard_graphene.jl sharp`
Expected: script completes and writes a six-row table containing `U_00` through `U_05`.

- [ ] **Step 2: Record the verified results in the report**

```markdown
Update the report's problem setup and workflow sections to state:
- `Lz = 6.7`
- orbital centers at `z = 1.675`
- output now includes `U_00` through `U_05`
Add the computed sharp-screening values to the report.
```

- [ ] **Step 3: Run the full test suite**

Run: `julia --project=. test/runtests.jl`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add data/screened_hubbard_graphene.csv results/2026-03-26-screened-graphene-report.md
git commit -m "feat: record bilayer slab screened graphene interactions"
```

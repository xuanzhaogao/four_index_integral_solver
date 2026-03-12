# TKM3D Hubbard Convergence Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a dedicated `ltkm3dc` workflow that reproduces `U_00`, `U_01`, `U_02`, and `U_03` and records their tolerance convergence against the tightest TKM run.

**Architecture:** Add a small helper file for tolerance validation and convergence-row generation so the data-shaping logic is testable without running the expensive orbital workflow in every test. Keep the physics path in a new script that mirrors the existing FBCPoisson setup, but swaps the potential evaluation to `TKM3D.ltkm3dc` and writes a separate CSV.

**Tech Stack:** Julia, BoundaryIntegral, TKM3D, CSV, DataFrames, Test

---

### Task 1: Add a testable helper layer for tolerance sweeps

**Files:**
- Create: `codes/volume_integral/tkm3d_hubbard_utils.jl`
- Create: `codes/volume_integral/test/runtests.jl`
- Create: `codes/volume_integral/test/tkm3d_hubbard_utils_test.jl`

**Step 1: Write the failing test**

```julia
using Test
include(joinpath(@__DIR__, "..", "tkm3d_hubbard_utils.jl"))

@testset "tolerance sweep rows use tightest run as reference" begin
    samples = [
        (tol = 1.0e-2, U_raw = 8.0, U_ev = 80.0),
        (tol = 1.0e-4, U_raw = 9.0, U_ev = 90.0),
    ]
    rows = annotate_reference_errors("U_01", samples)
    @test rows[end].rel_err_raw == 0.0
    @test rows[1].rel_err_raw == abs(8.0 - 9.0) / 9.0
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/runtests.jl`
Expected: FAIL because `tkm3d_hubbard_utils.jl` and `annotate_reference_errors` do not exist yet.

**Step 3: Write minimal implementation**

```julia
function annotate_reference_errors(pair, samples)
    reference = samples[end]
    return [(
        pair = pair,
        tol = row.tol,
        U_raw = row.U_raw,
        U_ev = row.U_ev,
        rel_err_raw = abs(row.U_raw - reference.U_raw) / abs(reference.U_raw),
        rel_err_ev = abs(row.U_ev - reference.U_ev) / abs(reference.U_ev),
    ) for row in samples]
end
```

**Step 4: Run test to verify it passes**

Run: `julia --project=. test/runtests.jl`
Expected: PASS with the helper test green.

**Step 5: Commit**

```bash
git add codes/volume_integral/tkm3d_hubbard_utils.jl codes/volume_integral/test/runtests.jl codes/volume_integral/test/tkm3d_hubbard_utils_test.jl
git commit -m "feat: add TKM3D Hubbard sweep helpers"
```

### Task 2: Add the dedicated `ltkm3dc` Hubbard convergence script

**Files:**
- Modify: `codes/volume_integral/tkm3d_hubbard_utils.jl`
- Create: `codes/volume_integral/compute_bare_hubbard_graphene_tkm3d.jl`
- Test: `codes/volume_integral/test/tkm3d_hubbard_utils_test.jl`

**Step 1: Write the failing test**

```julia
@testset "invalid tolerances are rejected" begin
    @test_throws ArgumentError normalize_tolerances([1.0e-2, 1.0, 1.0e-4])
end
```

**Step 2: Run test to verify it fails**

Run: `julia --project=. test/runtests.jl`
Expected: FAIL because `normalize_tolerances` is not implemented yet.

**Step 3: Write minimal implementation**

```julia
function normalize_tolerances(tolerances)
    values = sort!(collect(Float64.(tolerances)); rev = true)
    isempty(values) && throw(ArgumentError("tolerances must not be empty"))
    all(0.0 .< values .< 1.0) || throw(ArgumentError("tolerances must lie in (0, 1)"))
    return values
end
```

Then implement `compute_bare_hubbard_graphene_tkm3d.jl` to:
- reuse the current XSF loading and `VolumeSource` construction,
- evaluate each pair with `TKM3D.ltkm3dc`,
- call `estimate_kcut3dc` for metadata,
- use the helper functions to build the convergence rows,
- write `data/hubbard_graphene_tkm3d.csv`.

**Step 4: Run test to verify it passes**

Run: `julia --project=. test/runtests.jl`
Expected: PASS with both helper test sets green.

**Step 5: Commit**

```bash
git add codes/volume_integral/tkm3d_hubbard_utils.jl codes/volume_integral/compute_bare_hubbard_graphene_tkm3d.jl codes/volume_integral/test/tkm3d_hubbard_utils_test.jl
git commit -m "feat: add TKM3D Hubbard convergence script"
```

### Task 3: Verify the new workflow end-to-end

**Files:**
- Modify: `codes/volume_integral/data/hubbard_graphene_tkm3d.csv`

**Step 1: Run the new script on a short tolerance sweep**

```bash
julia --project=. compute_bare_hubbard_graphene_tkm3d.jl
```

The script should default to a small ordered tolerance set first while the workflow is being validated.

**Step 2: Inspect the output**

Check that:
- the CSV is created,
- all four pairs appear for each tolerance,
- the tightest tolerance row has `rel_err_raw = 0` and `rel_err_ev = 0` for each pair,
- `kcut` metadata is present.

**Step 3: Run the helper tests again**

Run: `julia --project=. test/runtests.jl`
Expected: PASS after the script changes.

**Step 4: Commit**

```bash
git add codes/volume_integral/data/hubbard_graphene_tkm3d.csv
git commit -m "data: add initial TKM3D Hubbard convergence output"
```

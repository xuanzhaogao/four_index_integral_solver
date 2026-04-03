# Bilayer Reference Plot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a one-off figure comparing the `eps_in = 2.4` bilayer-slab softmix data to the freestanding BLG Table I reference values.

**Architecture:** Keep the existing shared comparison module unchanged and add one dedicated plotting script for the bilayer case. Store the BLG reference values directly in that script and cover them with a narrow regression in the existing test file.

**Tech Stack:** Julia, CairoMakie, CSV.jl, DataFrames.jl, Test

---

### Task 1: Add a failing regression for the BLG reference values

**Files:**
- Modify: `test/test_screened_hubbard_comparison.jl`

- [ ] **Step 1: Write the failing test**

```julia
@testset "bilayer Table I reference values are defined" begin
    refs = Dict(
        "U_00" => ...,
        "U_01" => ...,
        "U_02" => ...,
        "U_03" => ...,
        "U_04" => ...,
        "U_05" => ...,
    )

    @test length(refs) == 6
end
```

- [ ] **Step 2: Run the focused test to verify it fails**

Run: `julia --project=. test/test_screened_hubbard_comparison.jl`
Expected: FAIL because the bilayer reference values are not yet exposed by any code path.

- [ ] **Step 3: Commit**

```bash
git add test/test_screened_hubbard_comparison.jl
git commit -m "test: cover bilayer reference plot values"
```

### Task 2: Add the one-off bilayer plot script

**Files:**
- Create: `scripts/plot_bilayer_hubbard_vs_reference.jl`
- Modify: `test/test_screened_hubbard_comparison.jl`

- [ ] **Step 1: Write the minimal implementation**

```julia
using CairoMakie
using CSV, DataFrames

const PAIR_ORDER = ["U_00", "U_01", "U_02", "U_03", "U_04", "U_05"]
const BLG_TABLE_I_EV = Dict(...)
```

Read `data/screened_hubbard_graphene.csv`, filter `mode == "SoftMixInversePermittivity"`, sort by `bandwidth`, and draw one curve plus one dashed reference line per pair. Save to `figs/bilayer_hubbard_vs_reference.png`.

- [ ] **Step 2: Run the focused test to verify it passes**

Run: `julia --project=. test/test_screened_hubbard_comparison.jl`
Expected: PASS

- [ ] **Step 3: Generate the figure**

Run: `julia --project=. scripts/plot_bilayer_hubbard_vs_reference.jl`
Expected: writes `figs/bilayer_hubbard_vs_reference.png`

- [ ] **Step 4: Commit**

```bash
git add scripts/plot_bilayer_hubbard_vs_reference.jl test/test_screened_hubbard_comparison.jl figs/bilayer_hubbard_vs_reference.png
git commit -m "feat: add bilayer reference comparison plot"
```

# Bare Monolayer k-mesh Sweep Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run the existing bare 4-channel pipeline on three monolayer-graphene Wannier datasets at different k-meshes (k_161601, k_252501, k_323201), parse CoQui reference values from each dataset's `_coqui_thc_crpa.out`, and produce a combined CSV + dated markdown report that tests Malte Rösner's Madelung-correction hypothesis.

**Architecture:** One new orchestration script driving the existing `compute_all_bare_channels` over three dataset directories. A small inline parser extracts the four "bare interactions (orbital-average)" lines from each CoQui output. The hypothesis-test analysis lives in the report markdown, not the script. No changes to `MonolayerOrbitalLoader.jl` or `MonolayerBareIntegrals.jl`.

**Tech Stack:** Julia 1.11 via juliaup (per `using-julia` rules), shared top-level `Project.toml`, `BoundaryIntegral.jl`, `CSV.jl`, `DataFrames.jl`. Reference XSFs + CoQui outputs at `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_*_nb_144_c_15/`.

---

## File Structure

Paths relative to `/mnt/home/xgao1/work/four_index_integral_solver/codes/graphene/`.

**New files:**

- `monolayer/scripts/bare_monolayer_kmesh_sweep.jl` — sweep driver: loops over three datasets, calls `compute_all_bare_channels`, parses each CoQui output, writes one CSV.
- `monolayer/data/bare_monolayer_kmesh_sweep.csv` — 12 rows (3 k-meshes × 4 channels), produced by the script.
- `monolayer/results/2026-05-13-bare-monolayer-kmesh-sweep.md` — narrative report with Madelung-hypothesis verdict.

**Unchanged but exercised:**

- `monolayer/src/MonolayerOrbitalLoader.jl` (`centered_monolayer_sources_padded`)
- `monolayer/src/MonolayerBareIntegrals.jl` (`compute_all_bare_channels`)

---

## CoQui output format (verified during planning)

Every `_coqui_thc_crpa.out` file in the three target directories contains exactly one block of the form:

```
bare interactions (orbital-average):
  - intra-orbital = (17.434191, 0.000000) eV
  - inter-orbital = (8.839615, 0.000000) eV
  - Hund's coupling (spin-flip) = (0.130804, 0.000000) eV
  - Hund's coupling (pair-hopping) = (0.130804, 0.000000) eV
static screened interactions (orbital-average):
```

The parser reads the four lines after `bare interactions (orbital-average):` and stops at `static screened`. The real part of each tuple `(real, imag)` is the channel value in eV. Mapping:

- `intra-orbital` → `:onsite`
- `inter-orbital` → `:nn`
- `Hund's coupling (spin-flip)` → `:hund_sf`
- `Hund's coupling (pair-hopping)` → `:hund_ph`

---

## Task 1: Inspect new datasets and confirm pipeline compatibility

**Files:** none (read-only verification)

- [ ] **Step 1: Confirm all three XSF pairs are readable and share the 150×150×192 / same-lattice metadata**

```bash
for d in k_161601_nb_144_c_15 k_252501_nb_144_c_15 k_323201_nb_144_c_15; do
  echo "=== $d ==="
  for f in graphene_00001.xsf graphene_00002.xsf _coqui_thc_crpa.out; do
    p=/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/$d/$f
    [ -r "$p" ] && echo "  OK   $f" || echo "  MISS $f"
  done
  grep -A 1 "DATAGRID_3D_UNKNOWN" /mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/$d/graphene_00001.xsf | tail -1
done
```

Expected output: every line `OK`; every grid-line shows `   150   150   192`. If any file is missing or unreadable, stop and report — the plan assumes all six XSFs and three `_coqui_thc_crpa.out` files are accessible.

- [ ] **Step 2: Confirm the orbital-average block exists in all three CoQui files**

```bash
for d in k_161601_nb_144_c_15 k_252501_nb_144_c_15 k_323201_nb_144_c_15; do
  echo "=== $d ==="
  grep -A 5 "bare interactions (orbital-average)" /mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/$d/_coqui_thc_crpa.out
done
```

Expected: each block shows four `- intra/inter/Hund's …` lines. Note down the four real-part values for each k-mesh — they are needed in the report's hypothesis-check table.

Reference values you should see (verified at planning time):

| k-mesh    | intra-orbital | inter-orbital | hund_sf  | hund_ph |
|-----------|--------------:|--------------:|---------:|--------:|
| k_161601  | 17.434191     | 8.839615      | 0.130804 | 0.130804 |
| k_252501  | 17.437588     | 8.837615      | 0.130581 | 0.130581 |
| k_323201  | 17.403608     | 8.831873      | 0.131551 | 0.131551 |

If the values printed by Step 2 differ from this table, **stop** — the CoQui files have changed since the plan was written; re-confirm with the user before proceeding.

---

## Task 2: Write the CoQui parser as a standalone, testable function

**Files:**
- Create: `monolayer/scripts/bare_monolayer_kmesh_sweep.jl` (parser only at this stage)

This is a one-off script, so we keep the parser inside the script file rather than extracting a module. The parser is a small, pure function we can exercise from the REPL.

- [ ] **Step 1: Create the script skeleton with just the parser**

Create `monolayer/scripts/bare_monolayer_kmesh_sweep.jl`:

```julia
# Sweep bare 4-channel integrals across three monolayer Wannier datasets
# (k_161601, k_252501, k_323201) and test Malte's Madelung-correction hypothesis.

using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "src", "MonolayerOrbitalLoader.jl"))
include(joinpath(@__DIR__, "..", "src", "MonolayerBareIntegrals.jl"))
using .MonolayerOrbitalLoader
using .MonolayerBareIntegrals

"""
    parse_coqui_bare(path)

Parse the "bare interactions (orbital-average)" block of a CoQui
`_coqui_thc_crpa.out` file. Returns a Dict mapping the four channel symbols
to their real-part values in eV.

Stops parsing at the first `static screened` line. Throws if any of the four
expected lines are missing.
"""
function parse_coqui_bare(path::AbstractString)
    isfile(path) || throw(ArgumentError("CoQui file not found: $path"))
    result = Dict{Symbol, Float64}()
    in_block = false
    for line in eachline(path)
        if occursin("bare interactions (orbital-average)", line)
            in_block = true
            continue
        end
        in_block || continue
        if occursin("static screened", line)
            break
        end
        m = match(r"-\s*(intra-orbital|inter-orbital|Hund's coupling \(spin-flip\)|Hund's coupling \(pair-hopping\))\s*=\s*\(\s*([-+0-9.eE]+)\s*,", line)
        m === nothing && continue
        label, val_str = m.captures[1], m.captures[2]
        key = label == "intra-orbital"                       ? :onsite   :
              label == "inter-orbital"                       ? :nn       :
              label == "Hund's coupling (spin-flip)"         ? :hund_sf  :
              label == "Hund's coupling (pair-hopping)"      ? :hund_ph  :
              error("unreachable: $label")
        result[key] = parse(Float64, val_str)
    end
    for k in (:onsite, :nn, :hund_sf, :hund_ph)
        haskey(result, k) || error("CoQui parser failed to find $k in $path")
    end
    return result
end
```

- [ ] **Step 2: Exercise the parser at the REPL to confirm it returns the expected values**

Run:

```bash
julia --project=. -e '
include("monolayer/scripts/bare_monolayer_kmesh_sweep.jl")
for d in ("k_161601_nb_144_c_15", "k_252501_nb_144_c_15", "k_323201_nb_144_c_15")
    p = joinpath("/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer", d, "_coqui_thc_crpa.out")
    println(d, " => ", parse_coqui_bare(p))
end
'
```

Expected output (values match the table in Task 1 Step 2):

```
k_161601_nb_144_c_15 => Dict(:onsite => 17.434191, :nn => 8.839615, :hund_sf => 0.130804, :hund_ph => 0.130804)
k_252501_nb_144_c_15 => Dict(:onsite => 17.437588, :nn => 8.837615, :hund_sf => 0.130581, :hund_ph => 0.130581)
k_323201_nb_144_c_15 => Dict(:onsite => 17.403608, :nn => 8.831873, :hund_sf => 0.131551, :hund_ph => 0.131551)
```

If any value differs from this expected output, fix the parser before proceeding. If `include` triggers a Julia precompile (the first invocation in a fresh session can take 30-60 s on this codebase), wait and re-run.

- [ ] **Step 3: Commit the parser scaffold**

```bash
git add monolayer/scripts/bare_monolayer_kmesh_sweep.jl
git commit -m "feat(monolayer): CoQui bare-interaction parser scaffold for k-mesh sweep"
```

---

## Task 3: Wire the parser into a full sweep with one dataset

Before scaling to three datasets, prove the end-to-end pipeline with one. We start with `k_252501` (a new dataset — exercises every code path including parsing).

**Files:**
- Modify: `monolayer/scripts/bare_monolayer_kmesh_sweep.jl`

- [ ] **Step 1: Append the sweep skeleton + Madelung table + a single-dataset run**

Append to `monolayer/scripts/bare_monolayer_kmesh_sweep.jl`:

```julia
const REF_ROOT = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer"
const SOURCE_TOL = 1e-3
const VOLUME_TOL = 1e-3

# Madelung corrections supplied by Malte (2026-05 email). Defined for the
# onsite channel only — the long-wavelength correction does not apply to the
# short-range nn / Hund's channels.
const MADELUNG_REF_EV = Dict(
    "k_161601_nb_144_c_15" =>  1.041,
    "k_252501_nb_144_c_15" =>  0.1842,
    "k_323201_nb_144_c_15" => -0.290,
)

const DATASETS = (
    "k_161601_nb_144_c_15",
    "k_252501_nb_144_c_15",
    "k_323201_nb_144_c_15",
)

const OUT_CSV = joinpath(@__DIR__, "..", "data", "bare_monolayer_kmesh_sweep.csv")

"""
    sweep_dataset(kmesh)

Run the 4-channel bare integral on a single dataset directory and combine
with parsed CoQui references. Returns a vector of NamedTuples (one per channel)
ready for assembly into the output CSV.
"""
function sweep_dataset(kmesh::AbstractString)
    dir = joinpath(REF_ROOT, kmesh)
    xsf_1 = joinpath(dir, "graphene_00001.xsf")
    xsf_2 = joinpath(dir, "graphene_00002.xsf")
    coqui = joinpath(dir, "_coqui_thc_crpa.out")
    isfile(xsf_1) || error("missing XSF 1 in $dir")
    isfile(xsf_2) || error("missing XSF 2 in $dir")
    isfile(coqui) || error("missing CoQui output in $dir")

    println("--- $kmesh ---")
    rows = compute_all_bare_channels(;
        orbital_1 = xsf_1, orbital_2 = xsf_2,
        source_tol = SOURCE_TOL, volume_tol = VOLUME_TOL,
        mirror_pad_level = 0,
    )
    coqui_vals = parse_coqui_bare(coqui)
    madelung_onsite = get(MADELUNG_REF_EV, kmesh, missing)

    out = NamedTuple[]
    for r in rows
        coqui_v = coqui_vals[r.channel]
        diff_ev = coqui_v - r.u_ev
        madelung_ref = r.channel === :onsite ? madelung_onsite : missing
        madelung_residual = r.channel === :onsite ? diff_ev - madelung_ref : missing
        push!(out, (
            kmesh = kmesh,
            channel = String(r.channel),
            our_u_ev = r.u_ev,
            coqui_u_ev = coqui_v,
            diff_ev = diff_ev,
            madelung_ref_ev = madelung_ref,
            madelung_residual_ev = madelung_residual,
            our_u_raw = r.u_raw,
            Na = r.Na,
            Nb = r.Nb,
            n_source_points = r.n_source_points,
            n_target_points = r.n_target_points,
            tkm_kmax = r.tkm_kmax,
            source_tol = r.source_tol,
            volume_tol = r.volume_tol,
            mirror_pad_level = r.mirror_pad_level,
            shift_x = r.shift_x,
            shift_y = r.shift_y,
            shift_z = r.shift_z,
        ))
    end
    return out
end
```

- [ ] **Step 2: Smoke-run `sweep_dataset` on one new dataset from the REPL**

Run:

```bash
julia --project=. -e '
include("monolayer/scripts/bare_monolayer_kmesh_sweep.jl")
rows = sweep_dataset("k_252501_nb_144_c_15")
for r in rows
    println(rpad(r.channel, 8), "  our=", r.our_u_ev, "  coqui=", r.coqui_u_ev, "  diff=", r.diff_ev, "  mad=", r.madelung_ref_ev, "  res=", r.madelung_residual_ev)
end
'
```

This call invokes `compute_all_bare_channels` and may take ~10-30 minutes (volume integration is the same cost as the 2026-04-23 run). Be patient. Expected output: four lines with `our_u_ev` close to (per Madelung prediction)

- `onsite ≈ 17.4376 − 0.184 = 17.253` eV
- `nn ≈ ?` (no strong prediction; somewhere near 8.7-8.8 eV)
- `hund_sf ≈ 0.13` eV, `hund_ph ≈ 0.13` eV

`madelung_residual_ev` on the onsite row should be small (target `|residual| ≲ 0.05` eV). If it's large (say > 0.3 eV), do not "fix" anything — record the value, continue, and discuss in the report. The job of this script is to measure, not to confirm.

- [ ] **Step 3: Commit the dataset-sweep machinery**

```bash
git add monolayer/scripts/bare_monolayer_kmesh_sweep.jl
git commit -m "feat(monolayer): per-dataset sweep with Madelung-residual diagnostic"
```

---

## Task 4: Drive all three datasets, write the CSV

**Files:**
- Modify: `monolayer/scripts/bare_monolayer_kmesh_sweep.jl`

- [ ] **Step 1: Append `main` and a top-level `main()` call**

Append to `monolayer/scripts/bare_monolayer_kmesh_sweep.jl`:

```julia
function main()
    all_rows = NamedTuple[]
    for kmesh in DATASETS
        rows = sweep_dataset(kmesh)
        append!(all_rows, rows)
        # Show per-dataset summary table immediately for easy mid-run progress
        df_local = DataFrame(rows)
        show(df_local[:, [:kmesh, :channel, :our_u_ev, :coqui_u_ev, :diff_ev, :madelung_ref_ev, :madelung_residual_ev]];
             allrows = true, allcols = true)
        println()
    end
    mkpath(dirname(OUT_CSV))
    table = DataFrame(all_rows)
    CSV.write(OUT_CSV, table)
    println("Wrote $(nrow(table)) rows to $(OUT_CSV)")
    show(table[:, [:kmesh, :channel, :our_u_ev, :coqui_u_ev, :diff_ev, :madelung_ref_ev, :madelung_residual_ev]];
         allrows = true, allcols = true)
    println()
end

main()
```

- [ ] **Step 2: Run the full sweep**

```bash
cd /mnt/home/xgao1/work/four_index_integral_solver/codes/graphene
julia --project=. monolayer/scripts/bare_monolayer_kmesh_sweep.jl
```

This will run three datasets × 4 channels and may take 30-90 minutes. Capture stdout to a log:

```bash
mkdir -p monolayer/logs
julia --project=. monolayer/scripts/bare_monolayer_kmesh_sweep.jl 2>&1 | tee monolayer/logs/kmesh_sweep_$(date +%Y%m%d-%H%M%S).log
```

Expected: at the end a 12-row table is printed and `monolayer/data/bare_monolayer_kmesh_sweep.csv` exists.

Sanity checks before declaring success:

1. CSV has 12 rows.
2. For each k-mesh, the onsite `our_u_ev` is in the range `[14, 19]` eV (loose). If any is wildly outside this, stop — the pipeline misbehaved.
3. `hund_sf` and `hund_ph` agree per row to ~1e-6 eV (real-orbital identity, same as 2026-04-23 verdict).
4. The `k_161601` onsite row's `our_u_ev` should be near 16.40 eV (within ~0.02 eV of the 2026-04-23 value 16.4004). A larger drift is worth a note but not a blocker.

- [ ] **Step 3: Inspect the CSV**

```bash
column -t -s, monolayer/data/bare_monolayer_kmesh_sweep.csv | cut -c1-160 | head -20
```

- [ ] **Step 4: Commit the script + CSV together**

```bash
git add monolayer/scripts/bare_monolayer_kmesh_sweep.jl monolayer/data/bare_monolayer_kmesh_sweep.csv
git commit -m "feat(monolayer): run 3-kmesh bare-channel sweep + CSV"
```

---

## Task 5: Write the dated report

**Files:**
- Create: `monolayer/results/2026-05-13-bare-monolayer-kmesh-sweep.md`

- [ ] **Step 1: Pull the 12 numbers from the CSV**

Either read `monolayer/data/bare_monolayer_kmesh_sweep.csv` directly, or run:

```bash
julia --project=. -e '
using CSV, DataFrames
df = CSV.read("monolayer/data/bare_monolayer_kmesh_sweep.csv", DataFrame)
show(df[:, [:kmesh, :channel, :our_u_ev, :coqui_u_ev, :diff_ev, :madelung_ref_ev, :madelung_residual_ev]]; allrows=true, allcols=true)
println()
'
```

- [ ] **Step 2: Create the report**

Create `monolayer/results/2026-05-13-bare-monolayer-kmesh-sweep.md`. Replace every `<fill>` with the actual number from the CSV. Compute `verdict` per the rules in the spec (`|residual| ≤ 0.05 eV` for all three = Confirmed; `≤ 0.1 eV` = Partial; otherwise Rejected).

```markdown
# Bare Monolayer Graphene — k-mesh Sweep & Madelung Verification

Date: 2026-05-13

## Context

Following the 2026-04-23 bare-monolayer report, Malte Rösner identified that
the 1.034 eV gap between our direct Wannier integral (16.400 eV) and the CoQui
cRPA bare value (17.434 eV) at the original `k_161601` sampling is very close
to CoQui's long-wavelength Madelung correction for that k-mesh (1.041 eV).

To test the hypothesis

    CoQui_v0  =  direct_Wannier_v0  +  Madelung_correction(k-mesh)

Malte provided two additional datasets at different k/q samplings (= different
effective real-space supercell sizes) with sign-changing Madelung corrections.
Reference data:

- `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_161601_nb_144_c_15`  (Madelung = +1.041 eV)
- `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_252501_nb_144_c_15`  (Madelung = +0.1842 eV)
- `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15`  (Madelung = −0.290 eV)

## Method

Identical pipeline to the 2026-04-23 report:

- Loader: `centered_monolayer_sources_padded` with `mirror_pad_level = 0`.
- Integrator: `BI.TKM3D.ltkm3dc` direct `1/|r-r'|` over `BI.VolumeSource` quadratures.
- Tolerances: `source_tol = 1e-3`, `volume_tol = 1e-3`.
- CoQui reference: parsed from the "bare interactions (orbital-average)" block
  of each dataset's `_coqui_thc_crpa.out`.

All three datasets share the 150×150×192 grid on the same lattice; only the
underlying Wannier orbitals differ (different k-meshes → different MLWF).

## Results

### Full 12-row table

| k-mesh    | channel  | our `u_ev` | CoQui `u_ev` | `CoQui − ours` |
|-----------|----------|-----------:|-------------:|---------------:|
| k_161601  | onsite   | <fill>     | 17.434191    | <fill>         |
| k_161601  | nn       | <fill>     |  8.839615    | <fill>         |
| k_161601  | hund_sf  | <fill>     |  0.130804    | <fill>         |
| k_161601  | hund_ph  | <fill>     |  0.130804    | <fill>         |
| k_252501  | onsite   | <fill>     | 17.437588    | <fill>         |
| k_252501  | nn       | <fill>     |  8.837615    | <fill>         |
| k_252501  | hund_sf  | <fill>     |  0.130581    | <fill>         |
| k_252501  | hund_ph  | <fill>     |  0.130581    | <fill>         |
| k_323201  | onsite   | <fill>     | 17.403608    | <fill>         |
| k_323201  | nn       | <fill>     |  8.831873    | <fill>         |
| k_323201  | hund_sf  | <fill>     |  0.131551    | <fill>         |
| k_323201  | hund_ph  | <fill>     |  0.131551    | <fill>         |

### Madelung hypothesis test (onsite channel only)

| k-mesh    | our `v0` | CoQui `v0` | `CoQui − ours` | Madelung ref | residual |
|-----------|---------:|-----------:|---------------:|-------------:|---------:|
| k_161601  | <fill>   | 17.434191  | <fill>         |  1.0410      | <fill>   |
| k_252501  | <fill>   | 17.437588  | <fill>         |  0.1842      | <fill>   |
| k_323201  | <fill>   | 17.403608  | <fill>         | −0.2900      | <fill>   |

Verdict thresholds (from the spec):
- Confirmed: `|residual| ≤ 0.05 eV` for all three rows
- Partial:   `|residual| ≤ 0.10 eV` for all three rows
- Rejected:  otherwise

**Verdict:** <fill — one of Confirmed / Partial / Rejected, with one sentence pointing at the largest residual>.

## Reproducibility

- Source data: paths above (CoQui directories at `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/`).
- Driver script: `monolayer/scripts/bare_monolayer_kmesh_sweep.jl`.
- Raw CSV:       `monolayer/data/bare_monolayer_kmesh_sweep.csv`.
- Run log:       `monolayer/logs/kmesh_sweep_*.log`.

## Notes

- The k_161601 onsite value was re-computed in this sweep and equals <fill>, compared with 16.4004 from the 2026-04-23 CSV. Drift = <fill> eV (expected to be small; any drift > 0.01 eV is itself worth a note).
- `hund_sf` and `hund_ph` agree to machine precision in every row, reproducing the real-orbital identity `V_abba = V_aabb` already established at 2026-04-23.
- Mirror-pad sweep was deliberately skipped (the 2026-04-23 report concluded pad-1 is inconclusive on these XSFs).

## Next

- If verdict is **Confirmed** or **Partial**: write back to Malte with the residual table and ask whether reproducing his Madelung correction in our code (Gygi-Baldereschi style) is in scope. The XSF resolution question (3×3×1 regeneration) is moot if the residual is already at the Madelung-formula level of precision.
- If verdict is **Rejected**: investigate (1) whether the k-mesh changes are actually changing the Wannier function shape as much as expected, (2) whether the CoQui Madelung convention matches what was reported (sign, normalization). Defer the XSF-resolution question until this is sorted.
- The bilayer-monolayer follow-up Malte mentioned remains blocked on him producing a starting example.
```

- [ ] **Step 3: Fill in every `<fill>` from the CSV**

Open the file, edit each `<fill>` cell with the actual number from `monolayer/data/bare_monolayer_kmesh_sweep.csv`. Compute the verdict per the threshold rules.

- [ ] **Step 4: Commit the report**

```bash
git add monolayer/results/2026-05-13-bare-monolayer-kmesh-sweep.md
git commit -m "docs(monolayer): k-mesh sweep results + Madelung-hypothesis verdict"
```

---

## Task 6: Final verification

- [ ] **Step 1: Working tree clean, three new files staged**

```bash
git status
git log --oneline -10
```

Expected: working tree clean; four new commits on `main` (parser scaffold, per-dataset sweep, full-sweep CSV, report).

- [ ] **Step 2: Re-run the test suite to confirm no regressions in the existing modules**

```bash
julia --project=. test/runtests.jl
```

Expected: all `bilayer_slab` and `monolayer` tests still pass — we did not touch `MonolayerOrbitalLoader.jl` or `MonolayerBareIntegrals.jl`.

- [ ] **Step 3: Quick re-derivation check on the CSV**

```bash
julia --project=. -e '
using CSV, DataFrames
df = CSV.read("monolayer/data/bare_monolayer_kmesh_sweep.csv", DataFrame)
@assert nrow(df) == 12
for row in eachrow(df[df.channel .== "onsite", :])
    if !ismissing(row.madelung_residual_ev)
        println(row.kmesh, "  residual = ", row.madelung_residual_ev, " eV")
    end
end
'
```

Expected: three lines printed, one per k-mesh, each showing the onsite Madelung residual that drives the report verdict.

---

## Self-Review

**Spec coverage:**
- Spec § "Scope → In scope" — run 3 datasets × 4 channels at `pad=0`: Task 3 + Task 4.
- Spec § "Scope → In scope" — parse CoQui from `_coqui_*_crpa*.out`: Task 2 (parser) + Task 1 Step 2 (format check). Note: we converged on `_coqui_thc_crpa.out` (orbital-average summary, readable for all three datasets) as the single parsing target; `_coqui_crpa_loc.out` is not used. This is a tightening of the spec, not a divergence — the spec said "try each in order", and one suffices.
- Spec § "CSV schema" — all columns in Task 3 Step 1 `push!`.
- Spec § "Madelung-correction comparison framework" → covered by `madelung_ref_ev` and `madelung_residual_ev` columns in the CSV and the dedicated table in the report.
- Spec § "Output → New dated report + CSV": Task 4 (CSV) and Task 5 (report).

**Out of scope items remain unimplemented (intentional):**
- Mirror-pad sweep on new datasets.
- Madelung correction in our own code.
- Reply to Malte about XSF resolution (left to the "Next" section of the report, post-verdict).

**Placeholders:** `<fill>` markers appear only in the report template, which Task 5 Step 3 explicitly populates from the CSV before commit. No `TBD`/`TODO` left in script code.

**Type consistency:** `r.channel` is a `Symbol` from `bare_channel_integral`; we lookup `coqui_vals[r.channel]` (Dict keyed by symbol) and convert to `String` only for the CSV column. The Madelung table uses string keys, looked up by `kmesh` (also string). No naming drift between tasks.

**Files mentioned only once vs. referenced consistently:** `monolayer/scripts/bare_monolayer_kmesh_sweep.jl`, `monolayer/data/bare_monolayer_kmesh_sweep.csv`, `monolayer/results/2026-05-13-bare-monolayer-kmesh-sweep.md` — all three appear in the file-structure header, in the steps that create/use them, and in the report's "Reproducibility" section.

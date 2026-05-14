# Screened Monolayer LZ Sweep Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Find a slab thickness `LZ` (at fixed `eps_in = 2.4`) that drives the screened density channels (onsite, nn) within 5% of CoQui's `U_ijkl` at the k_323201 monolayer dataset, by running the bilayer dielectric solver at `LZ ∈ {2.0, 3.0, 3.35, 5.0, 8.0}` Å.

**Architecture:** Black-box reuse of the just-completed `monolayer/src/MonolayerScreenedSolve.jl` (parser, source builder, target_specs, run_record) and `bilayer_slab/src/ScreenedOrbitalSolve.jl` (`solve_screened_mode`). A new ~80-line driver loops over LZ values, reconstructs sources with `Z_CENTER = LZ/2` per LZ, runs **one density Sharp solve per LZ** (no soft modes, no Hund's), parses CoQui once, writes a 10-row CSV. A new Slurm submit script mirrors the k_323201 one with a shorter wall request.

**Tech Stack:** Julia 1.11 via juliaup, shared `Project.toml`. Rusty `ccm` / genoa, 48 cores, ~10–15 min wall. XSFs + CoQui ref at `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15/`.

---

## File Structure

Paths relative to `/mnt/home/xgao1/work/four_index_integral_solver/codes/graphene/`.

**New files:**

- `monolayer/scripts/screened_monolayer_lz_sweep.jl` — Julia driver (~80 lines).
- `monolayer/scripts/submit_screened_monolayer_lz_sweep.slurm` — Slurm batch script for Rusty `ccm`/genoa, 48 cores, 2 h wall.
- `monolayer/data/screened_monolayer_lz_sweep.csv` — 10 rows (5 LZ × 2 channels). Produced by the Slurm job.
- `monolayer/results/2026-05-13-screened-monolayer-lz-sweep-report.md` — dated report.

**Reused unchanged (black box):**

- `monolayer/src/MonolayerScreenedSolve.jl` — `parse_coqui_loc`, `monolayer_screened_sources`, `density_target_specs`, `screened_run_record`.
- `bilayer_slab/src/ScreenedOrbitalSolve.jl` — `solve_screened_mode`.
- `monolayer/src/MonolayerOrbitalLoader.jl` — orbital loaders.

No test additions needed — the helpers are already covered by `monolayer/test/test_monolayer_screened_solve.jl` (45 assertions). This plan only adds a driver script + submit file + analysis report.

---

## Task 1: Confirm prerequisites

**Files:** none (read-only)

- [ ] **Step 1: Verify the existing k_323201 work is in place**

```bash
cd /mnt/home/xgao1/work/four_index_integral_solver/codes/graphene
ls -1 monolayer/scripts/screened_monolayer_k323201.jl \
      monolayer/scripts/submit_screened_monolayer_k323201.slurm \
      monolayer/src/MonolayerScreenedSolve.jl \
      monolayer/data/screened_monolayer_k323201.csv \
      monolayer/results/2026-05-13-screened-monolayer-k323201-report.md
```

Expected: all 5 files listed (the k_323201 baseline). If any is missing, stop and report BLOCKED — the LZ sweep depends on those artifacts already existing.

- [ ] **Step 2: Confirm reference k_323201 CoQui values are still readable**

```bash
grep -E "^\s*0\s+(0\s+0\s+0|0\s+1\s+1)" /mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15/_coqui_crpa_loc.out
```

Expected:

```
  0  0  0  0   17.4029    9.7824
  0  0  1  1    8.8319    5.1678
```

If different, stop.

---

## Task 2: Write the LZ-sweep driver

**Files:**
- Create: `monolayer/scripts/screened_monolayer_lz_sweep.jl`

The driver mirrors the structure of `screened_monolayer_k323201.jl` but with a `LZ_VALUES` loop in `main()`, no SoftMix modes, no Hund's solve.

- [ ] **Step 1: Create the driver**

Create `monolayer/scripts/screened_monolayer_lz_sweep.jl`:

```julia
using BoundaryIntegral
import BoundaryIntegral as BI

using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "..", "bilayer_slab", "src", "ScreenedOrbitalSolve.jl"))
using .ScreenedOrbitalSolve

include(joinpath(@__DIR__, "..", "src", "MonolayerScreenedSolve.jl"))
using .MonolayerScreenedSolve

const REF_DIR = "/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15"
const XSF_1   = joinpath(REF_DIR, "graphene_00001.xsf")
const XSF_2   = joinpath(REF_DIR, "graphene_00002.xsf")
const COQUI   = joinpath(REF_DIR, "_coqui_crpa_loc.out")

# Slab dielectric model (fixed)
const EPS_IN  = 2.4
const EPS_OUT = 1.0
const L       = 90.0

# LZ values to sweep (Å). Brackets the original 3.35 baseline on both sides.
const LZ_VALUES = [2.0, 3.0, 3.35, 5.0, 8.0]

# Solver tolerances (reused from k_323201 protocol)
const SOURCE_TOL = 1e-3
const N_QUAD = 6
const EDGE_REFINE_LEVEL = 4
const RHS_TOL = 1e-3
const LHS_TOL = 1e-5
const GMRES_ATOL = 1e-5
const GMRES_RTOL = 1e-5
const MAX_ORDER = 64
const MAX_DEPTH = 12

const OUT_CSV = joinpath(@__DIR__, "..", "data", "screened_monolayer_lz_sweep.csv")

function _solve_kwargs()
    return (
        n_quad = N_QUAD,
        edge_refine_level = EDGE_REFINE_LEVEL,
        rhs_tol = RHS_TOL,
        lhs_tol = LHS_TOL,
        gmres_atol = GMRES_ATOL,
        gmres_rtol = GMRES_RTOL,
        max_order = MAX_ORDER,
        max_depth = MAX_DEPTH,
        volume_tol = RHS_TOL,
    )
end

function _record_kwargs(LZ, z_center)
    return (
        mode_label = "Sharp", bandwidth = missing,
        eps_in = EPS_IN, eps_out = EPS_OUT, Lz = LZ, L = L, z_center = z_center,
        source_tol = SOURCE_TOL, rhs_tol = RHS_TOL, lhs_tol = LHS_TOL,
        gmres_atol = GMRES_ATOL, gmres_rtol = GMRES_RTOL,
        n_quad = N_QUAD, edge_refine_level = EDGE_REFINE_LEVEL,
        max_order = MAX_ORDER, max_depth = MAX_DEPTH,
    )
end

function main()
    println("Parsing CoQui reference ...")
    coqui = parse_coqui_loc(COQUI)
    println("  onsite U_ijkl = ", coqui[:onsite].U_ijkl)
    println("  nn     U_ijkl = ", coqui[:nn].U_ijkl)

    rows = NamedTuple[]

    for LZ in LZ_VALUES
        z_center = LZ / 2
        println("--- LZ = ", LZ, " Å  (Z_CENTER = ", z_center, ") ---")
        src = monolayer_screened_sources(
            orbital_1 = XSF_1, orbital_2 = XSF_2,
            source_tol = SOURCE_TOL, z_center = z_center,
        )
        println("  Nphi1 = ", src.Nphi1, "  Nphi2 = ", src.Nphi2)

        density_specs = density_target_specs(src)
        result = solve_screened_mode(
            src.vs1, density_specs,
            L, L, LZ, EPS_IN, EPS_OUT, BI.SharpScreening();
            _solve_kwargs()...,
        )
        println("  residual = ", result.residual,
                "  n_interface_points = ", result.n_interface_points)

        for pr in result.pair_results
            ch_sym = pr.pair == "onsite" ? :onsite : :nn
            row0 = screened_run_record(pr, result; _record_kwargs(LZ, z_center)...)
            U_ref = coqui[ch_sym].U_ijkl
            row = merge(row0, (
                LZ = LZ,
                coqui_v_ijkl = coqui[ch_sym].v_ijkl,
                coqui_U_ijkl = U_ref,
                diff_ev = U_ref - row0.u_total_ev,
                rel_err_pct = 100 * (row0.u_total_ev - U_ref) / U_ref,
            ))
            push!(rows, row)
            println("    ", pr.pair, "  ours=", row0.u_total_ev,
                    "  CoQui=", U_ref,
                    "  rel_err_pct=", row.rel_err_pct)
        end
    end

    mkpath(dirname(OUT_CSV))
    table = DataFrame(rows)
    CSV.write(OUT_CSV, table)
    println("Wrote $(nrow(table)) rows to $(OUT_CSV)")
    summary_cols = [:LZ, :channel, :u_int_ev, :u_scatter_ev, :u_total_ev, :coqui_U_ijkl, :rel_err_pct, :sigma_residual]
    show(table[:, summary_cols]; allrows = true, allcols = true)
    println()
end

main()
```

- [ ] **Step 2: Quick syntax / load test**

This does NOT run any solves — just confirms the file parses, modules load, and `main` is defined. Should complete in <5 s after precompile.

```bash
julia --project=. -e '
include("monolayer/scripts/screened_monolayer_lz_sweep.jl") |> _ -> nothing
' 2>&1 | head -1
```

This will actually call `main()` because of the top-level `main()` invocation. To avoid running the full sweep at this step, instead just check that the file compiles cleanly:

```bash
julia --project=. -e '
contents = read("monolayer/scripts/screened_monolayer_lz_sweep.jl", String)
# Comment out the bottom main() call before evaluating, just to load the funcs.
import Base.Meta
parse_ok = try
    Meta.parseall(contents)
    true
catch err
    println("PARSE ERROR: ", err)
    false
end
println("parse_ok = ", parse_ok)
'
```

Expected: `parse_ok = true`.

- [ ] **Step 3: Commit the driver (no run yet)**

```bash
git add monolayer/scripts/screened_monolayer_lz_sweep.jl
git commit -m "$(cat <<'EOF'
feat(monolayer): LZ-sweep driver at fixed eps_in=2.4 (Sharp density only)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Write the Slurm submit script

**Files:**
- Create: `monolayer/scripts/submit_screened_monolayer_lz_sweep.slurm`

- [ ] **Step 1: Create the submit script**

Create `monolayer/scripts/submit_screened_monolayer_lz_sweep.slurm`:

```bash
#!/bin/bash
#
# Slurm submission script for the screened-monolayer LZ sweep at fixed eps_in=2.4.
# Submit from the graphene project root (path containing Project.toml):
#
#     cd /mnt/home/xgao1/work/four_index_integral_solver/codes/graphene
#     sbatch monolayer/scripts/submit_screened_monolayer_lz_sweep.slurm
#
# Wall-clock estimate: ~10-15 min for 5 Sharp density solves (the k_323201 sweep
# did 6 mode invocations in 6m33s). 2-hour limit gives generous margin in case
# a larger LZ blows up the adaptive-interface point count.

#SBATCH --job-name=lz_sweep_k323201
#SBATCH --output=monolayer/logs/lz_sweep_%j.out
#SBATCH --error=monolayer/logs/lz_sweep_%j.err
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=02:00:00

set -euo pipefail

# Threading: BoundaryIntegral's FMM3D and BLAS use OpenMP; Julia threading is
# also enabled in case any kernel path uses Threads.@threads.
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}

mkdir -p monolayer/logs

# Use the juliaup-managed Julia. Do NOT `module load julia`.
echo "host: $(hostname)"
echo "date: $(date -Is)"
echo "julia: $(which julia)"
julia --version
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS JULIA_NUM_THREADS=$JULIA_NUM_THREADS"
echo ""

srun --cpus-per-task=$SLURM_CPUS_PER_TASK \
    julia --project=. monolayer/scripts/screened_monolayer_lz_sweep.jl

echo ""
echo "completed at: $(date -Is)"
```

- [ ] **Step 2: Commit**

```bash
git add monolayer/scripts/submit_screened_monolayer_lz_sweep.slurm
git commit -m "$(cat <<'EOF'
feat(monolayer): Slurm submit script for LZ sweep

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Submit and wait for the cluster job

The Slurm `sbatch` command must be executed by the user, per the organization's cluster-policy LAW (controller and agents may not execute non-read-only Slurm commands).

- [ ] **Step 1: User submits**

Ask the user to run, in this session via the `!` prefix or in their own terminal:

```
! sbatch monolayer/scripts/submit_screened_monolayer_lz_sweep.slurm
```

Capture the job ID printed (e.g., `Submitted batch job XXXXXXX`).

- [ ] **Step 2: Monitor**

After submission, check status with read-only Slurm commands (these are allowed):

```bash
JOBID=<paste job id here>
squeue -j "$JOBID"
sacct -j "$JOBID" --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

Wait for `State = COMPLETED` and `ExitCode = 0:0`. Do not poll repeatedly — check once after a 10–15 min interval, or wait for the user to confirm completion.

- [ ] **Step 3: Verify outputs**

```bash
ls -la monolayer/logs/lz_sweep_${JOBID}.{out,err}
ls -la monolayer/data/screened_monolayer_lz_sweep.csv
tail -30 monolayer/logs/lz_sweep_${JOBID}.out
```

Expected:
- `.out` ends with "Wrote 10 rows to .../screened_monolayer_lz_sweep.csv" and a 10-row DataFrame summary.
- CSV exists, has 11 lines (header + 10 data rows).

Quick sanity check via Julia:

```bash
julia --project=. -e '
using CSV, DataFrames
df = CSV.read("monolayer/data/screened_monolayer_lz_sweep.csv", DataFrame)
@assert nrow(df) == 10 "expected 10 rows, got $(nrow(df))"
@assert all(df.sigma_residual .< 1e-3) "some sigma_residual not converged"
println("OK: ", nrow(df), " rows, max residual = ", maximum(df.sigma_residual))
show(df[:, [:LZ, :channel, :u_total_ev, :coqui_U_ijkl, :rel_err_pct, :sigma_residual]]; allrows=true, allcols=true)
println()
'
```

If any GMRES residual > 1e-3, note which LZ(s) and flag for the report. If `nrow ≠ 10`, stop and investigate.

- [ ] **Step 4: Commit the CSV**

```bash
git add monolayer/data/screened_monolayer_lz_sweep.csv
git commit -m "$(cat <<'EOF'
data(monolayer): LZ-sweep CSV — 5 LZ × 2 channels at eps_in=2.4

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Write the dated report

**Files:**
- Create: `monolayer/results/2026-05-13-screened-monolayer-lz-sweep-report.md`

- [ ] **Step 1: Extract the numbers from the CSV**

```bash
julia --project=. -e '
using CSV, DataFrames, Printf
df = CSV.read("monolayer/data/screened_monolayer_lz_sweep.csv", DataFrame)
println("LZ-vs-rel_err table:")
println(rpad("LZ", 6), rpad("onsite ours", 14), rpad("onsite rel_err%", 18),
        rpad("nn ours", 12), rpad("nn rel_err%", 14))
for LZ in sort!(unique(df.LZ))
    on = df[(df.LZ .== LZ) .& (df.channel .== "onsite"), :][1, :]
    nn = df[(df.LZ .== LZ) .& (df.channel .== "nn"),     :][1, :]
    @printf("%-6.2f%-14.4f%-18.3f%-12.4f%-14.3f\n",
            LZ, on.u_total_ev, on.rel_err_pct, nn.u_total_ev, nn.rel_err_pct)
end
println()
println("Best LZ by max(|rel_err_onsite|, |rel_err_nn|):")
best_LZ = nothing
best_max = Inf
for LZ in sort!(unique(df.LZ))
    on = df[(df.LZ .== LZ) .& (df.channel .== "onsite"), :][1, :]
    nn = df[(df.LZ .== LZ) .& (df.channel .== "nn"),     :][1, :]
    m = max(abs(on.rel_err_pct), abs(nn.rel_err_pct))
    println("  LZ=", LZ, "  max_abs_rel_err=", m, "%")
    if m < best_max
        best_max = m
        best_LZ = LZ
    end
end
println("  best LZ = ", best_LZ, " with max_abs_rel_err = ", best_max, "%")
'
```

Record the best LZ value and its max-abs-rel-err. Determine verdict per spec:
- Strong: best max ≤ 5%
- Acceptable: best max ≤ 15%
- Disagreement: otherwise

- [ ] **Step 2: Create the report**

Create `monolayer/results/2026-05-13-screened-monolayer-lz-sweep-report.md`. Replace every `<fill>` with the actual numbers; choose the verdict per the criteria.

```markdown
# Screened Monolayer Graphene — LZ Sweep at eps_in=2.4 — Report

Date: 2026-05-13

## Context

The 2026-05-13 screened-monolayer report at k_323201 (`LZ = 3.35 Å, eps_in = 2.4`) showed +20% over-estimation on density channels (`U_onsite`, `U_nn`) while the Hund's channel agreed within 4%. The user's directive: keep `eps_in = 2.4` fixed and sweep `LZ` instead, on the hypothesis that the slab thickness controls how much of the orbital tail sits inside the screened region. This report runs that sweep.

Reference: `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15/_coqui_crpa_loc.out`.

| channel | U_ijkl (eV) |
|---------|------------:|
| onsite  | 9.7824      |
| nn      | 5.1678      |

## Method

| parameter | value |
|-----------|------:|
| eps_in    | 2.4 (fixed) |
| eps_out   | 1.0 (fixed) |
| LZ values | 2.0, 3.0, 3.35, 5.0, 8.0 Å |
| Z_CENTER  | LZ/2 (symmetric) |
| L (in-plane) | 90 Å |
| mode | Sharp only |
| solver tolerances | identical to the 2026-05-13 k_323201 run |

Per-LZ procedure: rebuild source VolumeSources with `z_center = LZ/2`, run one Sharp density solve, evaluate onsite (target=vs1) and nn (target=vs2). Hund's channel skipped — already in Strong agreement at the baseline LZ=3.35.

Compute: Rusty `ccm`/genoa, 48 cores, Slurm job `<fill: jobid>`, `<fill: elapsed>` wall, `<fill: maxRSS>` peak RSS.

## Results

| LZ (Å) | onsite ours | onsite rel_err% | nn ours | nn rel_err% | u_scatter_onsite | u_scatter_nn |
|--------|------------:|----------------:|--------:|------------:|-----------------:|-------------:|
|  2.0   | <fill>      | <fill>          | <fill>  | <fill>      | <fill>           | <fill>       |
|  3.0   | <fill>      | <fill>          | <fill>  | <fill>      | <fill>           | <fill>       |
|  3.35  | <fill>      | <fill>          | <fill>  | <fill>      | <fill>           | <fill>       |
|  5.0   | <fill>      | <fill>          | <fill>  | <fill>      | <fill>           | <fill>       |
|  8.0   | <fill>      | <fill>          | <fill>  | <fill>      | <fill>           | <fill>       |

`u_total = u_int + u_scatter`. `rel_err = 100 * (ours - CoQui) / CoQui`.

## Verdict

Best LZ by `max(|rel_err_onsite|, |rel_err_nn|)`: **LZ = <fill> Å** with max-abs-rel-err **<fill>%**.

Spec thresholds:
- Strong: best max ≤ 5%
- Acceptable: best max ≤ 15%
- Disagreement: otherwise

**Verdict: <fill>**

## Interpretation

<fill: a few sentences describing the trends. Specifically address:
  1. Does rel_err monotonically decrease with LZ, or is there a minimum?
  2. Do onsite and nn agree on the best LZ, or do they prefer different values?
  3. Is the u_scatter (interface) term getting more negative (more screening) at larger LZ, or less?
  4. If the best LZ is at the edge of [2.0, 8.0], should the sweep be widened in a follow-up?
>

## Reproducibility

- Driver:        `monolayer/scripts/screened_monolayer_lz_sweep.jl`
- Submit script: `monolayer/scripts/submit_screened_monolayer_lz_sweep.slurm`
- CSV:           `monolayer/data/screened_monolayer_lz_sweep.csv`
- Run logs:      `monolayer/logs/lz_sweep_<fill: jobid>.{out,err}`

## Next

<fill: depending on the verdict:
  - If Strong: declare LZ calibrated at the best value; replicate at k_252501 and k_161601 with the same LZ.
  - If Acceptable: open the eps_in sweep at the best LZ to push remaining error inside 5%.
  - If Disagreement and best LZ is at the boundary: widen the sweep.
  - If Disagreement and minimum is interior: slab model structurally limited — recommend trying an eps(z) profile or a 2D-Lindhard correction.
>
```

- [ ] **Step 3: Commit**

```bash
git add monolayer/results/2026-05-13-screened-monolayer-lz-sweep-report.md
git commit -m "$(cat <<'EOF'
docs(monolayer): LZ-sweep results + verdict

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Final verification

- [ ] **Step 1: Working tree clean, expected commits in place**

```bash
git status
git log --oneline -8
```

Expected commits (newest first):

- `docs(monolayer): LZ-sweep results + verdict`
- `data(monolayer): LZ-sweep CSV — 5 LZ × 2 channels at eps_in=2.4`
- `feat(monolayer): Slurm submit script for LZ sweep`
- `feat(monolayer): LZ-sweep driver at fixed eps_in=2.4 (Sharp density only)`
- `spec: LZ sweep ...` (pre-existing)

- [ ] **Step 2: Test suite still green**

```bash
julia --project=. test/runtests.jl
```

Bash timeout 600000 ms. Expected: all bilayer_slab + monolayer tests pass (189 total from the prior baseline; this plan adds no tests so the number is unchanged).

- [ ] **Step 3: CSV sanity**

```bash
julia --project=. -e '
using CSV, DataFrames
df = CSV.read("monolayer/data/screened_monolayer_lz_sweep.csv", DataFrame)
@assert nrow(df) == 10
@assert all(df.sigma_residual .< 1e-3)
@assert sort!(unique(df.LZ)) == [2.0, 3.0, 3.35, 5.0, 8.0]
println("OK: ", nrow(df), " rows, all sigma_residual < 1e-3")
'
```

---

## Self-Review

**Spec coverage:**
- Spec § "Scope → In scope" — 5 LZ values, Sharp only, 2 channels each, CoQui parse-once: Task 2 driver covers this with `LZ_VALUES` const and the single mode call.
- Spec § "File layout" — `screened_monolayer_lz_sweep.jl` (Task 2), `submit_screened_monolayer_lz_sweep.slurm` (Task 3), CSV (Task 4 output), report (Task 5).
- Spec § "Slab parameters" — encoded as `const` declarations in the driver.
- Spec § "CSV schema" — Task 2's `push!(rows, merge(row0, (LZ=…, coqui_v_ijkl=…, coqui_U_ijkl=…, diff_ev=…, rel_err_pct=…)))` produces all required columns (the `row0` from `screened_run_record` carries the 28 standard columns; we add LZ + 4 CoQui-comparison columns on top).
- Spec § "Verdict criteria" — encoded in the Task 5 report template and the verdict-selection code in Task 5 Step 1.
- Spec § "Risks / loose ends" — risks #1 (small-LZ tail) handled by reporting whatever residual we see and dropping non-converged rows from the verdict; risk #2 (n_interface_points scaling) absorbed by the 2-hour Slurm wall.

**Placeholder scan:**
- `<fill>` markers appear only in the report template (Task 5 Step 2). Task 5 Step 1 (CSV extraction) feeds the values that populate them. No `TBD`/`TODO`/"appropriate"/"handle edge cases" anywhere in script code.

**Type consistency:**
- `parse_coqui_loc` → `Dict{Symbol, NamedTuple{(:v_ijkl, :U_ijkl), Tuple{Float64, Float64}}}` (same as the k_323201 plan); consumed as `coqui[:onsite]` / `coqui[:nn]` in Task 2 main loop. ✓
- `monolayer_screened_sources(; orbital_1, orbital_2, source_tol, z_center)` keyword names match the existing module's signature. ✓
- `screened_run_record(pair_result, solve_result; mode_label, bandwidth, eps_in, eps_out, Lz, L, z_center, source_tol, rhs_tol, lhs_tol, gmres_atol, gmres_rtol, n_quad, edge_refine_level, max_order, max_depth)` — `_record_kwargs(LZ, z_center)` produces exactly these 17 keys. ✓
- The CSV column merge in Task 2 adds `LZ, coqui_v_ijkl, coqui_U_ijkl, diff_ev, rel_err_pct` on top of the 28-field record. Task 5 Step 1 reads `LZ`, `channel`, `u_total_ev`, `rel_err_pct`, `sigma_residual` — all defined. ✓
- The Slurm submit script template mirrors the existing `submit_screened_monolayer_k323201.slurm` line-for-line apart from job name, output paths, and wall time. ✓

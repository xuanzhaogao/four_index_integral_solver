# Screened Monolayer eps_in Sweep Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Find an in-slab permittivity `eps_in ∈ [2.0, 3.5]` at fixed `LZ = 3.35 Å` that drives the screened density channels (onsite, nn) within 5% of CoQui's `U_ijkl` at `k_323201`, by running the bilayer dielectric solver at `eps_in ∈ {2.0, 2.4, 2.7, 3.0, 3.2, 3.5}`.

**Architecture:** Mirror of the just-completed LZ-sweep plan. Black-box reuse of `monolayer/src/MonolayerScreenedSolve.jl` + `bilayer_slab/src/ScreenedOrbitalSolve.jl`. The source `VolumeSource`s are built once outside the eps_in loop because geometry is fixed across rows; the screened-source rescaling and the interface response are rebuilt per eps_in inside `solve_screened_mode`. The `eps_in = 2.4` row provides a baseline cross-check against the 2026-05-13 k_323201 numbers.

**Tech Stack:** Julia 1.11 via juliaup, shared `Project.toml`. Rusty `ccm` / genoa, 48 cores, ~8 min wall.

---

## File Structure

Paths relative to `/mnt/home/xgao1/work/four_index_integral_solver/codes/graphene/`.

**New files:**

- `monolayer/scripts/screened_monolayer_epsin_sweep.jl` — Julia driver (~80 lines).
- `monolayer/scripts/submit_screened_monolayer_epsin_sweep.slurm` — Slurm batch script for Rusty `ccm`/genoa, 48 cores, 1 h wall.
- `monolayer/data/screened_monolayer_epsin_sweep.csv` — 12 rows (6 eps_in × 2 channels). Produced by the Slurm job.
- `monolayer/results/2026-05-14-screened-monolayer-epsin-sweep-report.md` — dated report.

**Reused unchanged (black box):**

- `monolayer/src/MonolayerScreenedSolve.jl` — `parse_coqui_loc`, `monolayer_screened_sources`, `density_target_specs`, `screened_run_record`.
- `bilayer_slab/src/ScreenedOrbitalSolve.jl` — `solve_screened_mode`.
- `monolayer/src/MonolayerOrbitalLoader.jl` — orbital loaders.

No test additions needed — the helpers are already covered by `monolayer/test/test_monolayer_screened_solve.jl` (45 assertions).

---

## Task 1: Confirm prerequisites

**Files:** none (read-only)

- [ ] **Step 1: Verify the existing baseline artifacts are in place**

```bash
cd /mnt/home/xgao1/work/four_index_integral_solver/codes/graphene
ls -1 monolayer/scripts/screened_monolayer_k323201.jl \
      monolayer/scripts/screened_monolayer_lz_sweep.jl \
      monolayer/src/MonolayerScreenedSolve.jl \
      monolayer/data/screened_monolayer_k323201.csv \
      monolayer/data/screened_monolayer_lz_sweep.csv
```

Expected: all 5 files listed. If any missing, stop and report BLOCKED.

- [ ] **Step 2: Capture the k_323201 baseline values for later cross-check**

```bash
julia --project=. -e '
using CSV, DataFrames
df = CSV.read("monolayer/data/screened_monolayer_k323201.csv", DataFrame)
sub = df[(df.mode .== "Sharp") .& (df.channel .∈ Ref(["onsite", "nn"])), [:channel, :u_total_ev]]
println(sub)
'
```

Expected output (these are the targets for the eps_in=2.4 cross-check in Task 4):

```
2×2 DataFrame
 Row │ channel  u_total_ev
─────┼─────────────────────
   1 │ onsite     11.7471
   2 │ nn          6.113
```

Record these two values (onsite ≈ 11.7471, nn ≈ 6.113) for the cross-check in Task 5.

---

## Task 2: Write the eps_in-sweep driver

**Files:**
- Create: `monolayer/scripts/screened_monolayer_epsin_sweep.jl`

The driver mirrors `screened_monolayer_lz_sweep.jl` but with a `EPS_IN_VALUES` loop in `main()`, `LZ`/`Z_CENTER`/`EPS_OUT` held as constants, and the source built once outside the loop (because geometry is now fixed across rows).

- [ ] **Step 1: Create the driver**

Create `monolayer/scripts/screened_monolayer_epsin_sweep.jl`:

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

# Geometry (fixed)
const LZ       = 3.35
const Z_CENTER = LZ / 2
const EPS_OUT  = 1.0
const L        = 90.0

# eps_in values to sweep. Includes 2.4 baseline for cross-check.
const EPS_IN_VALUES = [2.0, 2.4, 2.7, 3.0, 3.2, 3.5]

# Solver tolerances (reused from k_323201 / LZ-sweep protocols)
const SOURCE_TOL = 1e-3
const N_QUAD = 6
const EDGE_REFINE_LEVEL = 4
const RHS_TOL = 1e-3
const LHS_TOL = 1e-5
const GMRES_ATOL = 1e-5
const GMRES_RTOL = 1e-5
const MAX_ORDER = 64
const MAX_DEPTH = 12

const OUT_CSV = joinpath(@__DIR__, "..", "data", "screened_monolayer_epsin_sweep.csv")

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

function _record_kwargs(eps_in)
    return (
        mode_label = "Sharp", bandwidth = missing,
        eps_in = eps_in, eps_out = EPS_OUT, Lz = LZ, L = L, z_center = Z_CENTER,
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

    println("Building sources (geometry is fixed across the sweep) ...")
    src = monolayer_screened_sources(
        orbital_1 = XSF_1, orbital_2 = XSF_2,
        source_tol = SOURCE_TOL, z_center = Z_CENTER,
    )
    println("  Nphi1 = ", src.Nphi1, "  Nphi2 = ", src.Nphi2)

    density_specs = density_target_specs(src)

    rows = NamedTuple[]

    for eps_in in EPS_IN_VALUES
        println("--- eps_in = ", eps_in, " (LZ = ", LZ, ", Z_CENTER = ", Z_CENTER, ") ---")
        result = solve_screened_mode(
            src.vs1, density_specs,
            L, L, LZ, eps_in, EPS_OUT, BI.SharpScreening();
            _solve_kwargs()...,
        )
        println("  residual = ", result.residual,
                "  n_interface_points = ", result.n_interface_points)

        for pr in result.pair_results
            ch_sym = pr.pair == "onsite" ? :onsite : :nn
            row0 = screened_run_record(pr, result; _record_kwargs(eps_in)...)
            U_ref = coqui[ch_sym].U_ijkl
            row = merge(row0, (
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
    summary_cols = [:eps_in, :channel, :u_int_ev, :u_scatter_ev, :u_total_ev, :coqui_U_ijkl, :rel_err_pct, :sigma_residual]
    show(table[:, summary_cols]; allrows = true, allcols = true)
    println()
end

main()
```

- [ ] **Step 2: Confirm the file parses cleanly (do NOT run main)**

```bash
julia --project=. -e '
contents = read("monolayer/scripts/screened_monolayer_epsin_sweep.jl", String)
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

Use bash timeout 60000 ms. Expected: `parse_ok = true`.

- [ ] **Step 3: Commit the driver (no run yet)**

```bash
git add monolayer/scripts/screened_monolayer_epsin_sweep.jl
git commit -m "$(cat <<'EOF'
feat(monolayer): eps_in-sweep driver at fixed LZ=3.35 (Sharp density only)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Write the Slurm submit script

**Files:**
- Create: `monolayer/scripts/submit_screened_monolayer_epsin_sweep.slurm`

- [ ] **Step 1: Create the submit script**

Create `monolayer/scripts/submit_screened_monolayer_epsin_sweep.slurm`:

```bash
#!/bin/bash
#
# Slurm submission script for the screened-monolayer eps_in sweep at fixed LZ=3.35.
# Submit from the graphene project root (path containing Project.toml):
#
#     cd /mnt/home/xgao1/work/four_index_integral_solver/codes/graphene
#     sbatch monolayer/scripts/submit_screened_monolayer_epsin_sweep.slurm
#
# Wall-clock estimate: ~8 min for 6 Sharp density solves (the LZ sweep did 5
# solves in 6m40s; geometry is fixed across rows so source-building only
# happens once). 1-hour limit gives ample margin.

#SBATCH --job-name=epsin_sweep_k323201
#SBATCH --output=monolayer/logs/epsin_sweep_%j.out
#SBATCH --error=monolayer/logs/epsin_sweep_%j.err
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=01:00:00

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
    julia --project=. monolayer/scripts/screened_monolayer_epsin_sweep.jl

echo ""
echo "completed at: $(date -Is)"
```

- [ ] **Step 2: Commit**

```bash
git add monolayer/scripts/submit_screened_monolayer_epsin_sweep.slurm
git commit -m "$(cat <<'EOF'
feat(monolayer): Slurm submit script for eps_in sweep

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Submit and wait for the cluster job

Per organization cluster-policy LAW, the controller and any agent may not execute `sbatch`. The user submits.

- [ ] **Step 1: User submits**

Ask the user to run, in this session via the `!` prefix or in their own terminal:

```
! sbatch monolayer/scripts/submit_screened_monolayer_epsin_sweep.slurm
```

Capture the job ID from `Submitted batch job XXXXXXX`.

- [ ] **Step 2: Monitor**

```bash
JOBID=<paste job id here>
squeue -j "$JOBID"
sacct -j "$JOBID" --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

Wait for `State = COMPLETED` and `ExitCode = 0:0`. Do not poll repeatedly — check once after ~10 min, or wait for the user to confirm completion.

- [ ] **Step 3: Verify outputs**

```bash
ls -la monolayer/logs/epsin_sweep_${JOBID}.{out,err}
ls -la monolayer/data/screened_monolayer_epsin_sweep.csv
tail -30 monolayer/logs/epsin_sweep_${JOBID}.out
```

Expected:
- `.out` ends with `Wrote 12 rows to .../screened_monolayer_epsin_sweep.csv` and a 12-row DataFrame summary.
- CSV exists, 13 lines (header + 12 data rows).

Quick sanity + baseline cross-check via Julia:

```bash
julia --project=. -e '
using CSV, DataFrames
df = CSV.read("monolayer/data/screened_monolayer_epsin_sweep.csv", DataFrame)
@assert nrow(df) == 12 "expected 12 rows, got $(nrow(df))"
@assert all(df.sigma_residual .< 1e-3) "some sigma_residual not converged: max = $(maximum(df.sigma_residual))"
@assert sort!(unique(df.eps_in)) == [2.0, 2.4, 2.7, 3.0, 3.2, 3.5] "eps_in values mismatch"
# Cross-check eps_in=2.4 row vs k_323201 baseline.
sub = df[(df.eps_in .== 2.4), [:channel, :u_total_ev]]
println("eps_in=2.4 cross-check:")
println(sub)
on = first(sub[sub.channel .== "onsite", :u_total_ev])
nn = first(sub[sub.channel .== "nn",     :u_total_ev])
@assert abs(on - 11.7471) < 1e-3 "onsite at eps_in=2.4 ($on) differs from k_323201 baseline (11.7471) by more than 1 mev"
@assert abs(nn -  6.113)  < 1e-3 "nn at eps_in=2.4 ($nn) differs from k_323201 baseline (6.113) by more than 1 mev"
println("OK: ", nrow(df), " rows, max residual = ", maximum(df.sigma_residual),
        ", eps_in=2.4 baseline cross-check passed.")
show(df[:, [:eps_in, :channel, :u_total_ev, :coqui_U_ijkl, :rel_err_pct, :sigma_residual]]; allrows=true, allcols=true)
println()
'
```

If the eps_in=2.4 cross-check fails, the pipeline has regressed since the 2026-05-13 run. STOP and investigate. Otherwise proceed.

- [ ] **Step 4: Commit the CSV**

```bash
git add monolayer/data/screened_monolayer_epsin_sweep.csv
git commit -m "$(cat <<'EOF'
data(monolayer): eps_in-sweep CSV — 6 eps_in × 2 channels at LZ=3.35

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Write the dated report

**Files:**
- Create: `monolayer/results/2026-05-14-screened-monolayer-epsin-sweep-report.md`

- [ ] **Step 1: Extract numbers from the CSV**

```bash
julia --project=. -e '
using CSV, DataFrames, Printf
df = CSV.read("monolayer/data/screened_monolayer_epsin_sweep.csv", DataFrame)
println("eps_in vs rel_err table:")
println(rpad("eps_in", 8), rpad("onsite ours", 14), rpad("onsite rel_err%", 18),
        rpad("nn ours", 12), rpad("nn rel_err%", 14))
for eps_in in sort!(unique(df.eps_in))
    on = df[(df.eps_in .== eps_in) .& (df.channel .== "onsite"), :][1, :]
    nn = df[(df.eps_in .== eps_in) .& (df.channel .== "nn"),     :][1, :]
    @printf("%-8.2f%-14.4f%-18.3f%-12.4f%-14.3f\n",
            eps_in, on.u_total_ev, on.rel_err_pct, nn.u_total_ev, nn.rel_err_pct)
end
println()
println("Best eps_in by max(|rel_err_onsite|, |rel_err_nn|):")
best_eps = nothing
best_max = Inf
for eps_in in sort!(unique(df.eps_in))
    on = df[(df.eps_in .== eps_in) .& (df.channel .== "onsite"), :][1, :]
    nn = df[(df.eps_in .== eps_in) .& (df.channel .== "nn"),     :][1, :]
    m = max(abs(on.rel_err_pct), abs(nn.rel_err_pct))
    println("  eps_in=", eps_in, "  max_abs_rel_err=", m, "%")
    if m < best_max
        best_max = m
        best_eps = eps_in
    end
end
println("  best eps_in = ", best_eps, " with max_abs_rel_err = ", best_max, "%")
'
```

Record the best eps_in and its max-abs-rel-err. Determine verdict per spec:
- Strong: best max ≤ 5%
- Acceptable: best max ≤ 15%
- Disagreement: otherwise

- [ ] **Step 2: Create the report**

Create `monolayer/results/2026-05-14-screened-monolayer-epsin-sweep-report.md`. Replace every `<fill>` with the actual numbers; choose the verdict per criteria.

```markdown
# Screened Monolayer Graphene — eps_in Sweep at LZ=3.35 — Report

Date: 2026-05-14

## Context

The 2026-05-13 LZ sweep at fixed `eps_in = 2.4` showed monotonic improvement in CoQui-vs-ours agreement as `LZ` grew, with the best swept value LZ=8 giving onsite +15.7% / nn +10.2% — Acceptable for nn but on the boundary of Disagreement for onsite. The linear extrapolation suggested LZ alone could not reach the Strong (5%) band within physically reasonable slab thicknesses, motivating an eps_in sweep at fixed LZ.

This report runs that sweep, holding `LZ = 3.35 Å` (graphite interlayer spacing, the original 2026-05-13 baseline) and sweeping `eps_in ∈ {2.0, 2.4, 2.7, 3.0, 3.2, 3.5}`.

Reference: `/mnt/ceph/users/mroesner/Graphene/cRPA4RSGW/graphene/monolayer/k_323201_nb_144_c_15/_coqui_crpa_loc.out`.

| channel | U_ijkl (eV) |
|---------|------------:|
| onsite  | 9.7824      |
| nn      | 5.1678      |

## Method

| parameter | value |
|-----------|------:|
| LZ        | 3.35 Å (fixed) |
| Z_CENTER  | 1.675 Å (LZ/2, symmetric) |
| eps_out   | 1.0 (fixed) |
| eps_in values | 2.0, 2.4, 2.7, 3.0, 3.2, 3.5 |
| L (in-plane) | 90 Å |
| mode | Sharp only |
| solver tolerances | identical to the 2026-05-13 k_323201 / LZ-sweep runs |

Sources are built once outside the eps_in loop (geometry is fixed across rows). The screened-source rescaling and the BEM operator are rebuilt per eps_in inside `solve_screened_mode`.

Compute: Rusty `ccm`/genoa, 48 cores, Slurm job `<fill: jobid>`, `<fill: elapsed>` wall, `<fill: maxRSS>` peak RSS.

## Cross-check: eps_in=2.4 baseline reproduction

The `eps_in=2.4` row of this sweep should reproduce the 2026-05-13 k_323201 Sharp-mode numbers (onsite 11.7471, nn 6.113) to 4 decimals. Observed:

| channel | this sweep (eps_in=2.4) | k_323201 baseline | drift |
|---------|------------------------:|------------------:|------:|
| onsite  | <fill>                  | 11.7471           | <fill> |
| nn      | <fill>                  |  6.113            | <fill> |

<fill: confirm drift < 1e-3, or note discrepancy.>

## Results

| eps_in | onsite ours | onsite rel_err% | nn ours | nn rel_err% | u_scatter onsite | u_scatter nn |
|-------:|------------:|----------------:|--------:|------------:|-----------------:|-------------:|
| 2.0    | <fill>      | <fill>          | <fill>  | <fill>      | <fill>           | <fill>       |
| 2.4    | <fill>      | <fill>          | <fill>  | <fill>      | <fill>           | <fill>       |
| 2.7    | <fill>      | <fill>          | <fill>  | <fill>      | <fill>           | <fill>       |
| 3.0    | <fill>      | <fill>          | <fill>  | <fill>      | <fill>           | <fill>       |
| 3.2    | <fill>      | <fill>          | <fill>  | <fill>      | <fill>           | <fill>       |
| 3.5    | <fill>      | <fill>          | <fill>  | <fill>      | <fill>           | <fill>       |

`u_total = u_int + u_scatter`. `rel_err = 100 * (ours − CoQui) / CoQui`. All GMRES residuals: <fill: min> to <fill: max>.

## Verdict

Best `eps_in` by `max(|rel_err_onsite|, |rel_err_nn|)`:

| eps_in | max abs rel_err |
|-------:|----------------:|
| 2.0    | <fill>          |
| 2.4    | <fill>          |
| 2.7    | <fill>          |
| 3.0    | <fill>          |
| 3.2    | <fill>          |
| 3.5    | <fill>          |

**Best eps_in in the swept range: <fill> with max-abs-rel-err <fill>%.**

Spec thresholds:
- Strong: best max ≤ 5%
- Acceptable: best max ≤ 15%
- Disagreement: otherwise

**Verdict: <fill>** — one sentence summarizing whether a single eps_in works for both channels or whether onsite and nn prefer materially different values.

## Interpretation

<fill: a few sentences describing the trends. Specifically address:
  1. Does rel_err monotonically decrease with eps_in for both channels? Or is there an interior minimum?
  2. Do onsite and nn agree on the best eps_in, or do they prefer different values (separation > 0.3)?
  3. Are u_int and u_scatter both moving as expected (u_int decreases with eps_in; u_scatter grows in magnitude)?
  4. If the best eps_in is at the edge of [2.0, 3.5], should the sweep be widened in a follow-up?
  5. Does the structurally-best eps_in approach the "experimental" graphene in-plane permittivity (~2.4-3.0 from literature) or something else?
>

## Reproducibility

- Driver:        `monolayer/scripts/screened_monolayer_epsin_sweep.jl`
- Submit script: `monolayer/scripts/submit_screened_monolayer_epsin_sweep.slurm`
- CSV:           `monolayer/data/screened_monolayer_epsin_sweep.csv`
- Run logs:      `monolayer/logs/epsin_sweep_<fill: jobid>.{out,err}`

## Next

<fill: depending on the verdict:
  - If Strong: declare eps_in calibrated at the best value; replicate at k_252501 and k_161601.
  - If Acceptable: open joint (eps_in, LZ) calibration or widen the eps_in range.
  - If Disagreement and best eps_in at boundary: widen the sweep (e.g., extend to 4.0 or 5.0).
  - If Disagreement and minimum is interior, with onsite/nn preferring different eps_in:
    slab model is structurally limited — recommend eps(z) profile or 2D-Lindhard correction.
>

The bilayer-monolayer follow-up Malte mentioned remains blocked on him producing a starting example.
```

- [ ] **Step 3: Commit**

```bash
git add monolayer/results/2026-05-14-screened-monolayer-epsin-sweep-report.md
git commit -m "$(cat <<'EOF'
docs(monolayer): eps_in-sweep results + verdict

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

- `docs(monolayer): eps_in-sweep results + verdict`
- `data(monolayer): eps_in-sweep CSV — 6 eps_in × 2 channels at LZ=3.35`
- `feat(monolayer): Slurm submit script for eps_in sweep`
- `feat(monolayer): eps_in-sweep driver at fixed LZ=3.35 (Sharp density only)`
- `spec: eps_in sweep ...`

- [ ] **Step 2: Test suite still green**

```bash
julia --project=. test/runtests.jl
```

Bash timeout 600000 ms. Expected: all 189 tests pass. This plan adds no tests so the number is unchanged.

- [ ] **Step 3: CSV sanity**

```bash
julia --project=. -e '
using CSV, DataFrames
df = CSV.read("monolayer/data/screened_monolayer_epsin_sweep.csv", DataFrame)
@assert nrow(df) == 12
@assert all(df.sigma_residual .< 1e-3)
@assert sort!(unique(df.eps_in)) == [2.0, 2.4, 2.7, 3.0, 3.2, 3.5]
println("OK: ", nrow(df), " rows, all sigma_residual < 1e-3")
'
```

---

## Self-Review

**Spec coverage:**
- Spec § "Scope → In scope" — 6 eps_in × 1 mode × 2 channels, CoQui parse once, baseline cross-check: Task 2 driver covers this with `EPS_IN_VALUES` const and the per-eps_in loop. Task 4 Step 3 codifies the baseline cross-check assertion.
- Spec § "File layout" — Task 2 (driver), Task 3 (submit), Task 4 (CSV produced + committed), Task 5 (report).
- Spec § "Slab parameters" — `const` declarations in the driver header.
- Spec § "CSV schema" — Task 2's `push!(rows, merge(row0, (coqui_v_ijkl=…, coqui_U_ijkl=…, diff_ev=…, rel_err_pct=…)))` adds the comparison columns on top of the 28-field `screened_run_record`. `eps_in` is sourced from `_record_kwargs(eps_in)` and ends up in the row via `screened_run_record` directly.
- Spec § "Verdict criteria" — Task 5 Step 1 (max-of-abs computation) + Task 5 Step 2 (verdict line in the report template).
- Spec § "Risks / loose ends" — Risk #1 (GMRES at large eps_in) handled by Task 4's residual sanity check; Risk #2 (source identity) handled by building sources once outside the loop in the driver; Risk #3 (baseline cross-check) is an explicit assertion in Task 4 Step 3.

**Placeholder scan:**
- `<fill>` markers appear only in the report template (Task 5 Step 2). Task 5 Step 1's extraction script feeds the values that populate them.
- No `TBD`/`TODO`/"appropriate"/"handle edge cases" anywhere.

**Type consistency:**
- `parse_coqui_loc` → `Dict{Symbol, NamedTuple{(:v_ijkl, :U_ijkl), Tuple{Float64, Float64}}}` (same as the k_323201 / LZ-sweep plans); consumed via `coqui[:onsite]` / `coqui[:nn]`. ✓
- `monolayer_screened_sources(; orbital_1, orbital_2, source_tol, z_center)` keyword names match. ✓
- `screened_run_record(pair_result, solve_result; mode_label, bandwidth, eps_in, eps_out, Lz, L, z_center, source_tol, rhs_tol, lhs_tol, gmres_atol, gmres_rtol, n_quad, edge_refine_level, max_order, max_depth)` — `_record_kwargs(eps_in)` produces exactly these 17 keys. ✓
- The CSV column merge in Task 2 adds `coqui_v_ijkl, coqui_U_ijkl, diff_ev, rel_err_pct` on top of the 28-field record. Note `eps_in` is already in the 28-field record (it's one of the kwargs to `screened_run_record`), so it's NOT in the merge tuple. Task 5 Step 1 reads `eps_in, channel, u_total_ev, rel_err_pct, sigma_residual` — all defined. ✓
- The Slurm submit script mirrors the existing LZ-sweep `submit_screened_monolayer_lz_sweep.slurm` line-for-line apart from job name, output paths, and wall time. ✓

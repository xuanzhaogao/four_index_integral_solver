# Screened Monolayer Hund's at Calibrated (eps_in=3.5, LZ=3.35) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run one `BI.SharpScreening()` Hund's solve at the calibrated `(eps_in = 3.5, LZ = 3.35)` and report the U_hund value alongside (a) CoQui's reference 0.0937 eV and (b) the previous `eps_in = 2.4` Hund's value of 0.0902 eV.

**Architecture:** Black-box reuse of `monolayer/src/MonolayerScreenedSolve.jl` (parser, source builder, `hund_target_specs`, `screened_run_record`) and `bilayer_slab/src/ScreenedOrbitalSolve.jl` (`solve_screened_mode`). A short driver script runs one solve, writes a 1-row CSV. A short Slurm submit script runs it on Rusty `ccm`/genoa.

**Tech Stack:** Julia 1.11 via juliaup, shared `Project.toml`. Rusty `ccm`/genoa, 48 cores, ~2 min wall.

---

## File Structure

Paths relative to `/mnt/home/xgao1/work/four_index_integral_solver/codes/graphene/`.

**New files:**

- `monolayer/scripts/screened_monolayer_hund_calibrated.jl` — Julia driver (~70 lines).
- `monolayer/scripts/submit_screened_monolayer_hund_calibrated.slurm` — Slurm submit, 30 min wall.
- `monolayer/data/screened_monolayer_hund_calibrated.csv` — 1 data row. Produced by the Slurm job.
- `monolayer/results/2026-05-14-screened-monolayer-hund-calibrated-report.md` — brief dated report.

**Reused unchanged (black box):**

- `monolayer/src/MonolayerScreenedSolve.jl` — `parse_coqui_loc`, `monolayer_screened_sources`, `hund_target_specs`, `screened_run_record`.
- `bilayer_slab/src/ScreenedOrbitalSolve.jl` — `solve_screened_mode`.

No test additions.

---

## Task 1: Confirm prerequisites

**Files:** none (read-only)

- [ ] **Step 1: Verify the calibrated-eps_in artifacts and the 2026-05-13 baseline are present**

```bash
cd /mnt/home/xgao1/work/four_index_integral_solver/codes/graphene
ls -1 monolayer/scripts/screened_monolayer_epsin_sweep.jl \
      monolayer/src/MonolayerScreenedSolve.jl \
      monolayer/data/screened_monolayer_k323201.csv \
      monolayer/data/screened_monolayer_epsin_sweep.csv
```

Expected: all 4 files listed. If any missing, BLOCKED.

- [ ] **Step 2: Capture the baseline Hund's value from k_323201**

```bash
julia --project=. -e '
using CSV, DataFrames
df = CSV.read("monolayer/data/screened_monolayer_k323201.csv", DataFrame)
sub = df[(df.channel .== "hund") .& (df.mode .== "Sharp"), [:channel, :eps_in, :u_total_ev, :u_int_ev, :u_scatter_ev]]
println(sub)
'
```

Expected: one row with `u_total_ev ≈ 0.0902` at `eps_in = 2.4`. Record this exact value for the report's baseline comparison.

---

## Task 2: Write the driver script

**Files:**
- Create: `monolayer/scripts/screened_monolayer_hund_calibrated.jl`

- [ ] **Step 1: Create the driver**

Create `monolayer/scripts/screened_monolayer_hund_calibrated.jl`:

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

# Calibrated slab parameters from 2026-05-14 eps_in sweep
const LZ       = 3.35
const Z_CENTER = LZ / 2
const EPS_IN   = 3.5
const EPS_OUT  = 1.0
const L        = 90.0

const SOURCE_TOL = 1e-3
const N_QUAD = 6
const EDGE_REFINE_LEVEL = 4
const RHS_TOL = 1e-3
const LHS_TOL = 1e-5
const GMRES_ATOL = 1e-5
const GMRES_RTOL = 1e-5
const MAX_ORDER = 64
const MAX_DEPTH = 12

const OUT_CSV = joinpath(@__DIR__, "..", "data", "screened_monolayer_hund_calibrated.csv")

function main()
    println("Parsing CoQui reference ...")
    coqui = parse_coqui_loc(COQUI)
    println("  hund_sf U_ijkl = ", coqui[:hund_sf].U_ijkl)

    println("Building sources at z_center = ", Z_CENTER, " ...")
    src = monolayer_screened_sources(
        orbital_1 = XSF_1, orbital_2 = XSF_2,
        source_tol = SOURCE_TOL, z_center = Z_CENTER,
    )
    println("  Nphi1 = ", src.Nphi1, "  Nphi2 = ", src.Nphi2)

    hund_specs = hund_target_specs(src)
    println("--- Hund's solve  eps_in = ", EPS_IN, "  LZ = ", LZ, " ---")
    result = solve_screened_mode(
        src.vs_product, hund_specs,
        L, L, LZ, EPS_IN, EPS_OUT, BI.SharpScreening();
        n_quad = N_QUAD, edge_refine_level = EDGE_REFINE_LEVEL,
        rhs_tol = RHS_TOL, lhs_tol = LHS_TOL,
        gmres_atol = GMRES_ATOL, gmres_rtol = GMRES_RTOL,
        max_order = MAX_ORDER, max_depth = MAX_DEPTH,
        volume_tol = RHS_TOL,
    )
    println("  residual = ", result.residual,
            "  n_interface_points = ", result.n_interface_points)

    rows = NamedTuple[]
    for pr in result.pair_results
        row0 = screened_run_record(pr, result;
            mode_label = "Sharp", bandwidth = missing,
            eps_in = EPS_IN, eps_out = EPS_OUT, Lz = LZ, L = L, z_center = Z_CENTER,
            source_tol = SOURCE_TOL, rhs_tol = RHS_TOL, lhs_tol = LHS_TOL,
            gmres_atol = GMRES_ATOL, gmres_rtol = GMRES_RTOL,
            n_quad = N_QUAD, edge_refine_level = EDGE_REFINE_LEVEL,
            max_order = MAX_ORDER, max_depth = MAX_DEPTH,
        )
        U_ref = coqui[:hund_sf].U_ijkl
        row = merge(row0, (
            coqui_v_ijkl = coqui[:hund_sf].v_ijkl,
            coqui_U_ijkl = U_ref,
            diff_ev = U_ref - row0.u_total_ev,
            rel_err_pct = 100 * (row0.u_total_ev - U_ref) / U_ref,
        ))
        push!(rows, row)
        println("    ", pr.pair, "  ours=", row0.u_total_ev,
                "  CoQui=", U_ref,
                "  rel_err_pct=", row.rel_err_pct)
    end

    mkpath(dirname(OUT_CSV))
    table = DataFrame(rows)
    CSV.write(OUT_CSV, table)
    println("Wrote $(nrow(table)) rows to $(OUT_CSV)")
    summary_cols = [:channel, :eps_in, :Lz, :u_int_ev, :u_scatter_ev, :u_total_ev, :coqui_U_ijkl, :rel_err_pct, :sigma_residual]
    show(table[:, summary_cols]; allrows = true, allcols = true)
    println()
end

main()
```

- [ ] **Step 2: Parse-check**

```bash
julia --project=. -e '
contents = read("monolayer/scripts/screened_monolayer_hund_calibrated.jl", String)
import Base.Meta
parse_ok = try
    Meta.parseall(contents); true
catch err
    println("PARSE ERROR: ", err); false
end
println("parse_ok = ", parse_ok)
'
```

Bash timeout 60000 ms. Expected: `parse_ok = true`.

- [ ] **Step 3: Commit**

```bash
git add monolayer/scripts/screened_monolayer_hund_calibrated.jl
git commit -m "$(cat <<'EOF'
feat(monolayer): Hund's-channel driver at calibrated (eps_in=3.5, LZ=3.35)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Write the Slurm submit script

**Files:**
- Create: `monolayer/scripts/submit_screened_monolayer_hund_calibrated.slurm`

- [ ] **Step 1: Create the submit script**

Create `monolayer/scripts/submit_screened_monolayer_hund_calibrated.slurm`:

```bash
#!/bin/bash
#
# Slurm submission script for the calibrated-Hund's evaluation.
# Submit from the graphene project root:
#
#     cd /mnt/home/xgao1/work/four_index_integral_solver/codes/graphene
#     sbatch monolayer/scripts/submit_screened_monolayer_hund_calibrated.slurm
#
# Wall-clock estimate: ~2 min for one Sharp Hund's solve. 30-min limit is
# generous overhead.

#SBATCH --job-name=hund_calibrated
#SBATCH --output=monolayer/logs/hund_calibrated_%j.out
#SBATCH --error=monolayer/logs/hund_calibrated_%j.err
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=00:30:00

set -euo pipefail

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}

mkdir -p monolayer/logs

echo "host: $(hostname)"
echo "date: $(date -Is)"
echo "julia: $(which julia)"
julia --version
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS JULIA_NUM_THREADS=$JULIA_NUM_THREADS"
echo ""

srun --cpus-per-task=$SLURM_CPUS_PER_TASK \
    julia --project=. monolayer/scripts/screened_monolayer_hund_calibrated.jl

echo ""
echo "completed at: $(date -Is)"
```

- [ ] **Step 2: Commit**

```bash
git add monolayer/scripts/submit_screened_monolayer_hund_calibrated.slurm
git commit -m "$(cat <<'EOF'
feat(monolayer): Slurm submit script for calibrated Hund's evaluation

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: User submits, controller monitors

Per organization cluster-policy LAW, the controller cannot execute `sbatch`. The user submits.

- [ ] **Step 1: User submits**

Ask the user to run:

```
! sbatch monolayer/scripts/submit_screened_monolayer_hund_calibrated.slurm
```

Capture the job ID.

- [ ] **Step 2: Wait, then check status**

After waiting ~2-5 min (or after the user confirms completion):

```bash
JOBID=<paste here>
sacct -j "$JOBID" --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

Expected: `State = COMPLETED`, `ExitCode = 0:0`.

- [ ] **Step 3: Verify outputs**

```bash
ls -la monolayer/logs/hund_calibrated_${JOBID}.{out,err}
ls -la monolayer/data/screened_monolayer_hund_calibrated.csv
tail -20 monolayer/logs/hund_calibrated_${JOBID}.out
```

Expected: stdout ends with `Wrote 1 rows to ...`, CSV is 2 lines (header + 1 data row).

Quick sanity check:

```bash
julia --project=. -e '
using CSV, DataFrames
df = CSV.read("monolayer/data/screened_monolayer_hund_calibrated.csv", DataFrame)
@assert nrow(df) == 1 "expected 1 row, got $(nrow(df))"
@assert df.channel[1] == "hund"
@assert df.eps_in[1] == 3.5
@assert df.Lz[1] == 3.35
@assert df.sigma_residual[1] < 1e-3 "GMRES did not converge: residual = $(df.sigma_residual[1])"
println("OK: 1 row, residual = ", df.sigma_residual[1],
        ", u_total_ev = ", df.u_total_ev[1],
        ", CoQui = ", df.coqui_U_ijkl[1],
        ", rel_err_pct = ", df.rel_err_pct[1], "%")
'
```

If `sigma_residual > 1e-3`, GMRES failed to converge — note in the report, but the value is the value (the spec acknowledged this risk).

- [ ] **Step 4: Commit the CSV**

```bash
git add monolayer/data/screened_monolayer_hund_calibrated.csv
git commit -m "$(cat <<'EOF'
data(monolayer): Hund's CSV at calibrated (eps_in=3.5, LZ=3.35) — 1 row

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Write the brief report

**Files:**
- Create: `monolayer/results/2026-05-14-screened-monolayer-hund-calibrated-report.md`

- [ ] **Step 1: Extract values from the new CSV + the k_323201 baseline**

```bash
julia --project=. -e '
using CSV, DataFrames
new_df = CSV.read("monolayer/data/screened_monolayer_hund_calibrated.csv", DataFrame)
new_row = new_df[1, :]
println("Calibrated (eps_in=3.5):")
println("  u_int_ev       = ", new_row.u_int_ev)
println("  u_scatter_ev   = ", new_row.u_scatter_ev)
println("  u_total_ev     = ", new_row.u_total_ev)
println("  coqui_U_ijkl   = ", new_row.coqui_U_ijkl)
println("  rel_err_pct    = ", new_row.rel_err_pct, "%")
println("  sigma_residual = ", new_row.sigma_residual)
println("  n_interface_points = ", new_row.n_interface_points)
println()
base_df = CSV.read("monolayer/data/screened_monolayer_k323201.csv", DataFrame)
base_hund = base_df[(base_df.channel .== "hund") .& (base_df.mode .== "Sharp"), :][1, :]
println("Baseline (eps_in=2.4, from k_323201 CSV):")
println("  u_int_ev       = ", base_hund.u_int_ev)
println("  u_scatter_ev   = ", base_hund.u_scatter_ev)
println("  u_total_ev     = ", base_hund.u_total_ev)
'
```

Record these values.

- [ ] **Step 2: Create the report**

Create `monolayer/results/2026-05-14-screened-monolayer-hund-calibrated-report.md`. Replace `<fill>` with the actual values from Step 1; choose the verdict per Strong/Acceptable/Disagreement.

```markdown
# Screened Monolayer Hund's at Calibrated (eps_in=3.5, LZ=3.35) — Report

Date: 2026-05-14

## Context

The 2026-05-14 eps_in sweep calibrated `(eps_in = 3.5, LZ = 3.35 Å)` as the slab-dielectric parameters that drive the screened density channels (onsite, nn) to within ~1.2 % of CoQui at `k_323201`. The Hund's channel was deferred in that sweep — it was already in Strong band at the original `eps_in = 2.4` (ours 0.0902 vs CoQui 0.0937, −3.7 %). This report measures U_hund at the new calibrated parameters.

Reference: `_coqui_crpa_loc.out` at k_323201.

| channel  | U_ijkl (eV) |
|----------|------------:|
| hund_sf  | 0.0937      |
| hund_ph  | 0.0937      |

## Method

| parameter | value |
|-----------|------:|
| LZ        | 3.35 Å (fixed) |
| Z_CENTER  | 1.675 Å (LZ/2, symmetric) |
| eps_in    | 3.5 (calibrated) |
| eps_out   | 1.0 |
| L (in-plane) | 90 Å |
| mode | Sharp only |
| source | `phi_1 · phi_2` (signed product VolumeSource) |
| solver tolerances | identical to k_323201 / eps_in-sweep runs |

Compute: Rusty `ccm`/genoa, 48 cores, Slurm job `<fill: jobid>`, `<fill: elapsed>` wall, `<fill: maxRSS>` peak RSS. GMRES residual `<fill>` (target < 10⁻³).

## Result

Comparison of Hund's U across the two eps_in values:

| eps_in | ours `U_hund` (eV) | u_int (eV) | u_scatter (eV) | CoQui U (eV) | rel_err % |
|-------:|-------------------:|-----------:|---------------:|-------------:|----------:|
| 2.4    |             0.0902 | 0.0930     | −0.00276       | 0.0937       | −3.7      |
| 3.5    | <fill>             | <fill>     | <fill>         | 0.0937       | <fill>    |

`u_total = u_int + u_scatter`. `rel_err = 100 * (ours − CoQui) / CoQui`.

## Verdict

Spec thresholds:
- Strong: `|rel_err| ≤ 5 %`
- Acceptable: `|rel_err| ≤ 15 %`
- Disagreement: otherwise

**Verdict: <fill>**

## Interpretation

<fill: a few sentences. Specifically address:
  - Did increasing eps_in from 2.4 to 3.5 shift U_hund by a small (~few %) or large (~tens of %) amount?
  - The expectation was small — Hund's source is a near-zero-charge dipole-like distribution that couples weakly to the long-range monopole part of the dielectric response, so the bulk of the additional eps_in screening should not reach Hund's. Did the result confirm this?
  - Compare u_scatter at eps_in=3.5 vs eps_in=2.4. If u_scatter grew by less than the density channels' u_scatter scaled, that quantitatively supports the dipole-weak-coupling story.
  - Net: does the calibrated parameter set (eps_in=3.5, LZ=3.35) simultaneously give Strong agreement on ALL three channels — onsite, nn, and Hund's?
>

## Reproducibility

- Driver:        `monolayer/scripts/screened_monolayer_hund_calibrated.jl`
- Submit script: `monolayer/scripts/submit_screened_monolayer_hund_calibrated.slurm`
- CSV:           `monolayer/data/screened_monolayer_hund_calibrated.csv`
- Run logs:      `monolayer/logs/hund_calibrated_<fill: jobid>.{out,err}`

## Next

<fill: depending on the verdict:
  - If Strong: declare calibration complete on all 4 channels; ready to proceed to k-mesh replication or write up the full slab-model story for Malte.
  - If Acceptable: note the slight worsening from eps_in=2.4; check whether onsite/nn agreement at eps_in=3.5 is worth the small Hund's degradation. If yes, hold the calibration. If no, consider an eps_in slightly below 3.5 (e.g., 3.4) as a compromise.
  - If Disagreement: investigate. The Hund's source's dipole-like character shouldn't be eps_in-sensitive at this magnitude; an unexpected shift would indicate either GMRES instability at higher eps_in or a model issue.
>

The bilayer-monolayer follow-up Malte mentioned remains blocked on him producing a starting example.
```

- [ ] **Step 3: Commit**

```bash
git add monolayer/results/2026-05-14-screened-monolayer-hund-calibrated-report.md
git commit -m "$(cat <<'EOF'
docs(monolayer): Hund's report at calibrated (eps_in=3.5, LZ=3.35)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Final verification

- [ ] **Step 1: Working tree clean, expected commits in place**

```bash
git status
git log --oneline -6
```

Expected commits (newest first):

- `docs(monolayer): Hund's report at calibrated (eps_in=3.5, LZ=3.35)`
- `data(monolayer): Hund's CSV at calibrated (eps_in=3.5, LZ=3.35) — 1 row`
- `feat(monolayer): Slurm submit script for calibrated Hund's evaluation`
- `feat(monolayer): Hund's-channel driver at calibrated (eps_in=3.5, LZ=3.35)`
- `spec: Hund's-channel evaluation at calibrated (eps_in=3.5, LZ=3.35)`

- [ ] **Step 2: Test suite still green**

```bash
julia --project=. test/runtests.jl
```

Bash timeout 600000 ms. Expected: all 189 tests still pass (no test additions).

- [ ] **Step 3: CSV sanity**

```bash
julia --project=. -e '
using CSV, DataFrames
df = CSV.read("monolayer/data/screened_monolayer_hund_calibrated.csv", DataFrame)
@assert nrow(df) == 1
@assert df.channel[1] == "hund"
@assert df.eps_in[1] == 3.5
@assert df.Lz[1] == 3.35
@assert df.sigma_residual[1] < 1e-3
println("OK")
'
```

---

## Self-Review

**Spec coverage:**
- Spec § "Scope → In scope" — 1 solve, 1 row CSV, brief report: Task 2 (driver), Task 4 (run + CSV), Task 5 (report).
- Spec § "File layout" — all 4 files created or produced by an explicit task.
- Spec § "Slab parameters" — encoded as `const` declarations in the driver.
- Spec § "CSV schema" — single-row CSV with the same columns as the eps_in-sweep CSV (28 standard + 4 comparison fields).
- Spec § "Verdict criteria" — encoded in the report template (Task 5 Step 2).
- Spec § "Risks" — Risk #1 (GMRES at higher eps_in) handled by the residual sanity check in Task 4 Step 3.

**Placeholder scan:**
- `<fill>` markers appear only in the report template (Task 5 Step 2). Task 5 Step 1 extracts the values that populate them.
- No `TBD`/`TODO` in script or Slurm.

**Type consistency:**
- `parse_coqui_loc` → `Dict{Symbol, NamedTuple{(:v_ijkl, :U_ijkl), Tuple{Float64, Float64}}}` — consumed via `coqui[:hund_sf]` (same as previous plans). ✓
- `monolayer_screened_sources(; orbital_1, orbital_2, source_tol, z_center)` — keyword names match. ✓
- `hund_target_specs(src)` returns a 1-entry vector with fields `pair, target_vs, Na, Nb`. ✓
- `screened_run_record(pair_result, solve_result; mode_label, bandwidth, eps_in, eps_out, Lz, L, z_center, source_tol, rhs_tol, lhs_tol, gmres_atol, gmres_rtol, n_quad, edge_refine_level, max_order, max_depth)` — all 17 kwargs supplied directly in the driver. ✓
- CSV-extension fields (`coqui_v_ijkl`, `coqui_U_ijkl`, `diff_ev`, `rel_err_pct`) merged onto the record. The CSV sanity check in Task 4 reads `channel, eps_in, Lz, sigma_residual, u_total_ev, coqui_U_ijkl, rel_err_pct` — all present. ✓

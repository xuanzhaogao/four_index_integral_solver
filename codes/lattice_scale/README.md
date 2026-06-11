# lattice_scale — multi-node four-index campaign

Spec: BoundaryIntegral.jl `docs/superpowers/specs/2026-06-10-multinode-lattice-campaign-design.md`.
Plan: BoundaryIntegral.jl `docs/superpowers/plans/2026-06-10-multinode-lattice-campaign.md`.

This env `dev`s `BoundaryIntegral.jl` (branch multi_rhs). Build the env once:
`julia --project -e 'using Pkg; Pkg.instantiate()'`.

## Pipeline

`prepare` (manifest) → `solve` (per-batch σ + truncated ρ, multi-node) → `consolidate`
(shared target set + ρ store) → `eval` (per-batch Φ at the shared targets + all-pair
contraction, multi-node) → `assemble` (dense V + V/Vᵀ symmetry diagnostic).

Status is derived from files on ceph (atomic writes), so any crash/walltime kill is
recovered by **resubmitting** — completed batches are skipped.
`driver.jl <c>.toml status` shows pending counts.

## Run order

**One-time setup:** submit Slurm jobs from the `codes/lattice_scale` directory so that
the relative paths `logs/` and `jobscripts/` resolve correctly (Slurm opens `--output`
before the script body runs).  `logs/` is tracked via `logs/.gitkeep`, so it exists
after a fresh clone; if you clean-clone onto a new machine run
`cd codes/lattice_scale` once before any `sbatch`.

1. `julia --project driver.jl campaigns/<c>.toml prepare`        (login/workstation; writes manifest)
2. Pilot ONE batch on a node — measure before scaling:
   `julia --project driver.jl campaigns/<c>.toml solve --only 1`  (via ssh to an interactive node)
   Note `t_solve`, `dof`, and the batch file size from the stats; size the campaign before step 3.
3. From `codes/lattice_scale`:
   `sbatch --nodes=<M> jobscripts/solve.sbatch campaigns/<c>.toml`   (**you** submit)
4. `julia --project driver.jl campaigns/<c>.toml consolidate`        (single node)
5. Pilot ONE eval — the u_inc/near-correction cost at the full target set is the
   campaign's biggest unknown; measure before committing nodes:
   `julia --project driver.jl campaigns/<c>.toml eval --only 1`
6. From `codes/lattice_scale`:
   `sbatch --nodes=<M> jobscripts/eval.sbatch campaigns/<c>.toml`     (**you** submit)
7. `julia --project driver.jl campaigns/<c>.toml assemble` → `V_full.jls` + `report.txt`

Crash/walltime recovery: just resubmit step 3 or 6 — status is file-derived and
completed batches are skipped.

## Worker topology

One worker per node (`--ntasks-per-node=1`, `--cpus-per-task=96` on `ccm`/`genoa`);
each task uses the whole node's cores via OMP/OpenBLAS threads (pinned in the sbatch
scripts). The FMM saturates ~32–64 cores, so node-sized tasks are the right grain.

## Slurm note

The sbatch scripts are templates — **submit them yourself** (`sbatch ...`). The driver
detects `SLURM_JOB_ID`/`SLURM_NTASKS` and spawns one Julia worker per task via
`SlurmClusterManager.SlurmManager`. See https://wiki.flatironinstitute.org/SCC/Software/Slurm.

## Performance / pilot-watch notes

Measured during the mandatory one-batch pilots (run order steps 2 and 5), these set node counts:

- **Eval scattered potential is K separate FMMs.** `eval_batch`'s `evaluate_batch_potential`
  builds the corrected layer-potential map (FMM + hcubature near-correction) **once** for the
  target set, then applies it per source column — i.e. K separate `nd=1` FMMs over the full
  target set, not one batched `nd=K` FMM. The expensive near-correction *setup* is amortized
  once; the per-column FMM cost is expected (a batched `nd=K` corrected `pottrg` is a future
  optimization, not a regression). The pilot's eval timing reflects this.
- **Interface post-refinement vs the target cloud** (`laplace3d_pottrg_fmm3d_corrected_hcubature`)
  can grow `num_points(interface)` when ~1e7 targets blanket the domain — the biggest eval-cost
  unknown. Watch the `num of sources: N → M` log line in the pilot.
- **`consolidate` loads all BatchResults at once** (~20 GB for ~200 batches; fine on a 1.5 TB
  node). It only reads gidx/weights/densities; σ/interface are loaded but unused.

## Re-preparing invalidates downstream artifacts

`prepare` guards against silently reusing a manifest built with different `nx/ny/cutoff/
n_centers_per_batch` (via `manifest.params`). If you intentionally re-`prepare` a campaign with
changed geometry, also delete the stale `targets.jls`, `rho_store.jls`, `batches/`, `V/`, and
`V_full.jls` under the campaign root — they are tied to the previous pair set.

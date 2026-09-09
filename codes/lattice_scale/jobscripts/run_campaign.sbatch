#!/bin/bash
# Full campaign in ONE allocation — run once and walk away.
#   SUBMIT MANUALLY:  sbatch --nodes=<M> jobscripts/run_campaign.sbatch <campaign.toml>
#
# Runs all five phases in order inside a single Slurm allocation:
#   prepare (head) -> solve (SlurmManager: 1 worker/node) -> consolidate (head)
#   -> eval (SlurmManager) -> assemble (head).
# The two heavy phases (solve, eval) spawn one Julia worker per allocated node;
# the three light phases run on the head. M-1 nodes sit briefly idle during the
# light phases — negligible vs the hours of solve/eval.
#
# Restartable: status is file-derived. If the job is walltime-killed or a phase
# fails, just resubmit the SAME command — finished batches are skipped, and the
# consolidate/assemble guards halt the script if solve/eval are incomplete (so it
# won't run eval on a half-done solve). set -e makes that halt stop the run.
#
# --time must cover solve + eval together (each can be hours at 10x10). Tune it.
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=24:00:00
#SBATCH --job-name=fi_campaign
#SBATCH --output=logs/campaign_%j.out

set -euo pipefail
TOML=${1:?usage: sbatch --nodes=<M> jobscripts/run_campaign.sbatch <campaign.toml>}

# pin BOTH thread pools (unpinned threads silently corrupt timings — known footgun)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# OpenBLAS gets ONE thread, not $SLURM_CPUS_PER_TASK. The hot stages are FMM3D and FINUFFT,
# both with their own pools and almost no BLAS; an equally sized OpenBLAS pool oversubscribes
# the node against them. Measured on the Sec. 5.2 evaluation benchmark: 2.95x slower at 64
# threads and 1.56x at 96 with OPENBLAS_NUM_THREADS=N vs 1.
export OPENBLAS_NUM_THREADS=1
export JULIA_GLUE_THREADS=8

# juliaup launcher (on PATH on compute nodes); the lattice_scale manifest is resolved
# with Julia 1.12 = the juliaup `release` default. If your juliaup default ever moves
# off 1.12, pin with `julia +1.12` instead.
JULIA=/mnt/home/xgao1/.juliaup/bin/julia

cd /mnt/home/xgao1/work/four_index_integral_solver/codes/lattice_scale
mkdir -p logs

run_phase () {   # $1 = phase name, $2 = "-t8" for the parallel phases or "" for head-only
    echo "==================== PHASE: $1  ($(date '+%F %T')) ===================="
    $JULIA --project $2 driver.jl "$TOML" "$1"
    echo "-------------------- $1 done ($(date '+%F %T')) --------------------"
}

run_phase prepare       ""
run_phase solve         "-t $JULIA_GLUE_THREADS"
run_phase consolidate   ""
run_phase eval          "-t $JULIA_GLUE_THREADS"
run_phase assemble      ""

echo "==================== CAMPAIGN COMPLETE ($(date '+%F %T')) ===================="
$JULIA --project driver.jl "$TOML" status

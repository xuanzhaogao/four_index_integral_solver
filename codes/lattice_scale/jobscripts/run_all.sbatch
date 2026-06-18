#!/bin/bash
# FULL four-index campaign in ONE job: prepare -> solve -> consolidate -> eval -> assemble,
# all on the allocated nodes. The solve/eval phases spawn one SlurmManager worker per task
# (--ntasks-per-node=1, one 96-core task per node), so M nodes process the batches M-way
# parallel; prepare/consolidate/assemble run on the head process.
#
# RESUMABLE: if the job hits the walltime or dies, just resubmit the SAME line — completed
# solve/eval batches are skipped (status is derived from files on ceph), and consolidate/
# assemble are idempotent. So a too-short --time is recovered by one more submit.
#
# SUBMIT (you):  sbatch --nodes=10 jobscripts/run_all.sbatch campaigns/lattice_10x10_multicube.toml
#
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=16:00:00
#SBATCH --job-name=fi_all
#SBATCH --output=logs/all_%j.out

set -euo pipefail
TOML=${1:?usage: sbatch --nodes=<M> jobscripts/run_all.sbatch <campaign.toml>}

# pin BOTH thread pools (unpinned threads silently corrupt timings — known footgun)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export JULIA_GLUE_THREADS=8
JULIA=/mnt/home/xgao1/.juliaup/bin/julia      # juliaup launcher (Julia 1.12 = juliaup release default)

cd /mnt/home/xgao1/work/four_index_integral_solver/codes/lattice_scale
mkdir -p logs

for PHASE in prepare solve consolidate eval assemble; do
    echo "=== $PHASE  $(date) ==="
    $JULIA --project -t $JULIA_GLUE_THREADS driver.jl "$TOML" "$PHASE"
done
echo "=== CAMPAIGN DONE  $(date)  ->  $(dirname "$TOML")/../<root>/V_full_eV.jls ==="

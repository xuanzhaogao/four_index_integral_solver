#!/bin/bash
# Eval phase. SUBMIT MANUALLY:  sbatch --nodes=<M> jobscripts/eval.sbatch <campaign.toml>
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=08:00:00
#SBATCH --job-name=fi_eval
#SBATCH --output=logs/eval_%j.out

set -euo pipefail
TOML=${1:?usage: sbatch jobscripts/eval.sbatch <campaign.toml>}

# pin BOTH thread pools (unpinned threads silently corrupt timings — known footgun)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export JULIA_GLUE_THREADS=8

# juliaup launcher (on PATH on compute nodes); the lattice_scale manifest is resolved
# with Julia 1.12 = the juliaup `release` default. If your juliaup default ever moves
# off 1.12, pin with `julia +1.12` instead.
JULIA=/mnt/home/xgao1/.juliaup/bin/julia

cd /mnt/home/xgao1/work/four_index_integral_solver/codes/lattice_scale
mkdir -p logs
$JULIA --project -t $JULIA_GLUE_THREADS driver.jl "$TOML" eval

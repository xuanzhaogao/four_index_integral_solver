#!/bin/bash
# Solve phase. SUBMIT MANUALLY:  sbatch --nodes=<M> jobscripts/solve.sbatch <campaign.toml>
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=08:00:00
#SBATCH --job-name=fi_solve
#SBATCH --output=logs/solve_%j.out

set -euo pipefail
TOML=${1:?usage: sbatch jobscripts/solve.sbatch <campaign.toml>}

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
$JULIA --project -t $JULIA_GLUE_THREADS driver.jl "$TOML" solve

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
# OpenBLAS gets ONE thread. The 2.95x-at-64-threads figure this once cited comes from the
# Sec. 5.2 benchmark -- a SINGLE process at -t 64/96, where an equally sized BLAS pool
# oversubscribes the node against FMM3D and FINUFFT. It does NOT describe this pipeline: job
# 7012357 timed the campaign's own evaluation at 8 and 96 Julia threads crossed with OPENBLAS
# 1 and 96, and all four cells landed within 3% (231.9-239.6 s). Keeping 1 because it is
# right for the single-process case and measured harmless here.
export OPENBLAS_NUM_THREADS=1
export JULIA_GLUE_THREADS=8

# juliaup launcher (on PATH on compute nodes); the codes/article manifest (found by the bare
# --project from this directory) is resolved with the juliaup `release` default.
JULIA=/mnt/home/xgao1/.juliaup/bin/julia

cd /mnt/home/xgao1/work/four_index_integral_solver/codes/article/sec53_lattice
mkdir -p logs
$JULIA --project -t $JULIA_GLUE_THREADS driver.jl "$TOML" eval

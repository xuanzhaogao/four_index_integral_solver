#!/bin/bash
#
# Slurm submission script for the bilayer bare-Coulomb sweep.
# Submit from the graphene project root:
#
#     cd /mnt/home/xgao1/work/four_index_integral_solver/codes/graphene
#     sbatch bilayer/scripts/submit_bare_bilayer_graphene.slurm
#
# Wall-clock estimate: 10 unique orbital pairs × ~30 s/pair (bare TKM3D
# integral, source_tol = 1e-3). 60-min limit is generous overhead.

#SBATCH --job-name=bare_bilayer
#SBATCH --output=bilayer/logs/bare_bilayer_%j.out
#SBATCH --error=bilayer/logs/bare_bilayer_%j.err
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=01:00:00

set -euo pipefail

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-48}

mkdir -p bilayer/logs

echo "host: $(hostname)"
echo "date: $(date -Is)"
echo "julia: $(which julia)"
julia --version
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS JULIA_NUM_THREADS=$JULIA_NUM_THREADS"
echo ""

srun --cpus-per-task=$SLURM_CPUS_PER_TASK \
    julia --project=. bilayer/scripts/bare_bilayer_graphene.jl

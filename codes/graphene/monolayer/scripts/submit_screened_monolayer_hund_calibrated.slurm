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

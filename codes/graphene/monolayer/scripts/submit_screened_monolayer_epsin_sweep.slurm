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

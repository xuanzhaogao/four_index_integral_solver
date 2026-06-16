#!/bin/bash
# exp65 multi-RHS benchmark — center orbital + lattice-neighbor pair densities,
# block solve vs sequential, full cutoff sweep (K up to ~19). One exclusive 96-core
# genoa node. Eval runs at 96 threads ON PURPOSE (documents the FINUFFT/FFTW thread
# pathology); four_index_matrix is K separate ~50 s TKM calls, so the largest-K eval is
# ~16 min and the full 5-cutoff sweep is ~60-90 min -> 2 h walltime.
# Submit yourself:  sbatch slurm/pvf_multi_rhs.sbatch
#SBATCH --job-name=pvf_multi_rhs
#SBATCH --partition=gen
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_multi_rhs_%j.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_multi_rhs_%j.err

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
WT=/mnt/home/xgao1/codes/BoundaryIntegral.jl
J=/mnt/home/xgao1/.juliaup/bin/julia          # juliaup binary (never module load julia)
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96

echo "=== node $(hostname)  cores $(nproc)  $(date) ==="
$J --project=$WT $NR/exp65_orbital_bench/scripts/run_multi_rhs.jl
echo "=== DONE $(date) ==="

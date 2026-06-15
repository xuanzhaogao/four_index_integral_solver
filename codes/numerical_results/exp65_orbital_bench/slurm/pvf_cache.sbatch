#!/bin/bash
# PVF benchmark — cache_fft mode (precomputed fine-grid FFT + interp-only
# evaluation) vs the standard field path. One exclusive 96-core genoa node.
#SBATCH --job-name=pvf_cache
#SBATCH --partition=gen
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_cache_%j.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_cache_%j.err

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
WT=/mnt/home/xgao1/codes/BoundaryIntegral.jl
J=/mnt/home/xgao1/.juliaup/bin/julia          # juliaup binary (never module load julia)
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96

echo "=== node $(hostname)  cores $(nproc)  $(date) ==="
$J --project=$WT $NR/exp65_orbital_bench/scripts/bench_cache_fft.jl
echo "=== DONE $(date) ==="

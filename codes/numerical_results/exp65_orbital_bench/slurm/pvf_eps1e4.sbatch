#!/bin/bash
# Old vs new u_int with both at eps = 1e-4 (bench_eps1e4.jl) — separates the
# reuse benefit from the tolerance benefit.
#SBATCH --job-name=pvf_eps1e4
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_eps1e4_%j.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_eps1e4_%j.err

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
WT=/mnt/home/xgao1/codes/BoundaryIntegral.jl
J=/mnt/home/xgao1/.juliaup/bin/julia          # juliaup binary (never module load julia)
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96

echo "=== node $(hostname)  cores $(nproc)  $(date) ==="
$J --project=$WT $NR/exp65_orbital_bench/scripts/bench_eps1e4.jl
echo "=== DONE $(date) ==="

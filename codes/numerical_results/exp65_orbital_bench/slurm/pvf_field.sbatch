#!/bin/bash
# PVF benchmark — warm old-vs-new A/B (build + RHS + u_int, VolumeSource path
# included). One exclusive 96-core genoa node. Submit together with the other
# two via the one-liner in the README header of this dir.
#SBATCH --job-name=pvf_field
#SBATCH --partition=gen
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_field_%j.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_field_%j.err

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
WT=/mnt/home/xgao1/codes/BoundaryIntegral.jl
J=/mnt/home/xgao1/.juliaup/bin/julia          # juliaup binary (never module load julia)
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96

echo "=== node $(hostname)  cores $(nproc)  $(date) ==="
$J --project=$WT $NR/exp65_orbital_bench/scripts/bench_field_path.jl
echo "=== DONE $(date) ==="

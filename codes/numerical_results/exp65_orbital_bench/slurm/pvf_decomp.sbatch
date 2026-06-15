#!/bin/bash
# PVF benchmark — clean decomposition of the old ltkm3dc u_int call
# (type-1 / scale / type-2 / forced upsampfac). On an EXCLUSIVE node the
# decomposed pieces should sum to ~= the full call (on contended worker7020
# they summed to ~2x, which flagged the contention).
#SBATCH --job-name=pvf_decomp
#SBATCH --partition=gen
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_decomp_%j.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_decomp_%j.err

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
WT=/mnt/home/xgao1/codes/BoundaryIntegral.jl
J=/mnt/home/xgao1/.juliaup/bin/julia          # juliaup binary (never module load julia)
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96

echo "=== node $(hostname)  cores $(nproc)  $(date) ==="
$J --project=$WT $NR/exp65_orbital_bench/scripts/diag_ltkm_decomp.jl
echo "=== DONE $(date) ==="

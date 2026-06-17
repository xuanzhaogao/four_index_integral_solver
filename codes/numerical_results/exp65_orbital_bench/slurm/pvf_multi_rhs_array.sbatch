#!/bin/bash
# exp65 multi-RHS benchmark as a SLURM JOB ARRAY — one task per neighbor cutoff (5 tasks),
# each on its OWN exclusive 96-core genoa node, run in parallel. Each task runs a SINGLE
# cutoff (ORBBENCH_CUTOFF) and writes its own per-cutoff record:
#   data/raw/multi_rhs_K<K>_cut<cutoff>.jls
# No sequential baseline; the block path is the production path. FINUFFT FFT is capped to 16
# threads IN-FUNCTION (TKM3D ltkm3dc) to dodge the 96-thread FFTW plan pathology, while
# OMP=96 keeps the FMM (lfmm3d) full-width.
#
# Submit (you):   sbatch slurm/pvf_multi_rhs_array.sbatch
# After ALL 5 tasks finish, gather + plot on a login node / workstation:
#   /mnt/home/xgao1/.juliaup/bin/julia --project=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results \
#       /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/scripts/plot_scaling.jl
#
# Submitted to the ccm partition on genoa nodes (per the SCC cluster guidance).
#SBATCH --job-name=pvf_mrhs
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --array=0-4
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_mrhs_%A_%a.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_mrhs_%A_%a.err

CUTOFFS=(0.0 1.5 2.5 2.9 3.9)          # array index -> neighbor cutoff (K = 1, 4, 10, 13, 19)
export ORBBENCH_CUTOFF=${CUTOFFS[$SLURM_ARRAY_TASK_ID]}

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
WT=/mnt/home/xgao1/codes/BoundaryIntegral.jl
J=/mnt/home/xgao1/.juliaup/bin/julia   # juliaup binary (never module load julia)
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96              # full FMM threads; FINUFFT FFT capped in-function (TKM3D)

echo "=== array task $SLURM_ARRAY_TASK_ID  cutoff=$ORBBENCH_CUTOFF  node $(hostname)  $(date) ==="
$J --project=$WT $NR/exp65_orbital_bench/scripts/run_multi_rhs.jl
echo "=== DONE cutoff=$ORBBENCH_CUTOFF  $(date) ==="

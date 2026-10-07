#!/bin/bash
# exp66 multicube multi-RHS as a SLURM JOB ARRAY — one task per neighbor cutoff (5 tasks),
# each on its OWN exclusive 96-core genoa node, run in parallel. Each task runs a SINGLE
# cutoff (MULTICUBE_CUTOFF) and writes its own per-cutoff record:
#   data/raw/multicube_K<K>_cut<cutoff>.jls
# Substrate = two 90x90x90 cubes (eps 11.9 | 3.9, Si/SiO2) + slab (eps 10); source = the
# central graphene p_z orbital + neighbors -> K pair densities, block-solved (exp65 pipeline).
# 96 threads everywhere: TKM3D's in-function FINUFFT nthreads cap (16) avoids the FFTW
# 96-thread plan pathology in BOTH the interface build and the eval (no DUCC override needed).
#
# Submit (you):   sbatch sec52_performance/slurm/run_multicube_array.sbatch
# After ALL 5 tasks finish, gather the per-cutoff .jls (cat the .out files for the summary).
#
#SBATCH --job-name=multicube
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --array=0-4
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/article/sec52_performance/logs/slurm_multicube_%A_%a.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/article/sec52_performance/logs/slurm_multicube_%A_%a.err

CUTOFFS=(0.0 1.5 2.5 2.9 3.9)          # array index -> neighbor cutoff (K = 1, 4, 10, 13, 19)
export MULTICUBE_CUTOFF=${CUTOFFS[$SLURM_ARRAY_TASK_ID]}

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/article
J=/mnt/home/xgao1/.juliaup/bin/julia   # juliaup binary (never module load julia)
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96              # full threads for FMM; TKM3D caps FINUFFT at 16 internally

echo "=== array task $SLURM_ARRAY_TASK_ID  cutoff=$MULTICUBE_CUTOFF  node $(hostname)  $(date) ==="
$J --project=$NR $NR/sec52_performance/scripts/run_multicube.jl
echo "=== DONE cutoff=$MULTICUBE_CUTOFF  $(date) ==="

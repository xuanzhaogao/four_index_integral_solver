#!/bin/bash
# exp66 multicube multi-RHS — SINGLE cutoff (K=1, onsite only) on one exclusive 96-core
# genoa node. Two 90x90x90 substrate cubes (eps 11.9 | 3.9, Si/SiO2) + slab (eps 10);
# source = the central graphene p_z orbital (+ neighbors within MULTICUBE_CUTOFF), as exp65.
# For the full K-sweep as parallel array tasks use run_multicube_array.sbatch instead.
#
# Runs at full 96 threads: TKM3D's in-function FINUFFT nthreads cap (default 16) avoids
# the FFTW 96-thread plan pathology in BOTH the interface build and the eval (no DUCC).
#
# Submit (you):  sbatch sec52_performance/slurm/run_multicube.sbatch
#SBATCH --job-name=multicube
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/article/sec52_performance/logs/slurm_multicube_%j.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/article/sec52_performance/logs/slurm_multicube_%j.err

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/article
J=/mnt/home/xgao1/.juliaup/bin/julia       # juliaup binary (never module load julia)
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96                   # full threads for the FMM; TKM3D caps FINUFFT at 16 internally
export MULTICUBE_CUTOFF=0.0                 # single K=1 (onsite only); raise for more neighbors

echo "=== node $(hostname)  cores $(nproc)  $(date) ==="
$J --project=$NR $NR/sec52_performance/scripts/run_multicube.jl
echo "=== DONE $(date) ==="

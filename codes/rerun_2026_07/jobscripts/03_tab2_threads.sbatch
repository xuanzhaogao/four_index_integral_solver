#!/bin/bash
# RERUN 03 -- Table 2 (strong scaling of the block solve, article Sec. 6.3).
#
# One exclusive node per thread count, run as an array. Threads pinned:
# unpinned, the OS packs low thread counts onto a single CCD and flattens the
# 1->2 step (0.95x unpinned vs 1.23x pinned on identical hardware).
#
# Timing-only itmax=1: the full block-solve time is reconstructed offline by
# plot_scaling.jl as t_RHS + t_1iter * niter_ref, with niter_ref read from the
# converged K=1 record that job 04 produces. RUN JOB 04 BEFORE PLOTTING -- with
# the preconditioner on, the published niter_ref = 36 is stale.
#
# Writes exp66_multicube/data_v2/raw_threads/.
#
# Submit (you):  sbatch 03_tab2_threads.sbatch
#
#SBATCH --job-name=rr03_threadpin
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --array=0-7
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/logs/rr03_threadpin_%A_%a.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/logs/rr03_threadpin_%A_%a.err

set -uo pipefail
NTHREADS_LIST=(1 2 4 8 16 32 64 96)
N=${NTHREADS_LIST[$SLURM_ARRAY_TASK_ID]}

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
J=/mnt/home/xgao1/.juliaup/bin/julia
source /mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/jobscripts/_provenance.sh

export RERUN_TAG=_v2
export PRECONDITION=1
export CORRECT_EDGES=1
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export JULIA_NUM_THREADS=$N
export OMP_NUM_THREADS=$N
export OPENBLAS_NUM_THREADS=$N
export MULTICUBE_CUTOFF=0.0          # K=1
export MULTICUBE_GMRES_ITMAX=1       # timing-only single iteration
export MULTICUBE_RUNTAG=pin_nt${N}
provenance

echo "=== array task $SLURM_ARRAY_TASK_ID  nthreads=$N  $(date '+%F %T') ==="
echo "OMP_PROC_BIND=$OMP_PROC_BIND  OMP_PLACES=$OMP_PLACES"
$J --project=$NR --startup-file=no $NR/exp66_multicube/scripts/run_multicube.jl
echo "=== DONE nthreads=$N  $(date '+%F %T') ==="

#!/bin/bash
# exp66 multicube THREAD-SCALING sweep (left panel of the redesigned multi-RHS figure):
# K=1 (onsite only, MULTICUBE_CUTOFF=0.0) on the SAME multicube system as the K-sweep
# (run_multicube_array.sbatch), run at nthreads in {1,2,4,8,16,32,64,96} as a SLURM JOB
# ARRAY — one task per thread count, each on its OWN exclusive 96-core genoa node (only
# nthreads of the 96 cores are used; the rest sit idle so there is no cross-run contention
# and the low-thread-count timings are honest).
#
# TIMING-ONLY block solve: MULTICUBE_GMRES_ITMAX=1 runs a single block-GMRES iteration
# instead of converging (which would take hours at nthreads=1). The full block-solve time
# is reconstructed OFFLINE in the plot script as t_1iter * niter_ref, where niter_ref is the
# GMRES iteration count from the real converged K=1 run at 96 threads already on disk
# (data/raw/multicube_K1_cut0.0.jls, from run_multicube_array.sbatch) — GMRES's iteration
# count is a property of the operator/preconditioning, not of the thread count, so this
# reference transfers across the sweep. precompute/pottrg/eval are timed for real at every
# thread count (they are cheap enough not to need the same trick).
#
# Each task writes data/raw_threads/multicube_nt<N>.jls (MULTICUBE_RUNTAG=nt<N> keeps this
# completely separate from the K-sweep's data/raw/ records and their gather regex).
#
# Submit (you):   sbatch exp66_multicube/slurm/run_multicube_threads_array.sbatch
# After all 8 tasks finish: julia --project=. plot_scripts/plot_multicube_scaling.jl
#
#SBATCH --job-name=multicube_threads
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --array=0-7
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp66_multicube/logs/slurm_multicube_threads_%A_%a.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp66_multicube/logs/slurm_multicube_threads_%A_%a.err

NTHREADS_LIST=(1 2 4 8 16 32 64 96)     # array index -> thread count
N=${NTHREADS_LIST[$SLURM_ARRAY_TASK_ID]}

export MULTICUBE_CUTOFF=0.0             # K=1 (onsite only) -- same system as the K-sweep
export MULTICUBE_GMRES_ITMAX=1          # timing-only: 1 block-GMRES iteration, not converged
export MULTICUBE_RUNTAG=nt${N}          # -> data/raw_threads/multicube_nt<N>.jls

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
J=/mnt/home/xgao1/.juliaup/bin/julia    # juliaup binary (never module load julia)
export JULIA_NUM_THREADS=$N
export OMP_NUM_THREADS=$N               # FMM3D threads
export OPENBLAS_NUM_THREADS=$N          # pin BLAS too -- unpinned threads corrupt timings

echo "=== array task $SLURM_ARRAY_TASK_ID  nthreads=$N  node $(hostname)  $(date) ==="
$J --project=$NR $NR/exp66_multicube/scripts/run_multicube.jl
echo "=== DONE nthreads=$N  $(date) ==="

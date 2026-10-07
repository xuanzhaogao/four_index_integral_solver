#!/bin/bash
# Pinned thread-scaling sweep for the multi-RHS block solve (K=1) as a SLURM JOB ARRAY -- one
# task per thread count, each on its OWN exclusive genoa node, run in parallel. This is the
# earlier run_multicube_threads_array.sbatch WITH thread pinning added: unpinned, the OS packed
# low thread counts onto a single CCD (shared L3 + memory-bandwidth link), which flattened the
# 1->2 step (0.95x unpinned vs 1.23x pinned on identical HW). OMP_PROC_BIND=spread /
# OMP_PLACES=cores spread the N OpenMP threads across cores -> across CCDs/sockets at low N, so
# every thread count is placement-clean. (Node-to-node variance is negligible -- the 1-thread
# cost matched to <0.5% across two different genoa nodes -- so a per-node array is fine here.)
#
# Timing-only itmax=1: validated that the single-iteration block-GMRES time equals the steady
# per-iteration time, i.e. no warm-up/JIT tail. Full block-solve time is reconstructed offline
# as t_RHS + t_1iter * niter_ref (niter from the converged K=1 @ 96t run in data/raw/).
#
# Writes data/raw_threads/multicube_pin_nt<N>.jls (distinct from the unpinned nt<N> records;
# the plot script prefers the pinned set when present).
#
# Submit (you):  sbatch sec52_performance/slurm/run_multicube_threads_pinned.sbatch
# After all 8 tasks finish: julia --project=. sec52_performance/plot_scaling.jl
#
#SBATCH --job-name=mc_threadpin
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --array=0-7
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/article/sec52_performance/logs/slurm_multicube_pin_%A_%a.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/article/sec52_performance/logs/slurm_multicube_pin_%A_%a.err

set -uo pipefail
NTHREADS_LIST=(1 2 4 8 16 32 64 96)     # array index -> thread count
N=${NTHREADS_LIST[$SLURM_ARRAY_TASK_ID]}

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/article
J=/mnt/home/xgao1/.juliaup/bin/julia    # juliaup binary (never module load julia)

# pin + spread OpenMP threads across cores (=> across CCDs/sockets at low counts)
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export JULIA_NUM_THREADS=$N
export OMP_NUM_THREADS=$N
export OPENBLAS_NUM_THREADS=$N
export MULTICUBE_CUTOFF=0.0             # K=1, same system as the K-sweep
export MULTICUBE_GMRES_ITMAX=1          # timing-only single iteration (no JIT/warm-up tail)
export MULTICUBE_RUNTAG=pin_nt${N}

echo "=== array task $SLURM_ARRAY_TASK_ID  nthreads=$N  node $(hostname)  $(date) ==="
echo "OMP_PROC_BIND=$OMP_PROC_BIND  OMP_PLACES=$OMP_PLACES"
$J --project=$NR --startup-file=no $NR/sec52_performance/scripts/run_multicube.jl
echo "=== DONE nthreads=$N  $(date) ==="

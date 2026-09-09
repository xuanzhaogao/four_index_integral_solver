#!/bin/bash
# Sec. 5.2 pinned thread-scaling sweep (K=1) on the Sec. 5.3 SYSTEM. Same as
# rerun_2026_07/jobscripts/03_tab2_threads.sbatch except MULTICUBE_GEOM=sec53 (x3 Si|SiO2
# cubes, graphene slab at eps = 2.4). One task per thread count, each on its own exclusive
# genoa node. Threads pinned: unpinned, the OS packs low thread counts onto a single CCD and
# flattens the 1->2 step (0.95x unpinned vs 1.23x pinned on identical hardware).
#
# Timing-only itmax=1: the full block-solve time is reconstructed offline by plot_scaling.jl
# as t_RHS + t_1iter * niter_ref, with niter_ref read from the converged K=1 record that
# run_multicube_array_sec53.sbatch produces. RUN THE K-SWEEP FIRST -- niter changes with the
# slab permittivity (lower contrast at eps = 2.4), so the _v2 reference is stale here.
#
# Writes exp66_multicube/data_sec53/raw_threads/multicube_pin_nt<N>.jls.
#
# Node image pinned to rocky8 (genoa,rocky8 + --reservation=rocky8): those nodes sit inside
# the `rocky8` SPEC_NODES reservation, so the reservation is required to schedule on them.
# One image for the whole sweep is deliberate -- the _v2 thread sweep split across images
# (nt 1/2/4/8 on worker7136/7205/7206/7207 = rocky8, nt 16/32/64/96 on
# worker7297/7298/7299/7391 = rocky9), so its speedup-vs-1-thread compares two platforms.
#
# SUBMIT-TIME TRAP: a stale SBATCH_RESERVATION=rocky9 lives in the environment (that
# reservation was deleted once rocky9 became the default), and per sbatch(1) "environment
# variables will override any options set in a batch script" -- so it defeats the
# --reservation directive below and the job fails with "Requested reservation is invalid".
# Override it on the command line (which beats the environment):
#   sbatch --reservation=rocky8 <this script>
# or clear it first:
#   unset SBATCH_RESERVATION SALLOC_RESERVATION SLURM_RESERVATION
#
# Submit (you):  sbatch --reservation=rocky8 exp66_multicube/slurm/run_multicube_threads_pinned_sec53.sbatch
#
#SBATCH --job-name=mc53_threadpin
#SBATCH --partition=ccm
#SBATCH --constraint=genoa,rocky8
#SBATCH --reservation=rocky8
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --array=0-7
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp66_multicube/logs/slurm_mc53_pin_%A_%a.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp66_multicube/logs/slurm_mc53_pin_%A_%a.err

set -uo pipefail
NTHREADS_LIST=(1 2 4 8 16 32 64 96)     # array index -> thread count
N=${NTHREADS_LIST[$SLURM_ARRAY_TASK_ID]}

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
J=/mnt/home/xgao1/.juliaup/bin/julia    # juliaup binary (never module load julia)
source /mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/jobscripts/_provenance.sh

# Overridable: `export` here would otherwise defeat --export=ALL,RERUN_TAG=... and send
# every rerun to the same tree, silently overwriting earlier records under identical
# filenames.
export RERUN_TAG=${RERUN_TAG:-_sec53}
export MULTICUBE_GEOM=sec53
# COMBINED solve + evaluation measurement. MULTICUBE_EVAL=model evaluates the batch at the
# FIXED target set of the Sec. 5.3 model (union of all 198 orbitals' quadrature points,
# N_p = 7,168,390) after a real solve, so one record carries prepare, block solve AND the
# evaluation the ERI assembly actually performs -- same interface, same accuracy, same job.
# Previously the two were separate runs (an EVAL_ONLY=1 sweep at a different accuracy), so a
# table combining them mixed discretizations.
export MULTICUBE_EVAL=model
export MULTICUBE_TARGET_CACHE=$NR/exp66_multicube/data_sec53evalpilot/targets_model.jls
export PRECONDITION=1
export CORRECT_EDGES=1
# pin + spread OpenMP threads across cores (=> across CCDs/sockets at low counts)
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export JULIA_NUM_THREADS=$N
export OMP_NUM_THREADS=$N
export OPENBLAS_NUM_THREADS=$N
export MULTICUBE_CUTOFF=0.0             # K=1, same system as the K-sweep
export MULTICUBE_GMRES_ITMAX=1          # timing-only single iteration (no JIT/warm-up tail)
export MULTICUBE_RUNTAG=pin_nt${N}
provenance

[ -f "$MULTICUBE_TARGET_CACHE" ] || { echo "FATAL: target cache missing: $MULTICUBE_TARGET_CACHE"; exit 1; }

echo "=== array task $SLURM_ARRAY_TASK_ID  nthreads=$N  $(date '+%F %T') ==="
echo "OMP_PROC_BIND=$OMP_PROC_BIND  OMP_PLACES=$OMP_PLACES"
$J --project=$NR --startup-file=no $NR/exp66_multicube/scripts/run_multicube.jl
echo "=== DONE nthreads=$N  $(date '+%F %T') ==="

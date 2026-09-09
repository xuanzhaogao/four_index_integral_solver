#!/bin/bash
# EVALUATION benchmark, THREAD sweep -- the eval counterpart of Table 2's strong scaling.
# K = 1 fixed, one thread count per task, each on its own exclusive rocky8 node.
#
# NO SOLVE. MULTICUBE_EVAL_ONLY=1 builds the interface (pottrg and the near corrections need
# it, and N fixes the matvec cost) but skips the corrected LHS operator, the RHS assembly and
# block GMRES, filling sigma with a fixed-seed random field instead. Evaluation cost depends on
# sigma's SHAPE, not its values, so the timing is exact while the solve -- which at 1 thread
# costs more than the measurement -- is never paid for. Physics outputs of these records are
# NaN by construction.
#
# Targets are the FIXED model set: the union of the quadrature points of all 198 Sec. 5.3
# orbitals, N_p = 7,168,390. Fixed is the point -- an earlier attempt used the K-batch's own
# grid, which grows with K, and its flat per-source time was two effects cancelling.
#
# Writes data_sec53eval/raw_threads/multicube_evnt<N>.jls.
#
# Submit (you):  sbatch --reservation=rocky8 exp66_multicube/slurm/run_evalbench_threads_sec53.sbatch
#   (--reservation must be on the COMMAND LINE: a stale SBATCH_RESERVATION=rocky9 in the
#    environment overrides script directives, per sbatch(1))
#
#SBATCH --job-name=mc53_evthreads
#SBATCH --partition=ccm
#SBATCH --constraint=genoa,rocky8
#SBATCH --reservation=rocky8
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --array=0-7
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp66_multicube/logs/slurm_mc53_evthreads_%A_%a.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp66_multicube/logs/slurm_mc53_evthreads_%A_%a.err

set -uo pipefail
NTHREADS_LIST=(1 2 4 8 16 32 64 96)
N=${NTHREADS_LIST[$SLURM_ARRAY_TASK_ID]}

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
J=/mnt/home/xgao1/.juliaup/bin/julia
source /mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/jobscripts/_provenance.sh

export RERUN_TAG=${RERUN_TAG:-_sec53eval}   # override to keep pre/post-fix trees apart
export MULTICUBE_GEOM=sec53
export MULTICUBE_EVAL=model
export MULTICUBE_EVAL_ONLY=1
export MULTICUBE_TARGET_CACHE=$NR/exp66_multicube/data_sec53evalpilot/targets_model.jls
export MULTICUBE_CUTOFF=0.0             # K = 1
# RUNTAG set after PIN is known (below)
export PRECONDITION=1
export CORRECT_EDGES=1
# Thread placement. EVALBENCH_PIN=1 (default) spreads threads across cores, which is what the
# SOLVE thread sweep needs: unpinned, the OS packs low thread counts onto one CCD and flattens
# the 1->2 step. For the EVALUATION it appears to do the opposite -- at 96 threads, K=1,
# N_p = 7.17M the pinned run took 150.5 s against 96.6 s (K sweep) and 93.2 s (pilot 7005239),
# both unpinned and agreeing to 4%. The far-field lfmm3d over 7.17M targets is memory-bound and
# is deliberately left on the full OpenMP pool (only FINUFFT is capped, volume_field.jl:118),
# so spreading threads away from their data is the suspect. EVALBENCH_PIN=0 runs the control.
PIN=${EVALBENCH_PIN:-0}   # measured: no effect on the eval (<=3%); on for the SOLVE sweep only
if [ "$PIN" = "1" ]; then
  export OMP_PROC_BIND=spread
  export OMP_PLACES=cores
  TAGSUF=""
else
  unset OMP_PROC_BIND OMP_PLACES
  TAGSUF="_nopin"
fi
export JULIA_NUM_THREADS=$N
export OMP_NUM_THREADS=$N
# EVALBENCH_OPENBLAS overrides the OpenBLAS pool. It is set to $N here because this script
# was modelled on the SOLVE thread sweep, where BLAS carries real work. The evaluation has
# almost no BLAS, and at 96 threads the four measurements split exactly on this variable:
# 93.2 s (pilot) and 96.6 s (K sweep), both with it UNSET, against 150.5 s and 154.7 s with it
# at 96 -- consistent with OpenBLAS oversubscribing alongside FMM3D's own 96 OpenMP threads.
# Set EVALBENCH_OPENBLAS=1 to test that; leave unset for the historical behaviour.
# Default 1, NOT $N. The evaluation has almost no BLAS, and OpenBLAS spawning $N
# threads alongside FMM3D's own pool cost 2.95x at 64 threads and 1.56x at 96 --
# measured. $N was inherited from the SOLVE sweep, where BLAS carries real work.
export OPENBLAS_NUM_THREADS=${EVALBENCH_OPENBLAS:-1}
TAGSUF="${TAGSUF}${EVALBENCH_OPENBLAS:+_ob$EVALBENCH_OPENBLAS}"
export MULTICUBE_RUNTAG=evnt${N}${TAGSUF}
provenance

[ -f "$MULTICUBE_TARGET_CACHE" ] || { echo "FATAL: target cache missing: $MULTICUBE_TARGET_CACHE"; exit 1; }

echo "=== eval thread sweep  nthreads=$N  pin=$PIN  $(date '+%F %T') ==="
$J --project=$NR --startup-file=no $NR/exp66_multicube/scripts/run_multicube.jl
echo "=== DONE nthreads=$N  $(date '+%F %T') ==="

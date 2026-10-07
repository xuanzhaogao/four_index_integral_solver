#!/bin/bash
# Sec. 5.2 K-sweep on the Sec. 5.3 SYSTEM -- one task per neighbor cutoff, each on its own
# exclusive genoa node. Same as rerun_2026_07/jobscripts/04_tab3_ksweep.sbatch except for the
# substrate: MULTICUBE_GEOM=sec53 replaces the L = 90, eps_slab = 10 demo substrate with the
# x3 converged Si|SiO2 cubes (L = 270) and the graphene slab at its cRPA eps = 2.4, i.e. the
# same system the article's lattice section reports. Discretization and the sweep are unchanged.
#
# Writes sec52_performance/data_sec53/raw/ -- data/ (published) and data_v2/ (July rerun) are
# left byte-identical for comparison.
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
# Submit (you):  sbatch --reservation=rocky8 sec52_performance/slurm/run_multicube_array_sec53.sbatch
# After all 7 tasks finish:
#   RERUN_TAG=_sec53 julia --project=. sec52_performance/plot_scaling.jl
#
#SBATCH --job-name=mc53_ksweep
#SBATCH --partition=ccm
#SBATCH --constraint=genoa,rocky8
#SBATCH --reservation=rocky8
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --array=0-9
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/article/sec52_performance/logs/slurm_mc53_%A_%a.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/article/sec52_performance/logs/slurm_mc53_%A_%a.err

set -uo pipefail
CUTOFFS=(0.0 1.5 2.5 2.9 3.9 4.5 5.0 5.5 6.0 6.5)   # -> K = 1, 4, 10, 13, 19, 25, 31, 37, 40, 46
# K = 37/40/46 used to live in a separate extraK script, added when the K <= 31 sweep had
# already run and the curve needed extending past the amortization peak. Kept apart they
# were two scripts to hold in sync for one sweep; K = 46 is the last point that fits the
# graphene neighbour shells below the node memory ceiling, so the list is complete here.
export MULTICUBE_CUTOFF=${CUTOFFS[$SLURM_ARRAY_TASK_ID]}

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/article
J=/mnt/home/xgao1/.juliaup/bin/julia    # juliaup binary (never module load julia)
source /mnt/home/xgao1/work/four_index_integral_solver/codes/article/common/_provenance.sh

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
export MULTICUBE_TARGET_CACHE=$NR/sec52_performance/data_sec53evalpilot/targets_model.jls
export PRECONDITION=1
export CORRECT_EDGES=1
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96               # TKM3D caps FINUFFT at 16 internally
provenance

[ -f "$MULTICUBE_TARGET_CACHE" ] || { echo "FATAL: target cache missing: $MULTICUBE_TARGET_CACHE"; exit 1; }

echo "=== array task $SLURM_ARRAY_TASK_ID  cutoff=$MULTICUBE_CUTOFF  $(date '+%F %T') ==="
$J --project=$NR --startup-file=no $NR/sec52_performance/scripts/run_multicube.jl
echo "=== DONE cutoff=$MULTICUBE_CUTOFF  $(date '+%F %T') ==="

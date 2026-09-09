#!/bin/bash
# Sec. 5.2 K-sweep with the EVALUATION BENCHMARK (MULTICUBE_EVAL=model): K solved sources
# evaluated at the FIXED target set of the Sec. 5.3 model -- the union of the quadrature points
# of all 198 orbitals, N_p = 7,168,390. One task per batch size, own exclusive rocky8 node.
#
# Sized from pilot job 7005239 (K = 1 and 4 in one process), which measured
#     t_eval(K) ~= 81.6 + 11.6 K  s
# i.e. a shared, K-independent cost of 81.6 s (the corrected pottrg map over N_p targets plus
# the near/far setup) and 11.6 s per source. Per-source cost therefore falls as 81.6/K + 11.6:
# 93.2 s at K = 1, 32.0 s at K = 4, extrapolating to ~13.4 s at K = 46 -- a 7x reduction
# saturating at the field-evaluation floor. That is the multi-RHS amortization for the
# EVALUATION stage, and it only becomes visible with N_p held fixed.
#
# MEMORY IS NOT THE CONSTRAINT, contrary to the earlier worry. The pilot peaked at 41.7 GB
# (K = 1) and 119.0 GB (K = 4) against 42.8 and 118.7 for the row-protocol sweep: pottrg over
# 7.17M targets costs 0.7 GB. Peak RSS stays solve-dominated (26.2 GB per RHS), so K = 46 lands
# near the row sweep's 1223 GB of a 1502 GiB node and the full sweep is reachable.
#
# Expected wall time at K = 46: ~19 s precompute + ~1664 s solve + ~620 s eval ~= 38 min.
#
# THE TARGET CACHE IS REUSED, NOT REBUILT. Building it costs ~165 s and 86 GB; ten tasks each
# building it would waste that ten times AND race on the same output path. MULTICUBE_TARGET_CACHE
# points every task at the pilot's copy, which is read-only here.
#
# Writes exp66_multicube/data_sec53model/raw/. The row-protocol sweep in data_sec53/, the July
# rerun in data_v2/ and the submitted data/ are all untouched.
#
# Node image pinned to rocky8 to match every other sec53 measurement. SUBMIT-TIME TRAP: a stale
# SBATCH_RESERVATION=rocky9 in the environment overrides script directives (sbatch(1):
# "environment variables will override any options set in a batch script"), so pass the
# reservation on the command line:
#
# Submit (you):  sbatch --reservation=rocky8 exp66_multicube/slurm/run_multicube_array_modeleval_sec53.sbatch
# After it finishes:
#   RERUN_TAG=_sec53model julia --project=. exp66_multicube/scripts/article_tables.jl
#
#SBATCH --job-name=mc53_modeleval
#SBATCH --partition=ccm
#SBATCH --constraint=genoa,rocky8
#SBATCH --reservation=rocky8
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --array=0-9
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp66_multicube/logs/slurm_mc53_model_%A_%a.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp66_multicube/logs/slurm_mc53_model_%A_%a.err

set -uo pipefail
CUTOFFS=(0.0 1.5 2.5 2.9 3.9 4.5 5.0 5.5 6.0 6.5)   # -> K = 1, 4, 10, 13, 19, 25, 31, 37, 40, 46
export MULTICUBE_CUTOFF=${CUTOFFS[$SLURM_ARRAY_TASK_ID]}

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
J=/mnt/home/xgao1/.juliaup/bin/julia    # juliaup binary (never module load julia)
source /mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/jobscripts/_provenance.sh

export RERUN_TAG=_sec53model
export MULTICUBE_GEOM=sec53
export MULTICUBE_EVAL=model
export MULTICUBE_TARGET_CACHE=$NR/exp66_multicube/data_sec53evalpilot/targets_model.jls
export PRECONDITION=1
export CORRECT_EDGES=1
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96               # TKM3D caps FINUFFT at 16 internally
provenance

[ -f "$MULTICUBE_TARGET_CACHE" ] || { echo "FATAL: target cache missing: $MULTICUBE_TARGET_CACHE"; exit 1; }

echo "=== array task $SLURM_ARRAY_TASK_ID  cutoff=$MULTICUBE_CUTOFF  eval=model  $(date '+%F %T') ==="
$J --project=$NR --startup-file=no $NR/exp66_multicube/scripts/run_multicube.jl
echo "=== DONE cutoff=$MULTICUBE_CUTOFF  $(date '+%F %T') ==="

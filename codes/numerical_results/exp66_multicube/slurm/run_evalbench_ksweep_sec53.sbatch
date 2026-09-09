#!/bin/bash
# EVALUATION benchmark, K sweep -- the eval counterpart of Table 3. 96 threads, one batch size
# per task, each on its own exclusive rocky8 node.
#
# NO SOLVE. MULTICUBE_EVAL_ONLY=1 builds the interface but skips the corrected LHS operator,
# the RHS assembly and block GMRES, filling sigma with a fixed-seed random field. Evaluation
# cost depends on sigma's SHAPE, not its values, so the timing is exact and the solve is never
# paid for. Physics outputs of these records are NaN by construction.
#
# Targets are the FIXED model set: all 198 Sec. 5.3 orbitals' quadrature points,
# N_p = 7,168,390, the SAME for every K. That is what makes t_eval/K a per-source cost.
#
# From pilot 7005239, t_eval(K) ~= 81.6 + 11.6 K s at 96 threads: a shared, K-independent cost
# (the corrected pottrg map over N_p targets) plus a per-source field cost. Per-source time
# should therefore fall as 81.6/K + 11.6 -- 93 s at K = 1 to ~13 s at K = 46 -- and this sweep
# is what measures that curve. Expect ~11 min for the largest task.
#
# Writes data_sec53eval/raw_threads/multicube_evK<cutoff>.jls.
#
# Submit (you):  sbatch --reservation=rocky8 exp66_multicube/slurm/run_evalbench_ksweep_sec53.sbatch
#   (--reservation must be on the COMMAND LINE, per sbatch(1))
#
#SBATCH --job-name=mc53_evksweep
#SBATCH --partition=ccm
#SBATCH --constraint=genoa,rocky8
#SBATCH --reservation=rocky8
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --array=0-9
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp66_multicube/logs/slurm_mc53_evksweep_%A_%a.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp66_multicube/logs/slurm_mc53_evksweep_%A_%a.err

set -uo pipefail
CUTOFFS=(0.0 1.5 2.5 2.9 3.9 4.5 5.0 5.5 6.0 6.5)   # -> K = 1, 4, 10, 13, 19, 25, 31, 37, 40, 46
CUT=${CUTOFFS[$SLURM_ARRAY_TASK_ID]}

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
J=/mnt/home/xgao1/.juliaup/bin/julia
source /mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/jobscripts/_provenance.sh

export RERUN_TAG=${RERUN_TAG:-_sec53eval}   # override to keep pre/post-fix trees apart
export MULTICUBE_GEOM=sec53
export MULTICUBE_EVAL=model
export MULTICUBE_EVAL_ONLY=1
export MULTICUBE_TARGET_CACHE=$NR/exp66_multicube/data_sec53evalpilot/targets_model.jls
export MULTICUBE_CUTOFF=$CUT
export MULTICUBE_RUNTAG=evK${CUT}
export PRECONDITION=1
export CORRECT_EDGES=1
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96               # TKM3D caps FINUFFT at 16 internally
provenance

[ -f "$MULTICUBE_TARGET_CACHE" ] || { echo "FATAL: target cache missing: $MULTICUBE_TARGET_CACHE"; exit 1; }

echo "=== eval K sweep  cutoff=$CUT  $(date '+%F %T') ==="
$J --project=$NR --startup-file=no $NR/exp66_multicube/scripts/run_multicube.jl
echo "=== DONE cutoff=$CUT  $(date '+%F %T') ==="

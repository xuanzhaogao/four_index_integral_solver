#!/bin/bash
# RERUN 04 -- Table 3 (block multiple-right-hand-side solve, article Sec. 6.3).
#
# One task per neighbour cutoff; K = 1, 4, 10, 13, 19, 25, 31. The published run
# split these across two submissions (K<=19, then 25 and 31); they are combined
# here into one array so every point comes from one code state.
#
# MEMORY: K=31 peaked at 538 GB resident in the published run. These tasks need
# a full-memory exclusive node; do not co-schedule.
#
# Also produces the converged K=1 record that job 03's reconstruction reads as
# niter_ref, so run this one even if you only care about Table 2.
#
# Writes exp66_multicube/data_v2/raw/.
#
# Submit (you):  sbatch 04_tab3_ksweep.sbatch
#
#SBATCH --job-name=rr04_ksweep
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --array=0-6
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/logs/rr04_ksweep_%A_%a.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/logs/rr04_ksweep_%A_%a.err

set -uo pipefail
CUTOFFS=(0.0 1.5 2.5 2.9 3.9 4.5 5.0)   # -> K = 1, 4, 10, 13, 19, 25, 31
export MULTICUBE_CUTOFF=${CUTOFFS[$SLURM_ARRAY_TASK_ID]}

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
J=/mnt/home/xgao1/.juliaup/bin/julia
source /mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/jobscripts/_provenance.sh

export RERUN_TAG=_v2
export PRECONDITION=1
export CORRECT_EDGES=1
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96            # TKM3D caps FINUFFT at 16 internally
provenance

echo "=== array task $SLURM_ARRAY_TASK_ID  cutoff=$MULTICUBE_CUTOFF  $(date '+%F %T') ==="
$J --project=$NR --startup-file=no $NR/exp66_multicube/scripts/run_multicube.jl
echo "=== DONE cutoff=$MULTICUBE_CUTOFF  $(date '+%F %T') ==="

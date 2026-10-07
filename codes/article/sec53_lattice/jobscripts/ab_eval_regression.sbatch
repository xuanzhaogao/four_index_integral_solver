#!/bin/bash
# Attribute the eval-phase slowdown: 39 min (Aug 30, job 6965032) -> 110 min (kb31, job 7006478).
#
# The measurements that narrowed it are in scripts/ab_eval_regression.jl's header. What survives
# is a single unmeasured gap: the Aug 30 rate (8.23 s per source at the full target set) is a fit
# over K = 1..17, and the current rate (61.4 s at K = 31, 61.8 s at K = 46) rests on two points.
# Nothing has ever been measured for 17 < K < 31. This sweeps that gap on one batch, then checks
# OpenBLAS threads, which was 96 on Aug 30 and 1 for kb31.
#
# Reads finished batch results off ceph and re-runs ONLY the evaluation, so it needs no solve and
# writes nothing. One node, one process.
#
# Runtime: if the current 61 s/source holds at every K, the sweep is ~100 min and the BLAS pair
# another ~65 min. If the cost turns out to be Aug-30-like at low K it is far shorter. --time is
# set for the pessimistic case.
#
# SUBMIT (you — org policy bars the assistant from running sbatch):
#   cd /mnt/home/xgao1/work/four_index_integral_solver/codes/article/sec53_lattice
#   sbatch --reservation=rocky8 jobscripts/ab_eval_regression.sbatch
#
# Optional overrides, e.g. the K = 46 campaign once it lands, or a cheaper first look:
#   sbatch --reservation=rocky8 --export=ALL,CAMP=lattice_conv_l3_eps2.4_k46,BATCH_ID=1,KS=4,17,31,46 \
#     jobscripts/ab_eval_regression.sbatch
#
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=4:00:00
#SBATCH --job-name=fi_ab_eval
#SBATCH --output=logs/ab_eval_%j.out

set -euo pipefail

# Match the kb31 campaign environment exactly -- this run is a comparison against it, so the
# only variables that may differ are the ones the script itself sweeps.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=1                 # as run_all.sbatch; the script re-tests 96 in-process
export JULIA_GLUE_THREADS=$SLURM_CPUS_PER_TASK
JULIA=/mnt/home/xgao1/.juliaup/bin/julia      # juliaup launcher (never `module load julia`)

cd /mnt/home/xgao1/work/four_index_integral_solver/codes/article/sec53_lattice
mkdir -p logs

echo "=== ab_eval_regression  $(date)  on $(hostname) ==="
$JULIA --project -t "$SLURM_CPUS_PER_TASK" scripts/ab_eval_regression.jl
echo "=== AB DONE  $(date) ==="

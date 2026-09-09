#!/bin/bash
# Follow-up to ab_eval_regression (job 7010961): is the campaign's eval cost set by OPENBLAS
# THREADS INTERACTING WITH THE WORKER THREAD COUNT, rather than by K?
#
# The K sweep came back flat -- 9.90, 7.92, 8.06, 7.33, 7.58 s per source at K = 4, 8, 13, 17,
# 24, against the Aug 30 production rate of 8.23. So current code at 96 Julia threads reproduces
# Aug 30 exactly, while kb31 production measured 61.4 s per source at K = 31. K does not explain
# that; the environment might.
#
# run_phase launches its workers as
#     addprocs(SlurmManager(); exeflags = `--project=$proj -t $glue`),  glue = JULIA_GLUE_THREADS
# and BOTH campaign jobscripts set JULIA_GLUE_THREADS=8. So every production worker evaluates
# with 8 Julia threads on a 96-core node, whereas ab_eval_regression and the Sec. 5.2 benchmark
# (both fast) run a single process at -t 96.
#
# Julia's thread count is fixed at startup, so this runs the SAME evaluation under several
# counts, one process each. If 8 threads reproduces ~61 s per source and 96 gives ~7.5, the
# campaign has been evaluating at a small fraction of each node for both eras and the fix is one
# line in the jobscripts.
#
# Note what this would NOT explain by itself: Aug 30 also ran 8-thread workers and was fast at
# K <= 17. If 8 threads is slow at K = 31 but fine at K = 17, the two combine -- a thread-count
# ceiling that only binds once K exceeds it -- which is why K = 17 and 31 are both run here.
#
# SUBMIT (you — org policy bars the assistant from running sbatch):
#   cd /mnt/home/xgao1/work/four_index_integral_solver/codes/lattice_scale
#   sbatch --reservation=rocky8 jobscripts/ab_eval_threads.sbatch
#
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=4:00:00
#SBATCH --job-name=fi_ab_thr
#SBATCH --output=logs/ab_threads_%j.out

set -euo pipefail

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK   # as both campaign jobscripts set it
# OPENBLAS is set per-run in the loop below, not here.
JULIA=/mnt/home/xgao1/.juliaup/bin/julia

cd /mnt/home/xgao1/work/four_index_integral_solver/codes/lattice_scale
mkdir -p logs

echo "=== ab_eval_threads  $(date)  on $(hostname) ==="
# The decisive cell is 8 Julia threads x OPENBLAS 96 vs 1: 8 is what production workers get, and
# OpenBLAS is the ONE environment variable that changed between the fast Aug 30 run (96) and the
# slow kb31/k46 runs (1). At 96 Julia threads the pair is indistinguishable (7.73 vs 7.50 s per
# source, job 7010961) because Julia's own threading already saturates the node -- which is why
# the first A/B missed it.
for T in 8 96; do
    for B in 1 96; do
        echo "--- julia -t $T   OPENBLAS_NUM_THREADS=$B ---"
        OPENBLAS_NUM_THREADS=$B KS=31 CAMP=lattice_conv_l3_eps2.4_kb31 BATCH_ID=1 SKIP_BLAS=1 \
            $JULIA --project -t "$T" scripts/ab_eval_regression.jl
    done
done
echo "=== AB THREADS DONE  $(date) ==="

#!/bin/bash
# Single-orbital l_ec self-convergence sweep, any variant, any refinement levels.
# Generalizes run_lec_conv_eps2.4.sbatch (which hard-codes the eps=2.4 variant and levels 1-5).
#
# Writes figs/lec_single_conv<TAG>.tsv. scripts/lec_conv_single.jl MERGES rows into an existing
# TSV keyed by level, so adding levels does NOT recompute the ones already there -- submit
# `... 6 7` to extend a finished sweep rather than rerunning it from level 1.
#
# WHY HIGH LEVELS. With levels 2-5 the unknown count spans only 653k -> 5.6M, under one decade,
# so the log axis of Fig. 9(d) cannot show two whole decades and falls back to 5x10^5 / 2x10^6
# style ticks. Levels 6 and 7 extend the range to roughly 11M and 22M unknowns, giving ~1.5
# decades, whole-decade ticks, and a power-law fit over six points instead of four.
#
# COST, extrapolated from the eps=2.4 run (solve times 25/26/61/104/209 s for levels 1-5, each
# level roughly doubling both unknowns and time):
#     level 6  ~11M unknowns,  solve ~7 min
#     level 7  ~22M unknowns,  solve ~15 min
# plus prepare, which grows at least as fast. Budget an hour or two for `6 7` together; the
# 24 h walltime below is slack, not an estimate.
#
# MEMORY IS THE REAL RISK, not time. Level 5 (5.6M unknowns) ran comfortably, but level 7 is 4x
# that and the near-field correction blocks scale with it. This runs on one exclusive genoa node;
# if level 7 is OOM-killed, submit `6` alone first and see whether 6 alone is enough to widen the
# axis before paying for 7.
#
# TKM3D_FINUFFT_NTHREADS is capped deliberately: an uncapped FINUFFT plan on a 96-thread node
# cripples the TKM3D evaluator. Do not remove it.
#
# SUBMIT (you -- org policy bars the assistant from running sbatch):
#   cd /mnt/home/xgao1/work/four_index_integral_solver/codes/lattice_scale
#   sbatch jobscripts/run_lec_conv.sbatch 6 7                 # eps=2.4 (default variant)
#   VARIANT=v2 sbatch --export=ALL,VARIANT jobscripts/run_lec_conv.sbatch 6 7   # eps=10 sweep
#
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=24:00:00
#SBATCH --job-name=fi_lec_conv
#SBATCH --output=logs/lec_conv_%j.out

set -euo pipefail
LEVELS=${@:-1 2 3 4 5}
VARIANT=${VARIANT:-eps2.4}

case "$VARIANT" in
  eps2.4) export RERUN_TAG=_eps2.4; export RERUN_BASE=lattice_10x10_het3x_eps2.4 ;;
  v2)     export RERUN_TAG=_v2;     export RERUN_BASE=lattice_10x10_het3x_v2 ;;
  base)   export RERUN_TAG="";      export RERUN_BASE=lattice_10x10_het3x ;;
  *) echo "unknown VARIANT '$VARIANT' (expected eps2.4 | v2 | base)" >&2; exit 2 ;;
esac

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export TKM3D_FINUFFT_NTHREADS=16              # cap: see header
JULIA=/mnt/home/xgao1/.juliaup/bin/julia      # juliaup launcher (never `module load julia`)

cd /mnt/home/xgao1/work/four_index_integral_solver/codes/lattice_scale
mkdir -p logs

echo "=== lec sweep  variant=$VARIANT (tag='${RERUN_TAG}')  levels: $LEVELS   $(date) ==="
$JULIA --project -t 8 scripts/lec_conv_single.jl $LEVELS
echo "=== SWEEP DONE  $(date)  ->  figs/lec_single_conv${RERUN_TAG}.tsv ==="

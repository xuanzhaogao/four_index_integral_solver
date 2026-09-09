#!/bin/bash
# Re-run ONLY eval + assemble on a campaign whose solve is already done, to benchmark the
# full-matrix evaluation (BI b41ee80) against the symmetry-restricted one it replaces.
#
# WHY. The symmetry restriction cut the potential evaluation 41% but paid for it with a Dict
# remapping every global target index into each batch's subset -- built per batch (up to 7.3M
# entries) and probed once per target per row, in a serial loop. Measured on kb31 batch 1:
# 235 s of evaluation inside a 1896 s t_phi. Evaluating every row restores direct indexing.
#
# PREDICTION for lattice_conv_l3_eps2.4_k46 (57 batches, K = 46, 7,375,550 targets), to be
# checked against the log rather than assumed:
#     potential  57 x (7.50 + 7.30*46)   = 19,574 s      [fit from the job 7010961 K sweep]
#     contraction 2623^2 x 0.393 ms      =  2,705 s      [rate back-derived from the Aug 30 run]
#     total 22,279 s = 6.19 node-hours -> 391 s per batch, 6 rounds on 10 nodes = ~39 min
# against 102.2 min measured with the triangle. Two falsifiable details beyond the total:
#   * per-batch t_phi should now be UNIFORM (~391 s). The old 1902 -> 143 s ramp across batch
#     ids was the shrinking row count; every batch now owns every row.
#   * assemble's max_rel_asym becomes an independent end-to-end check (nothing is mirrored),
#     so it should land near the 1.71e-5 the triangle run reported over its overlap only.
#
# BEFORE SUBMITTING, the previous V files must be out of the way -- eval_batch skips any batch
# whose V file is complete, so a rerun over them is a silent no-op. Already done for k46:
#     V -> V_triangle,  V_full_eV.jls -> V_full_eV_triangle.jls,  report.txt -> report_triangle.txt
# Do the same for any other campaign before pointing this at it.
#
# SUBMIT (you — org policy bars the assistant from running sbatch):
#   cd /mnt/home/xgao1/work/four_index_integral_solver/codes/lattice_scale
#   sbatch --reservation=rocky8 --nodes=10 jobscripts/eval_rerun.sbatch campaigns/lattice_conv_l3_eps2.4_k46.toml
#
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=04:00:00
#SBATCH --job-name=fi_eval_rerun
#SBATCH --output=logs/eval_rerun_%j.out

set -euo pipefail
TOML=${1:?usage: sbatch --nodes=<M> jobscripts/eval_rerun.sbatch <campaign.toml>}

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# Measured irrelevant for this path (job 7012357: 8/96 Julia threads x 1/96 OpenBLAS, all four
# within 3%). Held at 1 to match the run this benchmark is compared against.
export OPENBLAS_NUM_THREADS=1
export JULIA_GLUE_THREADS=8
JULIA=/mnt/home/xgao1/.juliaup/bin/julia

cd /mnt/home/xgao1/work/four_index_integral_solver/codes/lattice_scale
mkdir -p logs

for PHASE in eval assemble; do
    echo "=== $PHASE  $(date) ==="
    $JULIA --project -t $JULIA_GLUE_THREADS driver.jl "$TOML" "$PHASE"
done
echo "=== EVAL RERUN DONE  $(date) ==="

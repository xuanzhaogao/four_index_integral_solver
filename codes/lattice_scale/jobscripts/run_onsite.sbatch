#!/bin/bash
# Onsite-U probe campaign in ONE resumable job: prepare -> solve -> onsite_eval.
# ONSITE-ONLY (not the full V matrix): `solve` produces one σ per site (distributed, one
# SlurmManager worker per node), then `scripts/onsite_eval.jl` evaluates each site's Φ at
# its OWN support and contracts → U_i (in eV) into <root>/onsite_U.tsv. So the cost is
# 1792 solves + 1792 LOCAL evals — we deliberately SKIP consolidate/eval/assemble, which
# would build the full n×n matrix (Φ at the shared union of all supports). RESUMABLE:
# resubmit the same line — solved batches and finished onsite_U rows are skipped.
#
# WORKFLOW (you submit):
#   1. Pilot the cheap 1D line first (one node, ~28 sites):
#        sbatch --nodes=1  jobscripts/run_onsite.sbatch campaigns/onsite_line.toml
#      then: julia --project scripts/plot_onsite.jl campaigns/onsite_line.toml
#      Read the per-site solve/eval cost from the logs before scaling.
#   2. Full 2D map (1792 sites; size M for the SOLVE from the pilot):
#        sbatch --nodes=<M> jobscripts/run_onsite.sbatch campaigns/onsite_grid.toml
#      then: julia --project scripts/plot_onsite.jl campaigns/onsite_grid.toml
#
# NOTE: onsite_eval here runs on the head node (looping batches, each eval using the node's
# cores). For the big grid that leaves the other solve nodes idle during eval — to avoid
# that, run prepare+solve via this script, then onsite_eval.jl as a separate --nodes=1 job.
# To amortize the SOLVES, raise [batching] n_centers_per_batch (multi-RHS over nearby sites).
#
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=12:00:00
#SBATCH --job-name=fi_onsite
#SBATCH --output=logs/onsite_%j.out

set -euo pipefail
TOML=${1:?usage: sbatch --nodes=<M> jobscripts/run_onsite.sbatch <campaign.toml>}

# pin BOTH thread pools (unpinned threads silently corrupt timings — known footgun)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# OpenBLAS gets ONE thread. The 2.95x-at-64-threads figure this once cited comes from the
# Sec. 5.2 benchmark -- a SINGLE process at -t 64/96, where an equally sized BLAS pool
# oversubscribes the node against FMM3D and FINUFFT. It does NOT describe this pipeline: job
# 7012357 timed the campaign's own evaluation at 8 and 96 Julia threads crossed with OPENBLAS
# 1 and 96, and all four cells landed within 3% (231.9-239.6 s). Keeping 1 because it is
# right for the single-process case and measured harmless here.
export OPENBLAS_NUM_THREADS=1
export JULIA_GLUE_THREADS=8
JULIA=/mnt/home/xgao1/.juliaup/bin/julia      # juliaup launcher (Julia 1.12 = release default)

cd /mnt/home/xgao1/work/four_index_integral_solver/codes/lattice_scale
mkdir -p logs

for PHASE in prepare solve; do
    echo "=== $PHASE  $(date) ==="
    $JULIA --project -t $JULIA_GLUE_THREADS driver.jl "$TOML" "$PHASE"
done
echo "=== onsite_eval  $(date) ==="
$JULIA --project -t $JULIA_GLUE_THREADS scripts/onsite_eval.jl "$TOML"
echo "=== ONSITE CAMPAIGN DONE  $(date)  ->  <root>/onsite_U.tsv ==="

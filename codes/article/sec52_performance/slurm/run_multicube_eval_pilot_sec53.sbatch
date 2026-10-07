#!/bin/bash
# PILOT for the Sec. 5.2 evaluation benchmark: cost per source of evaluating K solved sources
# at the FIXED target set of the Sec. 5.3 model -- the union of the quadrature points of all
# 198 orbitals, N_p = 7,168,390 points.
#
# WHY A PILOT AND NOT THE SWEEP. N_p is 53x the row protocol's 136,580 targets and 11x the
# K = 46 batch grid, so the corrected pottrg map over N_p targets is a cost and a memory
# footprint nobody has measured. The row-protocol sweep already peaked at 1223 GB of a
# 1502 GiB node at K = 46. Sizing the full sweep from a guess would risk ten OOMs; this
# measures the eval's time and memory at small K, where the solve leaves the most headroom
# (K = 1 is 43 GB, K = 4 is 119 GB), and the sweep gets sized from the result.
#
# WHAT IT MEASURES. One evaluate_batch_potential call gives Phi[g, b] for all K sources at all
# N_p targets -- the K x N_p potential matrix the ERI assembly consumes, since every
# V[rho_a, rho_b] is a contraction of a column of Phi against a target density rho_a. With N_p
# FIXED, t_eval/K is a clean per-source cost. Two K values are run in ONE process so the
# amortization decomposes: the pottrg map over N_p targets is a single K-independent cost
# shared by all K columns, so
#     t_eval(K) ~ t_pottrg + K * t_field   =>   per-source ~ t_pottrg/K + t_field
# and the K = 1 / K = 4 pair separates t_pottrg from t_field. (The earlier "block" mode used
# the batch's own grid, which grows with K; its flat per-source time was two effects
# cancelling and measured nothing.)
#
# Unlike the small-target-set case, N_p spans 29 x 31 A while each source cluster is compact,
# so most targets are FAR from a given source and the nd=K batched far-field FMM inside
# evaluate_batch_potential should finally engage. That is the specific claim under test.
#
# Step 1 builds the target cache (~200 s, tens of GB, in its own process so the memory is
# released before the benchmark). It is idempotent: it exits immediately if the cache exists.
#
# Node image pinned to rocky8 to match every other sec53 measurement. SUBMIT-TIME TRAP: a
# stale SBATCH_RESERVATION=rocky9 in the environment overrides script directives (sbatch(1):
# "environment variables will override any options set in a batch script"), so pass the
# reservation on the command line:
#
# Submit (you):  sbatch --reservation=rocky8 sec52_performance/slurm/run_multicube_eval_pilot_sec53.sbatch
#
#SBATCH --job-name=mc53_evalpilot
#SBATCH --partition=ccm
#SBATCH --constraint=genoa,rocky8
#SBATCH --reservation=rocky8
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/article/sec52_performance/logs/slurm_mc53_evalpilot_%j.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/article/sec52_performance/logs/slurm_mc53_evalpilot_%j.err

set -uo pipefail
NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/article
J=/mnt/home/xgao1/.juliaup/bin/julia    # juliaup binary (never module load julia)
source /mnt/home/xgao1/work/four_index_integral_solver/codes/article/common/_provenance.sh

export RERUN_TAG=_sec53evalpilot
export MULTICUBE_GEOM=sec53
export PRECONDITION=1
export CORRECT_EDGES=1
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96               # TKM3D caps FINUFFT at 16 internally
provenance

echo "==== step 1: build the fixed 198-orbital target cache  $(date '+%F %T') ===="
$J --project=$NR --startup-file=no $NR/sec52_performance/scripts/prep_model_targets.jl || exit 1

echo "==== step 2: eval pilot at K = 1 and K = 4  $(date '+%F %T') ===="
export MULTICUBE_EVAL=model
export MULTICUBE_CUTOFFS=0.0,1.5        # K = 1 and K = 4, both in ONE process
$J --project=$NR --startup-file=no $NR/sec52_performance/scripts/run_multicube.jl

echo "==== DONE  $(date '+%F %T') ===="

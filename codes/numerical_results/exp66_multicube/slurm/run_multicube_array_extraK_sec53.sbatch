#!/bin/bash
# Sec. 5.2 K-sweep EXTENSION on the Sec. 5.3 system: K = 37, 40, 46, appended to the
# K = 1..31 points from run_multicube_array_sec53.sbatch. Same node type, same image, same
# discretization, same env -- only the cutoff differs, so panel (b) of the scaling figure can
# extend without a break in conditions.
#
# WHY these three and no more. Measured RSS over the K = 1..31 sweep is dead linear at
# 26.2 GB per RHS (RSS ~= 16.6 + 26.2 K GB; predicts 828.8 at K = 31 vs 829.5 measured), and
# the graphene neighbour shells jump 46 -> 58 with nothing between:
#
#     cutoff  5.0   5.5   6.0   6.5   6.8
#     K        31    37    40    46    58
#     RSS GB  830   986  1065  1222  1536   <- 1536 vs a 1538 GB node: certain OOM
#
# 1538 GB is the largest memory tier in ccm (genoa), so K = 46 (~1222 GB, 316 GB headroom) is
# the last point reachable on hardware comparable to the rest of the curve. Do NOT add 6.8.
#
# WHAT THIS TESTS. The per-RHS time is expected to stay flat, not turn up: the dominant work
# is niter*K FMM matvecs, strictly linear in K. Block-Arnoldi orthogonalization is the only
# K^2 term and it is negligible here -- ~2 N K^2 m^2 = 0.5 TFLOP at K = 31 and 1.2 TFLOP at
# K = 46, i.e. sub-second of BLAS-3 against a ~900-1800 s block solve. So the binding limit on
# batch size is memory CAPACITY, not a compute crossover; if per-RHS does rise it should do so
# only as the live heap approaches node capacity (GC/allocation pressure).
#
# Writes exp66_multicube/data_sec53/raw/ alongside the existing records.
#
# Node image pinned to rocky8 to match the K = 1..31 points (see the main sec53 script).
# SUBMIT-TIME TRAP: a stale SBATCH_RESERVATION=rocky9 in the environment overrides script
# directives (sbatch(1): "environment variables will override any options set in a batch
# script"), so pass the reservation on the command line:
#
# Submit (you):  sbatch --reservation=rocky8 exp66_multicube/slurm/run_multicube_array_extraK_sec53.sbatch
# After it finishes:  RERUN_TAG=_sec53 julia --project=. plot_scripts/plot_scaling.jl
#
#SBATCH --job-name=mc53_extraK
#SBATCH --partition=ccm
#SBATCH --constraint=genoa,rocky8
#SBATCH --reservation=rocky8
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --array=0-2
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp66_multicube/logs/slurm_mc53_extraK_%A_%a.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp66_multicube/logs/slurm_mc53_extraK_%A_%a.err

set -uo pipefail
CUTOFFS=(5.5 6.0 6.5)                   # array index -> neighbor cutoff (K = 37, 40, 46)
export MULTICUBE_CUTOFF=${CUTOFFS[$SLURM_ARRAY_TASK_ID]}

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
J=/mnt/home/xgao1/.juliaup/bin/julia    # juliaup binary (never module load julia)
source /mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/jobscripts/_provenance.sh

export RERUN_TAG=_sec53
export MULTICUBE_GEOM=sec53
export PRECONDITION=1
export CORRECT_EDGES=1
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96               # TKM3D caps FINUFFT at 16 internally
provenance

echo "=== array task $SLURM_ARRAY_TASK_ID  cutoff=$MULTICUBE_CUTOFF  $(date '+%F %T') ==="
$J --project=$NR --startup-file=no $NR/exp66_multicube/scripts/run_multicube.jl
echo "=== DONE cutoff=$MULTICUBE_CUTOFF  $(date '+%F %T') ==="

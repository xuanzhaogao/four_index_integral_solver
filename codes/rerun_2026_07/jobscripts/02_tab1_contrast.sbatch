#!/bin/bash
# RERUN 02 -- Table 1 (dielectric-contrast study, article Sec. 6.1).
#
# Why: the N_iter column is the point of this table, and the design note predicts
# the story changes from "roughly doubles toward the conductor limit" (2.4x) to
# "grows ~40%". The E_r column should be unchanged.
#
# Writes exp63_contrast/data_v2/.
#
# Submit (you):  sbatch 02_tab1_contrast.sbatch
#
#SBATCH --job-name=rr02_contrast
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=16:00:00
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/logs/rr02_contrast_%j.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/logs/rr02_contrast_%j.err

set -uo pipefail
NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
J=/mnt/home/xgao1/.juliaup/bin/julia
source /mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/jobscripts/_provenance.sh

export RERUN_TAG=_v2
export PRECONDITION=1
export GMRES_VERBOSE=0
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96
provenance

echo "==== contrast sweep  $(date '+%F %T') ===="
$J --project=$NR --startup-file=no $NR/exp63_contrast/scripts/run_contrast_ratio.jl
echo "==== analyze  $(date '+%F %T') ===="
$J --project=$NR --startup-file=no $NR/exp63_contrast/scripts/analyze_contrast.jl
echo "==== DONE  $(date '+%F %T') ===="

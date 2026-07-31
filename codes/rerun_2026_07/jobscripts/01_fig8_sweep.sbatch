#!/bin/bash
# RERUN 01 -- Fig. 8 (three-material benchmark, article Sec. 6.1).
#
# Why: the published run had the edge correction ON but no preconditioner. The
# preconditioner does not move sigma, so panel (b)'s V-error curves should be
# reproduced; panel (c)'s GMRES counts are what change.
#
# Writes exp61_convergence/data_v2/ -- the published data/ is left untouched.
#
# Submit (you):  sbatch 01_fig8_sweep.sbatch
#
#SBATCH --job-name=rr01_fig8
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=12:00:00
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/logs/rr01_fig8_%j.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/logs/rr01_fig8_%j.err

set -uo pipefail
NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
J=/mnt/home/xgao1/.juliaup/bin/julia
source /mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/jobscripts/_provenance.sh

export RERUN_TAG=_v2
export PRECONDITION=1
export GMRES_VERBOSE=0
export JULIA_NUM_THREADS=96
provenance

# Main p x r sweep + the p=8 reference. Edge correction on (Harness default),
# preconditioner on.
echo "==== fig1 sweep (edges on, precond on)  $(date '+%F %T') ===="
export OMP_NUM_THREADS=96
$J --project=$NR --startup-file=no $NR/exp61_convergence/scripts/run_fig1_sweep.jl

# No-edge ablation: the dashed curves. OMP capped at 16 to dodge the high-thread
# FFTW-plan pathology, exactly as the published run did.
echo "==== fig1 no-edge ablation (edges off, precond on)  $(date '+%F %T') ===="
export OMP_NUM_THREADS=16
$J --project=$NR --startup-file=no $NR/exp61_convergence/scripts/run_fig1_noedges.jl

echo "==== analyze  $(date '+%F %T') ===="
export OMP_NUM_THREADS=96
$J --project=$NR --startup-file=no $NR/exp61_convergence/scripts/analyze_fig1.jl
echo "==== DONE  $(date '+%F %T') ===="

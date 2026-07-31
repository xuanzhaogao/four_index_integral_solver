#!/bin/bash
# RERUN 07 -- single-orbital l_ec convergence (article Sec. 6.4: Fig. 10 panel (c)).
#
# Panel (c) is the paper's ONLY accuracy claim for the tensor, and it is the one
# most compromised by the missing edge correction: it is a self-convergence
# study, so all five levels shared the same uncorrected operator and the same
# bias. The error fell along a clean line and Richardson extrapolated smoothly --
# to the wrong limit. Rerunning it is the only way to know what the real
# convergence rate and U_inf are.
#
# Runs levels 1..5 sequentially on one node (prepare + single-batch solve per
# level). Each level writes its own small ceph campaign root; RERUN_TAG suffixes
# both the campaign name and the output TSV, so the published
# figs/lec_single_conv.tsv is untouched.
#
# COST: published solve times were 78, 126, 317, 626, 1242 s for levels 1-5,
# plus prepare and eval per level. Level 5 is 5.56M unknowns.
#
# Submit (you):  sbatch 07_lec_single.sbatch
#
#SBATCH --job-name=rr07_lec
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=24:00:00
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/logs/rr07_lec_%j.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/logs/rr07_lec_%j.err

set -euo pipefail
LS=/mnt/home/xgao1/work/four_index_integral_solver/codes/lattice_scale
J=/mnt/home/xgao1/.juliaup/bin/julia
source /mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/jobscripts/_provenance.sh

export RERUN_TAG=_v2
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export JULIA_GLUE_THREADS=8
export TKM3D_FINUFFT_NTHREADS=16
provenance

cd $LS
echo "==== lec_conv_single levels 1-5  ($(date '+%F %T')) ===="
$J --project -t $JULIA_GLUE_THREADS scripts/lec_conv_single.jl 1 2 3 4 5
echo "==== DONE  ($(date '+%F %T')) ===="

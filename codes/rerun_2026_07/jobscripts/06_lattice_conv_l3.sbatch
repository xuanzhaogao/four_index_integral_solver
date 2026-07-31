#!/bin/bash
# RERUN 06 -- lattice_conv_l3 (article Sec. 6.4: Fig. 10 panels (a) and (b)).
#
#
# This is the finest campaign (edge_refine_level 3, l_ec = 1.14 A, all solver
# tolerances 1e-5) and the source of the tensor plotted in Fig. 10(a),(b) and of
# every on-site U quoted in Sec. 6.4. The
# published run had correct_edges = false AND precondition = false, because
# solve_dielectric_lattice_batch defaulted both to false and the campaign path
# never overrode them. BI now defaults both to true, and solve_batch_core still
# passes neither -- so this rerun picks them up with no TOML change.
#
# PREREQUISITE: generate the retargeted campaign first, or this will write into
# the published ceph root and destroy the submitted data:
#   julia --project=codes/lattice_scale codes/rerun_2026_07/scripts/gen_rerun_campaigns.jl
#
# RESUMABLE: phase status is file-derived; resubmit the identical line to continue.
#
# COST: the published run was 8 h 20 m on 10 nodes (83 node-hours) at 67-71
# GMRES iterations and 1.35M dof per batch -- the most expensive job here. With
# both fixes the iteration count should fall sharply, so expect LESS wall-clock
# despite the added near-pair assembly. Walltime is set generously anyway.
#
# Submit (you):  sbatch --nodes=10 06_lattice_conv_l3.sbatch
#
#SBATCH --job-name=rr06_conv_l3
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=48:00:00
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/logs/rr06_conv_l3_%j.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/logs/rr06_conv_l3_%j.err

set -euo pipefail
TOML=campaigns/lattice_conv_l3_v2.toml

LS=/mnt/home/xgao1/work/four_index_integral_solver/codes/lattice_scale
J=/mnt/home/xgao1/.juliaup/bin/julia
source /mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/jobscripts/_provenance.sh

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export JULIA_GLUE_THREADS=8
export TKM3D_FINUFFT_NTHREADS=16
provenance

cd $LS
test -f "$TOML" || { echo "MISSING $TOML -- run gen_rerun_campaigns.jl first"; exit 1; }
mkdir -p logs   # driver.jl writes per-phase logs relative to cwd

run_phase () { echo "==== $1  ($(date '+%F %T')) ===="; $J --project $2 driver.jl "$TOML" "$1"; }

echo "==================== CAMPAIGN: $TOML  ($(date '+%F %T')) ===================="
run_phase prepare     ""
run_phase solve       "-t $JULIA_GLUE_THREADS"
run_phase consolidate ""
run_phase eval        "-t $JULIA_GLUE_THREADS"
run_phase assemble    ""
echo "==================== DONE: $TOML  ($(date '+%F %T')) ===================="
$J --project driver.jl "$TOML" status

#!/bin/bash
# RERUN 05 -- lattice_10x10_het3x (article Sec. 6.4: Table 4, and the geometry
# behind Fig. 9).
#
# This is the campaign whose per-stage banner timestamps ARE Table 4. The
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
# COST: the published run was 3 h 51 m on 10 nodes (38.5 node-hours) at 61-86
# GMRES iterations. With both fixes the iteration count should fall sharply, so
# expect LESS wall-clock despite the added near-pair assembly. Walltime is set
# generously anyway.
#
# Submit (you):  sbatch --nodes=10 05_lattice_het3x.sbatch
#
#SBATCH --job-name=rr05_het3x
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=48:00:00
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/logs/rr05_het3x_%j.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/rerun_2026_07/logs/rr05_het3x_%j.err

set -euo pipefail
TOML=campaigns/lattice_10x10_het3x_v2.toml

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

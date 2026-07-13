#!/bin/bash
# Full four-index V_ijkl benchmark on the ×3 heterojunction (campaigns/lattice_10x10_het3x.toml):
# real graphene honeycomb (198 orbitals) in a square window centered on the Si|SiO₂ junction,
# on the converged 270³ substrate cubes + tight ε=10 slab. Runs the whole pipeline in ONE
# allocation: prepare → solve → consolidate → eval → assemble → V_full_eV.jls + report.txt.
#
# The two heavy phases (solve, eval) spawn ONE SlurmManager worker per node (--ntasks-per-node=1,
# 96-core task/node), so M nodes process batches M-way parallel; the three light phases run on the
# head. RESUMABLE — status is file-derived on ceph, so resubmit the SAME line after a walltime
# kill and finished batches are skipped; consolidate/assemble are idempotent and guarded.
#
# SIZING (measured pilot, batch 1 on one genoa node — see logs/bench_pilot.log):
#   n_batches = <NB>,  t_solve/batch ≈ <TS>,  dof ≈ 801k,  batch file ≈ <SZ>.
#   With M nodes, solve ≈ ceil(NB/M)·t_solve. Pick M so solve+eval fit --time (start ~10 nodes).
#
# SUBMIT (you — org policy bars the assistant from running sbatch):
#   cd /mnt/home/xgao1/work/four_index_integral_solver/codes/lattice_scale
#   sbatch --nodes=<M> jobscripts/run_benchmark_het3x.sbatch
# Optionally pass a different campaign TOML as $1.
#
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=24:00:00
#SBATCH --job-name=fi_bench_het3x
#SBATCH --output=logs/bench_het3x_%j.out

set -euo pipefail
TOML=${1:-campaigns/lattice_10x10_het3x.toml}

# pin BOTH thread pools to the whole node (unpinned threads silently corrupt timings)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export JULIA_GLUE_THREADS=8
# cap TKM3D's internal FINUFFT off the 96-way FFTW-plan pathology (eval uses TKM3D; FMM3D keeps 96).
# TKM3D already defaults this to 16, but set it explicitly for safety.
export TKM3D_FINUFFT_NTHREADS=16
JULIA=/mnt/home/xgao1/.juliaup/bin/julia      # juliaup launcher (Julia 1.12 = release default)

cd /mnt/home/xgao1/work/four_index_integral_solver/codes/lattice_scale
mkdir -p logs

run_phase () {   # $1 = phase, $2 = extra julia flags ("-t 8" for the parallel phases, "" for head)
    echo "==================== PHASE: $1  ($(date '+%F %T')) ===================="
    $JULIA --project $2 driver.jl "$TOML" "$1"
    echo "-------------------- $1 done ($(date '+%F %T')) --------------------"
}

run_phase prepare       ""
run_phase solve         "-t $JULIA_GLUE_THREADS"
run_phase consolidate   ""
run_phase eval          "-t $JULIA_GLUE_THREADS"
run_phase assemble      ""

echo "==================== BENCHMARK COMPLETE ($(date '+%F %T')) ===================="
$JULIA --project driver.jl "$TOML" status

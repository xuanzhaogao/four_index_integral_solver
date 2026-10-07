#!/bin/bash
# l_ec self-convergence: run the FULL four-index pipeline (prepare -> solve -> consolidate ->
# eval -> assemble) for ONE convergence campaign passed as $1. Submit it once per l_ec so the two
# runs are independent allocations (run in parallel, each resumable on its own):
#   lattice_conv_l2.toml  (l_ec = 2.27 bohr, edge_refine_level 2)
#   lattice_conv_l3.toml  (l_ec = 1.14 bohr, edge_refine_level 3)
# Both: identical geometry / 198 orbitals / all solver tolerances 1e-5; only l_ec differs. Each
# writes <root>/V_full_eV.jls + report.txt; compare with scripts/compare_conv.jl.
#
# COST WARNING: rhs_tol tightened 1e-3 -> 1e-5, so interfaces are MUCH larger than the main
# benchmark (651k DOF / 3.2 h solve on 10 nodes at 1e-3) — expect longer per-batch solves,
# l_ec=1.14 the heavier. RESUMABLE: status is file-derived; resubmit the SAME line to continue.
#
# SUBMIT (you), two independent jobs:
#   cd .../codes/article/sec53_lattice
#   sbatch --nodes=<M> jobscripts/run_benchmark_conv.sbatch campaigns/lattice_conv_l2.toml
#   sbatch --nodes=<M> jobscripts/run_benchmark_conv.sbatch campaigns/lattice_conv_l3.toml
#
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=48:00:00
#SBATCH --job-name=fi_conv
#SBATCH --output=logs/conv_%j.out

set -euo pipefail
TOML=${1:?usage: sbatch --nodes=<M> jobscripts/run_benchmark_conv.sbatch <campaign.toml>}

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# OpenBLAS gets ONE thread. The 2.95x-at-64-threads figure this once cited comes from the
# Sec. 5.2 benchmark -- a SINGLE process at -t 64/96, where an equally sized BLAS pool
# oversubscribes the node against FMM3D and FINUFFT. It does NOT describe this pipeline: job
# 7012357 timed the campaign's own evaluation at 8 and 96 Julia threads crossed with OPENBLAS
# 1 and 96, and all four cells landed within 3% (231.9-239.6 s). Keeping 1 because it is
# right for the single-process case and measured harmless here.
export OPENBLAS_NUM_THREADS=1
export JULIA_GLUE_THREADS=8
export TKM3D_FINUFFT_NTHREADS=16
JULIA=/mnt/home/xgao1/.juliaup/bin/julia

cd /mnt/home/xgao1/work/four_index_integral_solver/codes/article/sec53_lattice
mkdir -p logs

run_phase () { echo "==== $1  ($(date '+%F %T')) ===="; $JULIA --project $2 driver.jl "$TOML" "$1"; }

echo "==================== CAMPAIGN: $TOML  ($(date '+%F %T')) ===================="
run_phase prepare     ""
run_phase solve       "-t $JULIA_GLUE_THREADS"
run_phase consolidate ""
run_phase eval        "-t $JULIA_GLUE_THREADS"
run_phase assemble    ""
echo "==================== DONE: $TOML  ($(date '+%F %T')) ===================="
$JULIA --project driver.jl "$TOML" status

#!/bin/bash
# Level-3 production campaign at the cRPA slab permittivity (campaigns/lattice_conv_l3_eps2.4.toml):
# the 198-orbital ×3 heterojunction lattice with the graphene slab at ε = 2.4 instead of 10.
#
# WHY THIS RUN EXISTS. The completed ε=2.4 run (lattice_10x10_het3x_eps2.4, job 6962885) is NOT a
# like-for-like replacement for the published Figure 9 data (lattice_conv_l3_v2). It differs in two
# ways that both loosen it:
#     edge_refine_level  2 (l_ec = 2.27 Å, N ≈ 0.65M)   vs  3 (l_ec = 1.14 Å, N ≈ 1.35M)
#     rhs_tol            1e-3                            vs  1e-5
# This campaign matches the published one on both, changing only the slab permittivity.
#
# Same five-phase pipeline as jobscripts/run_all.sbatch, with the campaign fixed. RESUMABLE:
# status is derived from files on ceph, so resubmitting the same line skips finished batches and
# consolidate/assemble are idempotent.
#
# SIZING. The level-2 ε=2.4 run took 1h32m on 10 nodes (198 batches). Level 3 roughly doubles the
# interface unknowns, and the tighter rhs_tol adds further cost, so budget several times that;
# --time below is deliberately generous since an overrun is recovered by one more submit.
#
# SUBMIT (you — org policy bars the assistant from running sbatch):
#   cd /mnt/home/xgao1/work/four_index_integral_solver/codes/lattice_scale
#   sbatch --nodes=10 jobscripts/run_conv_l3_eps2.4.sbatch
#
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=24:00:00
#SBATCH --job-name=fi_convl3_eps24
#SBATCH --output=logs/convl3_eps2.4_%j.out

set -euo pipefail
TOML=${1:-campaigns/lattice_conv_l3_eps2.4.toml}

# pin BOTH thread pools (unpinned threads silently corrupt timings — known footgun)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export JULIA_GLUE_THREADS=8
JULIA=/mnt/home/xgao1/.juliaup/bin/julia      # juliaup launcher (never `module load julia`)

cd /mnt/home/xgao1/work/four_index_integral_solver/codes/lattice_scale
mkdir -p logs

for PHASE in prepare solve consolidate eval assemble; do
    echo "=== $PHASE  $(date) ==="
    $JULIA --project -t $JULIA_GLUE_THREADS driver.jl "$TOML" "$PHASE"
done
echo "=== CAMPAIGN DONE  $(date)  ->  /mnt/ceph/users/xgao1/four_index/lattice_conv_l3_eps2.4/V_full_eV.jls ==="

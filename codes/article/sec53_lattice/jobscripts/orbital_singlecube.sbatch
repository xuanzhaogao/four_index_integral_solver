#!/bin/bash
# Orbital "1 slab on 1 cube" onsite-U reference, as a Slurm ARRAY of 6 independent cases:
#   scale ∈ {2, 3, 4}  ×  cube ε ∈ {3.9 (SiO2), 11.9 (Si)}.
# Each array task takes ONE exclusive genoa node and uses all 96 cores. Real graphene orbital
# source; observable = full screened onsite U in eV (scripts/orbital_singlecube.jl).
#
# SUBMIT (you — I don't run sbatch):
#   cd /mnt/home/xgao1/work/four_index_integral_solver/codes/article/sec53_lattice
#   sbatch jobscripts/orbital_singlecube.sbatch
# Each case writes figs/orbital_singlecube_eps<eps>_s<scale>.tsv and its own logs/orbsc_<A>_<a>.out.
# Resubmitting re-runs all tasks (solve_batch skips if its BatchResult already exists on ceph).
#
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=96
#SBATCH --time=12:00:00
#SBATCH --array=0-5
#SBATCH --job-name=orb_sc
#SBATCH --output=logs/orbsc_%A_%a.out

set -euo pipefail

# array index → (scale, eps): 0:(2,3.9) 1:(2,11.9) 2:(3,3.9) 3:(3,11.9) 4:(4,3.9) 5:(4,11.9)
scales=(2.0 2.0 3.0 3.0 4.0 4.0)
epses=( 3.9 11.9 3.9 11.9 3.9 11.9)
SCALE=${scales[$SLURM_ARRAY_TASK_ID]}
EPS=${epses[$SLURM_ARRAY_TASK_ID]}

# pin BOTH thread pools to the full node (known footgun: unpinned threads corrupt timings)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
# OpenBLAS gets ONE thread. The 2.95x-at-64-threads figure this once cited comes from the
# Sec. 5.2 benchmark -- a SINGLE process at -t 64/96, where an equally sized BLAS pool
# oversubscribes the node against FMM3D and FINUFFT. It does NOT describe this pipeline: job
# 7012357 timed the campaign's own evaluation at 8 and 96 Julia threads crossed with OPENBLAS
# 1 and 96, and all four cells landed within 3% (231.9-239.6 s). Keeping 1 because it is
# right for the single-process case and measured harmless here.
export OPENBLAS_NUM_THREADS=1
# cap TKM3D's internal FINUFFT away from the 96-way FFTW plan pathology (FMM3D still uses 96)
export TKM3D_FINUFFT_NTHREADS=16
JULIA=/mnt/home/xgao1/.juliaup/bin/julia

cd /mnt/home/xgao1/work/four_index_integral_solver/codes/article/sec53_lattice
mkdir -p logs

echo "=== orbital singlecube  scale=$SCALE  eps=$EPS  $(date) ==="
$JULIA --project -t "$SLURM_CPUS_PER_TASK" scripts/orbital_singlecube.jl "$SCALE" "$EPS"
echo "=== DONE  scale=$SCALE  eps=$EPS  $(date) ==="

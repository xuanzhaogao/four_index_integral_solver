#!/bin/bash
# Phase 2: one array task per OMP_NUM_THREADS, each on its OWN exclusive genoa node, timing the
# nd=64 FMM at N=461,376 (loads the prebuilt operator). Run AFTER the build:
#   JID=$(sbatch --parsable jobscripts/slurm_bench_nd64_build.slurm)
#   sbatch --dependency=afterok:$JID jobscripts/slurm_bench_nd64_array.slurm
# Each task prints a "ROW,omp,3,N,64,t_fmm,t_near" line to append to data/fmm_bench.csv.
#SBATCH --job-name=nd64_sweep
#SBATCH --array=0-7
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=00:40:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --partition=ccm
#SBATCH --constraint=genoa            # 96-core nodes (matches the OMP sweep up to 96)

cd "${SLURM_SUBMIT_DIR:-/mnt/home/xgao1/work/four_index_integral_solver/codes/multi_rhs}"
JULIA=$(ls -d ~/.julia/juliaup/julia-1.12*/bin/julia 2>/dev/null | sort -V | tail -1)   # manifest resolved with 1.12; 1.10 fails
[ -x "$JULIA" ] || JULIA=~/.juliaup/bin/julia
OMPS=(1 2 4 8 16 32 64 96)
t=${OMPS[$SLURM_ARRAY_TASK_ID]}
echo "host=$(hostname)  array_task=$SLURM_ARRAY_TASK_ID  OMP_NUM_THREADS=$t"
OMP_NUM_THREADS=$t OPENBLAS_NUM_THREADS=$t "$JULIA" --project bench_nd64.jl

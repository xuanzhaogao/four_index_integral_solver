#!/bin/bash
# Uniform-point FMM baseline as a job array: one task per OMP_NUM_THREADS, each on its own node
# (--exclusive, isolated timings). No build phase needed — each task regenerates the same seeded
# uniform points. Submit:  sbatch jobscripts/slurm_bench_fmm_uniform_array.slurm
#SBATCH --job-name=fmm_uni_sweep
#SBATCH --array=0-7
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --partition=ccm
#SBATCH --constraint=genoa            # 96-core nodes (matches the OMP sweep up to 96)

cd "${SLURM_SUBMIT_DIR:-/mnt/home/xgao1/work/four_index_integral_solver/codes/multi_rhs}"
JULIA=$(ls -d ~/.julia/juliaup/julia-1.12*/bin/julia 2>/dev/null | sort -V | tail -1)   # manifest resolved with 1.12; 1.10 fails
[ -x "$JULIA" ] || JULIA=~/.juliaup/bin/julia
OMPS=(1 2 4 8 16 32 64 96)
t=${OMPS[$SLURM_ARRAY_TASK_ID]}
echo "host=$(hostname)  array_task=$SLURM_ARRAY_TASK_ID  OMP_NUM_THREADS=$t"
OMP_NUM_THREADS=$t OPENBLAS_NUM_THREADS=$t "$JULIA" --project bench_fmm_uniform.jl

#!/bin/bash
# FMM3D throughput baseline on 3D-uniform points (same N/nd/threads as the slab benchmark).
# Submit:  sbatch jobscripts/slurm_bench_fmm_uniform.slurm
#SBATCH --job-name=fmm_uniform
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=01:30:00
#SBATCH --output=%x_%j.out
#SBATCH --partition=ccm
#SBATCH --constraint=genoa            # 96-core nodes (matches the OMP sweep up to 96)

cd "${SLURM_SUBMIT_DIR:-/mnt/home/xgao1/work/four_index_integral_solver/codes/multi_rhs}"
JULIA=$(ls -d ~/.julia/juliaup/julia-1.12*/bin/julia 2>/dev/null | sort -V | tail -1)   # manifest resolved with 1.12; 1.10 fails
[ -x "$JULIA" ] || JULIA=~/.juliaup/bin/julia
echo "host=$(hostname)  julia=$JULIA  cores=$(nproc)"

for t in 1 2 4 8 16 32 64 96; do
    echo "================ OMP_NUM_THREADS=$t ================"
    OMP_NUM_THREADS=$t OPENBLAS_NUM_THREADS=$t "$JULIA" --project bench_fmm_uniform.jl
done

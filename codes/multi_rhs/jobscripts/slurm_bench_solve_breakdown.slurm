#!/bin/bash
# Block-solve cost breakdown (FMM / near-correction matmul / glue / GMRES) at pinned threads,
# 2-shell graphene K=36, fmm_tol=1e-4.
# Submit:  sbatch jobscripts/slurm_bench_solve_breakdown.slurm
#SBATCH --job-name=solve_breakdown
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --partition=ccm
#SBATCH --constraint=genoa            # 96-core nodes (matches the OMP sweep up to 96)

cd "${SLURM_SUBMIT_DIR:-/mnt/home/xgao1/work/four_index_integral_solver/codes/multi_rhs}"
JULIA=$(ls -d ~/.julia/juliaup/julia-1.12*/bin/julia 2>/dev/null | sort -V | tail -1)   # manifest resolved with 1.12; 1.10 fails
[ -x "$JULIA" ] || JULIA=~/.juliaup/bin/julia
echo "host=$(hostname)  julia=$JULIA  cores=$(nproc)"

rm -f bench_solve_ops.jls
"$JULIA" --project bench_solve_breakdown.jl --build      # build + serialize once
for t in 1 16 96; do
    echo "================ OMP_NUM_THREADS=$t ================"
    OMP_NUM_THREADS=$t OPENBLAS_NUM_THREADS=$t "$JULIA" --project bench_solve_breakdown.jl
done

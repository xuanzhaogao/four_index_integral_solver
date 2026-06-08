#!/bin/bash
# Block-solve cost breakdown (FMM / near-correction matmul / glue / GMRES) at pinned threads.
# Run on worker7010:  bash jobbench_solve_breakdown.sh
set -e
cd "$(dirname "$0")/.."
JULIA=~/.juliaup/bin/julia
rm -f bench_solve_ops.jls
$JULIA --project bench_solve_breakdown.jl --build   # build once at default threads
for t in 1 16 96; do
    echo "================ OMP_NUM_THREADS=$t ================"
    OMP_NUM_THREADS=$t OPENBLAS_NUM_THREADS=$t $JULIA --project bench_solve_breakdown.jl
done

#!/bin/bash
# One-off driver: build the solve-breakdown ops, then run the instrumented block_gmres
# breakdown at pinned OMP_NUM_THREADS = 1 and 16. (Reconciles per_rhs 18x vs breakdown ~4x.)
set -e
cd /mnt/home/xgao1/work/four_index_integral_solver/codes/multi_rhs
JULIA=$(ls -d ~/.julia/juliaup/julia-1.12*/bin/julia 2>/dev/null | sort -V | tail -1)
echo "julia=$JULIA  host=$(hostname)"
rm -f bench_solve_ops.jls
"$JULIA" --project bench_solve_breakdown.jl --build
for t in 1 16; do
    echo "======== OMP_NUM_THREADS=$t ========"
    OMP_NUM_THREADS=$t OPENBLAS_NUM_THREADS=$t "$JULIA" --project bench_solve_breakdown.jl
done
echo ALLDONE

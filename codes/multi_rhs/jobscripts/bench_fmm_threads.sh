#!/bin/bash
# Sweep OMP_NUM_THREADS for the FMM/near-correction micro-benchmark.
# Run on worker7096:  bash jobbench_fmm_threads.sh
set -e
cd "$(dirname "$0")/.."
JULIA=~/.juliaup/bin/julia
rm -f bench_fmm_ops.jls          # rebuild operators fresh (first run builds + serializes)
for t in 1 2 4 8 16 32 64 96; do
    echo "================ OMP_NUM_THREADS=$t ================"
    OMP_NUM_THREADS=$t OPENBLAS_NUM_THREADS=$t $JULIA --project bench_fmm.jl
done

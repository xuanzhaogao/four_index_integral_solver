#!/bin/bash
# Sweep OMP_NUM_THREADS for the uniform-point FMM throughput baseline.
# Run on worker7096:  bash jobbench_fmm_uniform_threads.sh
set -e
cd "$(dirname "$0")/.."
# Use the direct julia version binary (not the ~/.juliaup/bin/julia launcher) to avoid juliaup
# config-lock contention when another juliaup-launched run is active on a machine sharing
# ~/.juliaup over GPFS.
JULIA=$(ls -d ~/.julia/juliaup/julia-*/bin/julia 2>/dev/null | head -1)
[ -x "$JULIA" ] || JULIA=~/.juliaup/bin/julia
for t in 1 2 4 8 16 32 64 96; do
    echo "================ OMP_NUM_THREADS=$t ================"
    OMP_NUM_THREADS=$t OPENBLAS_NUM_THREADS=$t $JULIA --project bench_fmm_uniform.jl
done

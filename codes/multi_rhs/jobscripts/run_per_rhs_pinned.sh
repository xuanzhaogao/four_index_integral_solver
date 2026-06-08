#!/bin/bash
# Re-run the per-RHS K-sweep with PINNED threads (OMP+OpenBLAS), so the K points are all at
# the SAME thread count — fixes the earlier unpinned run where small-K ran ~single-threaded
# and large-K silently multithreaded (the spurious "18x").
set -e
cd /mnt/home/xgao1/work/four_index_integral_solver/codes/multi_rhs
JULIA=$(ls -d ~/.julia/juliaup/julia-1.12*/bin/julia 2>/dev/null | sort -V | tail -1)
echo "julia=$JULIA  host=$(hostname)  OMP=${OMP_NUM_THREADS}"
OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 "$JULIA" --project bench_per_rhs.jl
echo ALLDONE

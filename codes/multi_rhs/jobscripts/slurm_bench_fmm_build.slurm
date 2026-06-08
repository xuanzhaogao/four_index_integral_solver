#!/bin/bash
# Phase 1: build the slab FMM operators (levels 0-3) ONCE and serialize, so the parallel
# job-array tasks only load + time (no build race). Submit then chain the array (see below).
#   JID=$(sbatch --parsable jobscripts/slurm_bench_fmm_build.slurm)
#   sbatch --dependency=afterok:$JID jobscripts/slurm_bench_fmm_array.slurm
#SBATCH --job-name=fmm_build
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=00:40:00
#SBATCH --output=%x_%j.out
#SBATCH --partition=ccm
#SBATCH --constraint=genoa            # 96-core nodes (matches the OMP sweep up to 96)

cd "${SLURM_SUBMIT_DIR:-/mnt/home/xgao1/work/four_index_integral_solver/codes/multi_rhs}"
JULIA=$(ls -d ~/.julia/juliaup/julia-1.12*/bin/julia 2>/dev/null | sort -V | tail -1)   # manifest resolved with 1.12; 1.10 fails
[ -x "$JULIA" ] || JULIA=~/.juliaup/bin/julia
echo "host=$(hostname)  julia=$JULIA"
# Precompile ONCE here, alone, so the parallel array tasks only LOAD a settled cache.
# (The prior run failed because an independently-submitted array raced this precompile over
#  GPFS and corrupted ~/.julia/compiled -> "UndefVarError: StaticData not defined".)
"$JULIA" --project -e 'using Pkg; Pkg.precompile()'
rm -f bench_fmm_ops.jls
"$JULIA" --project bench_fmm.jl --build

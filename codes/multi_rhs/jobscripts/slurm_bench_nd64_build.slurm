#!/bin/bash
# Phase 1: build the level-3 slab operator ONCE (N=461,376) and serialize, so the parallel
# nd=64 array tasks only load + time (no build/precompile race over GPFS). Submit then chain:
#   JID=$(sbatch --parsable jobscripts/slurm_bench_nd64_build.slurm)
#   sbatch --dependency=afterok:$JID jobscripts/slurm_bench_nd64_array.slurm
#SBATCH --job-name=nd64_build
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.out
#SBATCH --partition=ccm
#SBATCH --constraint=genoa            # 96-core nodes (matches the OMP sweep up to 96)

cd "${SLURM_SUBMIT_DIR:-/mnt/home/xgao1/work/four_index_integral_solver/codes/multi_rhs}"
JULIA=$(ls -d ~/.julia/juliaup/julia-1.12*/bin/julia 2>/dev/null | sort -V | tail -1)   # manifest resolved with 1.12; 1.10 fails
[ -x "$JULIA" ] || JULIA=~/.juliaup/bin/julia
echo "host=$(hostname)  julia=$JULIA"
"$JULIA" --project -e 'using Pkg; Pkg.precompile()'   # settle the cache alone, before the array
rm -f bench_nd64_ops.jls
"$JULIA" --project bench_nd64.jl --build

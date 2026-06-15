#!/bin/bash
# PrecomputedVolumeField benchmark suite — runs on a DEDICATED (exclusive)
# 96-core genoa node so timings are free of the contention that corrupted the
# interactive worker7020 runs (there the decomposed pieces summed to ~2x the
# full call). Submit with the one-liner in the header comment below.
#
#   sbatch /mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/slurm/pvf_bench.sbatch
#
# Partition: `gen` is Rusty's default general CPU partition (genoa, 96c, 1.5TB),
# open to all users. Swap `-p gen` for your center's partition (e.g. -p ccq)
# if you want higher priority. Julia: juliaup binary only (never module load).

#SBATCH --job-name=pvf_bench
#SBATCH --partition=gen
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=01:30:00
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_pvf_bench_%j.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_pvf_bench_%j.err

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
WT=/mnt/home/xgao1/codes/BoundaryIntegral.jl-wt-pvf
J=/mnt/home/xgao1/.juliaup/bin/julia          # juliaup-managed Julia (NOT module load julia)
LOGS=$NR/exp65_orbital_bench/logs

export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96

echo "=== node $(hostname)  cores $(nproc)  $(date) ==="
echo "julia: $($J --version)"

run() {  # run() <script> <logfile>
    local script=$1 log=$2
    echo "--- $(date '+%H:%M:%S')  start $script"
    if $J --project=$WT $NR/exp65_orbital_bench/scripts/$script > $LOGS/$log 2>&1; then
        echo "--- $(date '+%H:%M:%S')  done  $script  ->  $log"
        grep -aE 'build|u_int|RHS|speedup|SPEEDUP|FULL ltkm3dc|type-1|type-2|scale|sum of|DONE|peak RSS' $LOGS/$log | tail -25
    else
        echo "--- $(date '+%H:%M:%S')  FAILED $script (see $log) — continuing"
    fi
}

# 1. warm old-vs-new A/B (build + RHS + u_int), VolumeSource path included
run bench_field_path.jl    slurm_bench_field.log
# 2. cache_fft mode (precomputed fine-grid FFT + interp-only)
run bench_cache_fft.jl     slurm_bench_cache_fft.log
# 3. clean decomposition of the old ltkm3dc call (exclusive node => sum ~= full)
run diag_ltkm_decomp.jl    slurm_diag_ltkm_decomp.log

echo "=== ALL BENCHMARKS DONE $(date) ==="

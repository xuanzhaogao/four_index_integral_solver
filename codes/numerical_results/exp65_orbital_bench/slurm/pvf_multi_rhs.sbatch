#!/bin/bash
# exp65 multi-RHS benchmark — center orbital + lattice-neighbor pair densities,
# block solve vs sequential, full cutoff sweep (K up to ~19). One exclusive 96-core
# genoa node. Runs at full 96 threads: the FINUFFT 96-thread FFTW plan pathology (which
# crippled BOTH the interface build and the four_index_matrix eval) is fixed IN-FUNCTION in
# TKM3D's ltkm3dc (FINUFFT nthreads capped to 16, env TKM3D_FINUFFT_NTHREADS) — so the FFT
# avoids the pathology while the FMM (lfmm3d) keeps all 96 threads. Build dropped from
# 131-145 s to ~21 s at large K; eval from ~46 to ~12 s/K. Requires the TKM3D
# nthreads-cap fix (branch finufft-nthreads-cap) active in the TKM3D the env resolves.
# Submit yourself:  sbatch slurm/pvf_multi_rhs.sbatch
#SBATCH --job-name=pvf_multi_rhs
#SBATCH --partition=gen
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_multi_rhs_%j.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp65_orbital_bench/logs/slurm_multi_rhs_%j.err

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
WT=/mnt/home/xgao1/codes/BoundaryIntegral.jl
J=/mnt/home/xgao1/.juliaup/bin/julia          # juliaup binary (never module load julia)
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96           # full FMM threads; the FINUFFT FFT is capped in-function (TKM3D ltkm3dc)
# export TKM3D_FINUFFT_NTHREADS=16 # in-function FFT cap (TKM3D default 16); set 0 for FINUFFT's omp default

echo "=== node $(hostname)  cores $(nproc)  $(date) ==="
$J --project=$WT $NR/exp65_orbital_bench/scripts/run_multi_rhs.jl
echo "=== DONE $(date) ==="

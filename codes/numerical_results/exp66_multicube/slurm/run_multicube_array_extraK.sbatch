#!/bin/bash
# Two extra K-sweep points for the multi-RHS figure: K=25 (cutoff 4.5) and K=31 (cutoff 5.0),
# i.e. one and two more graphene neighbor shells beyond the current K=19 (cutoff 3.9). These
# are CONVERGED production solves (default itmax) at 96 threads on exclusive genoa nodes --
# identical conditions to the original K-sweep (run_multicube_array.sbatch) so the panel-(b)
# points stay mutually consistent. Each task writes data/raw/multicube_K<K>_cut<cutoff>.jls,
# which plot_scripts/plot_scaling.jl picks up automatically for panel (b) (K-sweep glob).
#
# Submit (you):  sbatch exp66_multicube/slurm/run_multicube_array_extraK.sbatch
# After both finish: julia --project=. plot_scripts/plot_scaling.jl
#
#SBATCH --job-name=multicube_extraK
#SBATCH --partition=ccm
#SBATCH --constraint=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --array=0-1
#SBATCH --output=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp66_multicube/logs/slurm_multicube_extraK_%A_%a.out
#SBATCH --error=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results/exp66_multicube/logs/slurm_multicube_extraK_%A_%a.err

CUTOFFS=(4.5 5.0)                        # array index -> neighbor cutoff (K = 25, 31)
export MULTICUBE_CUTOFF=${CUTOFFS[$SLURM_ARRAY_TASK_ID]}

NR=/mnt/home/xgao1/work/four_index_integral_solver/codes/numerical_results
J=/mnt/home/xgao1/.juliaup/bin/julia    # juliaup binary (never module load julia)
export JULIA_NUM_THREADS=96
export OMP_NUM_THREADS=96                # full threads for the FMM; TKM3D caps FINUFFT at 16 internally

echo "=== array task $SLURM_ARRAY_TASK_ID  cutoff=$MULTICUBE_CUTOFF  node $(hostname)  $(date) ==="
$J --project=$NR $NR/exp66_multicube/scripts/run_multicube.jl
echo "=== DONE cutoff=$MULTICUBE_CUTOFF  $(date) ==="

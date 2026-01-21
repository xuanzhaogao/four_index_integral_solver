#!/bin/bash

# set name of job
#SBATCH --job-name=convergence_corrected_far
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err

#SBATCH -p ccm
#SBATCH -N1
#SBATCH -C genoa

cd /mnt/home/xgao1/work/four_index_integral_solver/codes/box3d_corrected

julia --project=/mnt/home/xgao1/work/four_index_integral_solver/codes/box3d_corrected -t auto convergence_corrected_far.jl
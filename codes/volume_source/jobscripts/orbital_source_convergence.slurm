#!/bin/bash

# set name of job
#SBATCH --job-name=orbital
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err

#SBATCH -p ccm
#SBATCH -N1
#SBATCH -C genoa

cd /mnt/home/xgao1/work/four_index_integral_solver/codes/volume_source

julia --project=/mnt/home/xgao1/work/four_index_integral_solver/codes/volume_source -t auto orbital_source_convergence.jl
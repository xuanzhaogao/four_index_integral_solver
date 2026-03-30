#!/bin/bash

#SBATCH --job-name=plot_graphene_slab
#SBATCH --output=jobscripts/job.%J.out
#SBATCH --error=jobscripts/job.%J.err

#SBATCH -p ccm
#SBATCH -N1
#SBATCH -C genoa

cd /mnt/home/xgao1/work/four_index_integral_solver/codes/multibox3d_convergence

module load gcc fftw openmpi

julia --project=. scripts/plot_graphene_slab.jl

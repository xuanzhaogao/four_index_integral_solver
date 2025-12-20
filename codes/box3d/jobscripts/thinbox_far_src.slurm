#!/bin/bash

# set name of job
#SBATCH --job-name=thinbox_far_src
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err

#SBATCH -p ccm
#SBATCH -N1
#SBATCH -C genoa

cd /mnt/home/xgao1/work/four_index_integral_solver/codes/box3d

OMP_NUM_THREADS=96 julia --project=/mnt/home/xgao1/work/four_index_integral_solver/codes/box3d -t auto single_thin_box_far_src.jl
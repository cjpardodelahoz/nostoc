#!/bin/bash

#SBATCH --array=1-90
#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/targz_6240_6526_7535_%A_%a.out
#SBATCH --error=log/targz_6240_6526_7535_%A_%a.err
#SBATCH --partition=scavenger

S=$(cat scripts/sample_ids | sed -n ${SLURM_ARRAY_TASK_ID}p)

tar -czf analyses/reads/${S}.tar.gz analyses/reads/${S}
rm -rf analyses/reads/${S}

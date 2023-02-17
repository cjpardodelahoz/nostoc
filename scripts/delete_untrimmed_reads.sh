#!/bin/bash

#SBATCH --array=1-90
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --partition=scavenger

S=$(cat scripts/sample_ids_6240_6526_7535 | sed -n ${SLURM_ARRAY_TASK_ID}p)

rm analyses/reads/${S}/${S}*.fastq

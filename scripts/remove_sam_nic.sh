#!/bin/bash

#SBATCH --array=1-120
#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/remove_sam_nic_%A_%a.out
#SBATCH --error=log/remove_sam_nic_%A_%a.err
#SBATCH --partition=scavenger

S=$(cat scripts/sample_ids_nic | sed -n ${SLURM_ARRAY_TASK_ID}p)

rm analyses/assemblies/${S}/${S}.sam

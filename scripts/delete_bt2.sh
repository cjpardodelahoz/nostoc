#!/bin/bash

#SBATCH --array=1-32
#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/delete_bt2_%A_%a.out
#SBATCH --error=log/delete_bt2_%A_%a.err
#SBATCH --partition=scavenger

S=$(cat scripts/sample_ids | sed -n ${SLURM_ARRAY_TASK_ID}p)

rm analyses/assemblies/${S}/${S}*.bt2
rm analyses/assemblies/${S}/${S}_sorted.bam
rm analyses/assemblies/${S}/${S}_assembly_depths.txt

#!/bin/bash

#SBATCH --array=1-32
#SBATCH --mem-per-cpu=8G
#SBATCH -c 1
#SBATCH --error=log/merge_reads_7535.err
#SBATCH --output=log/merge_reads_7535.out
#SBATCH --partition=common

S=$(cat scripts/sample_ids | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Make sample specific directory
mkdir analyses/reads/${S}
# Merge reads and sort into folders
cat analyses/reads/${S}*R1* > analyses/reads/${S}/${S}_R1_all.fastq
cat analyses/reads/${S}*R2* > analyses/reads/${S}/${S}_R2_all.fastq
# Delete copied read files
rm analyses/reads/${S}*R1*
rm analyses/reads/${S}*R2*

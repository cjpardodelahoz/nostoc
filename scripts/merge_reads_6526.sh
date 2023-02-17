#!/bin/bash

#SBATCH --array=1-30
#SBATCH --mem-per-cpu=16G
#SBATCH -c 1
#SBATCH --error=log/merge_reads_6526_%A_%a.err
#SBATCH --output=log/merge_reads_6526_%A_%a.out
#SBATCH --partition=scavenger

S=$(cat scripts/sample_ids_6526 | sed -n ${SLURM_ARRAY_TASK_ID}p)

rm -rf analyses/reads/${S}

# Copy reads into analyses directory
cp data/reads/6526/${S}*.gz analyses/reads/
# Make sample specific directory
mkdir analyses/reads/${S}
# Unzip reads
#gunzip analyses/reads/${S}*R1*
#gunzip analyses/reads/${S}*R2*
# Merge reads and sort into folders
mv analyses/reads/${S}*R1* analyses/reads/${S}/${S}_R1_all.fastq.gz
mv analyses/reads/${S}*R2* analyses/reads/${S}/${S}_R2_all.fastq.gz
# Delete copied read files
#rm analyses/reads/${S}*R1*
#rm analyses/reads/${S}*R2*

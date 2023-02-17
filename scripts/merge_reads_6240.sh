#!/bin/bash

#SBATCH --array=1-28
#SBATCH --mem-per-cpu=16G
#SBATCH -c 1
#SBATCH --error=log/merge_reads_6240_%A_%a.err
#SBATCH --output=log/merge_reads_6240_%A_%a.out
#SBATCH --partition=scavenger

S=$(cat scripts/sample_ids_6240 | sed -n ${SLURM_ARRAY_TASK_ID}p)

#rm -rf analyses/reads/${S}

# Copy reads into analyses directory
cp data/reads/6240_1_2/${S}*.gz analyses/reads/
# Make sample specific directory
mkdir analyses/reads/${S}
# Unzip reads
#gunzip analyses/reads/${S}*R1*
#gunzip analyses/reads/${S}*R2*
# Merge reads and sort into folders
cat analyses/reads/${S}*R1* > analyses/reads/${S}/${S}_R1_all.fastq.gz
cat analyses/reads/${S}*R2* > analyses/reads/${S}/${S}_R2_all.fastq.gz
# Delete copied read files
rm analyses/reads/${S}*R1*
rm analyses/reads/${S}*R2*

#!/bin/bash

#SBATCH --array=1-11
#SBATCH --mem-per-cpu=16G
#SBATCH -c 1
#SBATCH --error=log/merge_reads_intermediares_%A_%a.err
#SBATCH --output=log/merge_reads_intermediares_%A_%a.out
#SBATCH --partition=scavenger

S=$(cat scripts/sample_ids_intermediares | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Make sample specific directory
mkdir analyses/reads/${S}
# Merge reads and sort into folders
cat data/reads/nic_reads/reads_genomes_intermediaires/${S}*R1* > \
analyses/reads/${S}/${S}_R1_all.fastq.gz
cat data/reads/nic_reads/reads_genomes_intermediaires/${S}*R2* > \
analyses/reads/${S}/${S}_R2_all.fastq.gz

#!/bin/bash

#SBATCH --array=1-2
#SBATCH --mem-per-cpu=4G
#SBATCH -c 1
#SBATCH --error=log/copy_vecherskii_%A_%a.err
#SBATCH --output=log/copy_vecherskii_%A_%a.out
#SBATCH --partition=scavenger

S=$(cat scripts/sample_ids_vecherskii | sed -n ${SLURM_ARRAY_TASK_ID}p)

mkdir -p analyses/reads/${S}
cp data/reads/vecherskii/${S}_1.fastq analyses/reads/${S}/${S}_R1_all.fastq
cp data/reads/vecherskii/${S}_2.fastq analyses/reads/${S}/${S}_R2_all.fastq

#!/bin/bash

#SBATCH --array=1-2
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/samtosortedbam_vecherskii_%A_%a.out
#SBATCH --error=log/samtosortedbam_vecherskii_%A_%a.err
#SBATCH --partition=scavenger

module load samtools/1.9

S=$(cat scripts/sample_ids_vecherskii | sed -n ${SLURM_ARRAY_TASK_ID}p)

samtools view --threads 8 -b analyses/assemblies/${S}/${S}.sam | \
samtools sort -@ 4 -m 16G -o analyses/assemblies/${S}/${S}_sorted.bam

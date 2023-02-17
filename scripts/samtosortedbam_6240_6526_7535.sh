#!/bin/bash

#SBATCH --array=90
#SBATCH --mem-per-cpu=12G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/samtosortedbam_6240_6526_7535_%A_%a.out
#SBATCH --error=log/samtosortedbam_6240_6526_7535_%A_%a.err
#SBATCH --partition=scavenger

module load samtools/1.9

S=$(cat scripts/sample_ids_6240_6526_7535 | sed -n ${SLURM_ARRAY_TASK_ID}p)

samtools view --threads 4 -b analyses/assemblies/${S}/${S}.sam | \
samtools sort -@ 4 -m 12G -o analyses/assemblies/${S}/${S}_sorted.bam

#rm analyses/assemblies/${S}/${S}.sam

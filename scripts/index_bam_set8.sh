#!/bin/bash

#SBATCH --array=1-127
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/index_bam_%A_%a.out
#SBATCH --error=log/index_bam_%A_%a.err
#SBATCH --partition=scavenger

module load samtools/1.9

S=$(cat scripts/sample_ids_set8 | sed -n ${SLURM_ARRAY_TASK_ID}p)

samtools index -b -@ 4 analyses/assemblies/${S}/${S}_sorted.bam \
 analyses/assemblies/${S}/${S}_sorted.bam.bai
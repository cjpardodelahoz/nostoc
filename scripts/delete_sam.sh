#!/bin/bash

#SBATCH --array=1-32
#SBATCH --mem-per-cpu=25G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/bowtie2_map_%A_%a.out
#SBATCH --error=log/bowtie2_map_%A_%a.err
#SBATCH --partition=scavenger

#module load Bowtie2/2.3.5.1

S=$(cat scripts/sample_ids | sed -n ${SLURM_ARRAY_TASK_ID}p)

#bowtie2 --sensitive-local -p 8 --seed 1234  -x analyses/assemblies/${S}/${S}_assembly_db \
#-1 analyses/reads/${S}/${S}_R1_paired.fq.gz -2 analyses/reads/${S}/${S}_R2_paired.fq.gz -S \
#analyses/assemblies/${S}/${S}.sam

rm analyses/assemblies/${S}/${S}.sam
rm analyses/assemblies/${S}/${S}.bam

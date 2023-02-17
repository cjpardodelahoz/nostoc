#!/bin/bash

#SBATCH --array=1-120
#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/bt2_map_nic_%A_%a.out
#SBATCH --error=log/bt2_map_nic_%A_%a.err
#SBATCH --partition=common

module load Bowtie2/2.3.5.1

S=$(cat scripts/sample_ids_nic | sed -n ${SLURM_ARRAY_TASK_ID}p)

bowtie2 --sensitive-local -p 8 --seed 1234  -x analyses/assemblies/${S}/${S}_assembly_db \
-1 analyses/assemblies/${S}/corrected/${S}_R1_paired.fq*cor*.gz \
-2 analyses/assemblies/${S}/corrected/${S}_R2_paired.fq*cor*.gz -S analyses/assemblies/${S}/${S}.sam


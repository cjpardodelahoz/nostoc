#!/bin/bash

#SBATCH --array=1-120
#SBATCH --mem-per-cpu=2G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/bt2_db_nic_%A_%a.out
#SBATCH --error=log/bt2_db_nic_%A_%a.err
#SBATCH --partition=common

module load Python/3.8.1
module load Bowtie2/2.3.5.1

S=$(cat scripts/sample_ids_nic | sed -n ${SLURM_ARRAY_TASK_ID}p)

bowtie2-build --seed 1234 --threads 8 analyses/assemblies/${S}/contigs.fasta \
analyses/assemblies/${S}/${S}_assembly_db

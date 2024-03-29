#!/bin/bash

#SBATCH --array=1-112
#SBATCH --mem-per-cpu=12G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=log/assemble_metagenomes_%A_%a.out
#SBATCH --error=log/assemble_metagenomes_%A_%a.err
#SBATCH --partition=common

# Load SPAdes module or path
module load Python/3.8.1
module load SPAdes/3.14.1
# Variable with sample name
S=$(cat misc_files/read_accessions.txt | cut -f 2 | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Run metaSPAdes
spades.py --meta --pe1-1 analyses/reads/${S}/${S}_R1_paired.fq.gz \
    --pe1-2 analyses/reads/${S}/${S}_R2_paired.fq.gz \
    -k 55,75,95 \
    -t 16 \
    -o analyses/assemblies/${S}

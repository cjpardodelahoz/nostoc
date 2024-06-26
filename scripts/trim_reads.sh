#!/bin/bash

#SBATCH --array=1-112
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/trimmo_%A_%a.out
#SBATCH --error=log/trimmo_%A_%a.err
#SBATCH --partition=scavenger

# Trimmomatic module or path
module load Trimmomatic/0.39
# Variable with sample name
S=$(cat misc_files/read_accessions.txt | cut -f 2 | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Run Trimmomatic
java -jar $TRIMMOMATIC PE -threads 4 \
-summary analyses/reads/${S}/${S}_trimm_summary.txt \
analyses/reads/${S}/*R1_all.fastq.gz analyses/reads/${S}/*R2_all.fastq.gz \
analyses/reads/${S}/${S}_R1_paired.fq.gz analyses/reads/${S}/${S}_R1_unpaired.fq.gz \
analyses/reads/${S}/${S}_R2_paired.fq.gz analyses/reads/${S}/${S}_R2_unpaired.fq.gz \
ILLUMINACLIP:misc_files/TruSeq3-PE.fa:2:30:15:4:true LEADING:20 TRAILING:20 \
SLIDINGWINDOW:10:20 MINLEN:75

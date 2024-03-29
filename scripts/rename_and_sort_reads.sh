#!/bin/bash

#SBATCH --array=1-112
#SBATCH --mem-per-cpu=16G
#SBATCH -c 1
#SBATCH --error=log/rename_and_sort_reads_%A_%a.err
#SBATCH --output=log/rename_and_sort_reads_%A_%a.out
#SBATCH --partition=scavenger

# Variable with SRA accession
sra_accession=$(sed -n ${SLURM_ARRAY_TASK_ID}p misc_files/read_accessions.txt | cut -f 1)
# Variable with sample name
sample_name=$(sed -n ${SLURM_ARRAY_TASK_ID}p misc_files/read_accessions.txt | cut -f 2)

# Rename and sort reads
mkdir -p analyses/reads/${sample_name}
mv analyses/reads/${sra_accession}_1.fastq analyses/reads/${sample_name}/${sample_name}_R1_all.fastq
mv analyses/reads/${sra_accession}_2.fastq analyses/reads/${sample_name}/${sample_name}_R2_all.fastq
# Gzip the files
gzip analyses/reads/${sample_name}/${sample_name}_R1_all.fastq
gzip analyses/reads/${sample_name}/${sample_name}_R2_all.fastq

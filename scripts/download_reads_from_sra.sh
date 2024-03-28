#!/bin/bash

#SBATCH --mem-per-cpu=16G
#SBATCH -c 1
#SBATCH --error=log/download_reads.err
#SBATCH --output=log/download_reads.out
#SBATCH --partition=scavenger

# Set path for SRA toolkits
export PATH=/hpc/group/bio1/carlos/apps/sratoolkit.3.0.0-centos_linux64/bin:$PATH
# Make directory for raw reads
mkdir -p analyses/reads
# Download raw reads from SRA
prefetch -O analyses/reads $(cat misc_files/read_accessions.txt | cut -f 1) 
fasterq-dump -O analyses/reads $(cat misc_files/read_accessions.txt | cut -f 1)

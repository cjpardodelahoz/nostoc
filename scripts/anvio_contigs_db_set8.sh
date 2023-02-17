#!/bin/bash

#SBATCH --array=1-127
#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/anvio_contigs_db_set8_%A_%a.out
#SBATCH --error=log/anvio_contigs_db_set8_%A_%a.err
#SBATCH --partition=scavenger

source $(conda info --base)/etc/profile.d/conda.sh
conda activate anvio-7.1

S=$(cat scripts/sample_ids_set8 | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Directory for anvio files
mkdir -p analyses/bins/anvio/${S}
# Generate contigs database. This calculates 4-mer freqs, and identifies ORFs
anvi-gen-contigs-database -f analyses/assemblies/${S}/contigs.fasta \
 -o analyses/bins/anvio/${S}/contigs.db \
 -n ${S} \
 --num-threads 8 \
 --split-length 0
# Run default HMMs
anvi-run-hmms -c analyses/bins/anvio/${S}/contigs.db \
 --num-threads 8
# Set up ncbi cogs - first run only 
# anvi-setup-ncbi-cogs
# Annotate genes with COGs
anvi-run-ncbi-cogs -c analyses/bins/anvio/${S}/contigs.db \
 --num-threads 8
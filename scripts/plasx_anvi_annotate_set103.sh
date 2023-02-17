#!/bin/bash

#SBATCH --array=1-151
#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/plasx_anvi_annotate_set103_%A_%a.out
#SBATCH --error=log/plasx_anvi_annotate_set103_%A_%a.err
#SBATCH --partition=scavenger

source $(conda info --base)/etc/profile.d/conda.sh
conda activate anvio-7.1

# Genome ID variable
S=$(cat scripts/genome_ids_set103 | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Output directory
PREFIX="analyses/plasmid_detection/set103/${S}"
# Path to genomes
GENOME_PATH="analyses/cyano_genomes/set103"
# Threads
THREADS=8
# Directory for anvio files
mkdir -p ${PREFIX}
# Simplyfi contig names
anvi-script-reformat-fasta ${GENOME_PATH}/${S} \
 -o ${PREFIX}/contigs.fasta \
 --simplify-names \
 --seq-type NT \
 -r ${PREFIX}/header_key.txt
# Create an anvio contigs database from the fasta file
# - The `-L 0` parameter ensures that contigs remain intact and aren't split
anvi-gen-contigs-database -L 0 -T ${THREADS} -f ${PREFIX}/contigs.fasta \
 -o ${PREFIX}/contigs.db
# Export gene calls (including amino acid sequences) to text file
anvi-export-gene-calls --gene-caller prodigal -c ${PREFIX}/contigs.db \
 -o ${PREFIX}/gene_calls.txt
# Annotate COGs
anvi-run-ncbi-cogs -T ${THREADS} --cog-version COG14 --cog-data-dir \
 scripts/COG_2014 -c ${PREFIX}/contigs.db
# Annotate Pfams
anvi-run-pfams -T ${THREADS} --pfam-data-dir scripts/Pfam_v32 \
 -c ${PREFIX}/contigs.db
# Export functions to text file
anvi-export-functions --annotation-sources COG14_FUNCTION,Pfam \
 -c ${PREFIX}/contigs.db -o ${PREFIX}/cogs_and_pfams.txt
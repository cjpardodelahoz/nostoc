#!/bin/bash

#SBATCH --array=1-21
#SBATCH --mem-per-cpu=2G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/add_extra_contigs_from_graph_set9_%A_%a.out
#SBATCH --error=log/add_extra_contigs_from_graph_set9_%A_%a.err
#SBATCH --partition=scavenger

source $(conda info --base)/etc/profile.d/conda.sh
conda activate anvio-7.1

# Variable with Sample id
S=$(cat scripts/sample_ids_for_extra_contigs_from_graph | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Remove the "+" and "-" from bandage
sed -i 's/+//' analyses/genome_qc/set8/delta_contigs/${S}*/extra_contigs_from_graph.fasta
sed -i 's/-//' analyses/genome_qc/set8/delta_contigs/${S}*/extra_contigs_from_graph.fasta
# Add extra_contigs_from_graph.fasta file to the curated bin
cat analyses/genome_qc/set8/delta_contigs/${S}*/extra_contigs_from_graph.fasta >> \
 analyses/cyano_genomes/set10/${S}*.fa

#!/bin/bash

#SBATCH --array=1-1517
#SBATCH --mem-per-cpu=2G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/ml_gene_trees_set103_%A_%a.out
#SBATCH --error=log/ml_gene_trees_set103_%A_%a.err
#SBATCH --partition=scavenger

module load IQ-TREE/1.6.12

seq=$(cat scripts/busco_ids_filtered_set103 | sed -n ${SLURM_ARRAY_TASK_ID}p)

iqtree -nt 1 -pre analyses/phylogenetics/set103/trees/single/${seq} \
-s analyses/phylogenetics/set103/alignments/single/${seq}_ng.fna \
-m MFP+MERGE -spp analyses/phylogenetics/set103/alignments/single/${seq}_ng_codon_partition \
-bb 1000 -bnni -safe
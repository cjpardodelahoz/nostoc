#!/bin/bash

#SBATCH --array=1-773
#SBATCH --mem-per-cpu=2G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/ml_gene_trees_set73.%A_%a.out
#SBATCH --error=log/ml_gene_trees_set73.%A_%a.err
#SBATCH --partition=scavenger

module load IQ-TREE/1.6.12

seq=$(cat scripts/busco_ids_set73 | sed -n ${SLURM_ARRAY_TASK_ID}p)

iqtree -nt 1 -pre analyses/phylogenetics/set73/trees/single/${seq} \
-s analyses/phylogenetics/set73/alignments/single/${seq}_aln.fna \
-m GTR+G -bb 1000
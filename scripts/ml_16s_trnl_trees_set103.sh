#!/bin/bash

#SBATCH --mem-per-cpu=2G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/ml_16s_trnl_trees_set103.out
#SBATCH --error=log/ml_16_trnl_trees_set103.err
#SBATCH --partition=scavenger

module load IQ-TREE/1.6.12

# 16s tree
iqtree -nt 1 -pre analyses/phylogenetics/set103/trees/single/16s \
-s analyses/phylogenetics/set103/alignments/single/16s_aln.fas \
-m MFP -bb 1000 -bnni -safe

# trnl tree
iqtree -nt 1 -pre analyses/phylogenetics/set103/trees/single/trnl \
-s analyses/phylogenetics/set103/alignments/single/trnl_aln.fas \
-m MFP -bb 1000 -bnni -safe
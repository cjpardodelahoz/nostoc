#!/bin/bash

#SBATCH --mem-per-cpu=2G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/ml_rbclx_focal_groups.out
#SBATCH --error=log/ml_rbclx_focal_groups.err
#SBATCH --partition=scavenger

# Load IQTree module
module load IQ-TREE/1.6.12
# Make output directory
mkdir -p analyses/species_delimitation/cooccurrence/trees/focal
# Phylogroup V and relatives tree
iqtree -nt 1 \
 -pre analyses/species_delimitation/cooccurrence/trees/focal/rbclx_set103_global_abmi_v \
 -s analyses/species_delimitation/cooccurrence/seqs/rbclx_set103_global_abmi_v_aln.fna \
 -m MFP -bb 1000 -bnni -safe
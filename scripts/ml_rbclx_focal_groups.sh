#!/bin/bash

#SBATCH --array=1-16
#SBATCH --mem-per-cpu=2G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/ml_rbclx_focal_groups_%A_%a.out
#SBATCH --error=log/ml_rbclx_focal_groups_%A_%a.err
#SBATCH --partition=scavenger

# Define veriable with sequences
aln=$(cat misc_files/edited_alns.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
prefix=$(echo ${aln##*/} | sed "s|_aln_edited||" | sed "s|.fna||")
# Load IQTree module
module load IQ-TREE/1.6.12
# Phylogroup V and relatives tree
iqtree -nt 1 \
 -pre analyses/species_delimitation/rbclx/clade_assignment/trees/focal/${prefix} \
 -s ${aln} \
 -m MFP -bb 1000 -bnni -safe
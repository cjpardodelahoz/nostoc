#!/bin/bash

#SBATCH --array=1-17
#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/mafft_rbclx_focal_groups_%A_%a.out
#SBATCH --error=log/mafft_rbclx_focal_groups_%A_%a.err
#SBATCH --partition=scavenger

# Define veriable with sequences
clade=$(cat misc_files/clade_list.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Path to MAFFT
export PATH=/hpc/home/cjp47/mafft-7.475-with-extensions/bin:${PATH}
# Align focal group seqs
mafft --retree 1 --maxiterate 0 --adjustdirection \
 analyses/species_delimitation/rbclx/clade_assignment/seqs/${clade}.fna > \
 analyses/species_delimitation/rbclx/clade_assignment/alignments/${clade}_aln.fna

#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/mafft_rbclx_focal_groups.out
#SBATCH --error=log/mafft_rbclx_focal_groups.err
#SBATCH --partition=scavenger

# Path to MAFFT
export PATH=/hpc/home/cjp47/mafft-7.475-with-extensions/bin:${PATH}
# Align Phylogroup V and allies seqs
mafft --retree 1 --maxiterate 0 --adjustdirection \
 analyses/species_delimitation/cooccurrence/seqs/rbclx_set103_global_abmi_v.fna > \
 analyses/species_delimitation/cooccurrence/seqs/rbclx_set103_global_abmi_v_aln.fna

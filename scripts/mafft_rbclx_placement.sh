#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/mafft_rbclx_placement.out
#SBATCH --error=log/mafft_rbclx_placement.err
#SBATCH --partition=scavenger

export PATH=/hpc/home/cjp47/mafft-7.475-with-extensions/bin:${PATH}


mafft --retree 1 --maxiterate 0 --adjustdirection \
 analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_set103_abmi_public.fna > \
 analyses/species_delimitation/rbclx/clade_assignment/alignments/rbclx_set103_abmi_public_aln.fna

#!/bin/bash

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/concatenate_rbclx.out
#SBATCH --error=log/concatenate_rbclx.err
#SBATCH --partition=scavenger

export PATH=/hpc/group/bio1/carlos/apps/AMAS/amas:${PATH}

AMAS.py concat -i analyses/species_delimitation/cooccurrence/seqs/*set103.fna \
 -f fasta -d dna -p analyses/species_delimitation/cooccurrence/seqs/partition.txt \
 --concat-out analyses/species_delimitation/cooccurrence/seqs/rbclx_set103.fna \
 --part-format raxml
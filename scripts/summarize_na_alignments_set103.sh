#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/summarize_na_alignments_set103.out
#SBATCH --error=log/summarize_na_alignments_set103.err
#SBATCH --partition=scavenger

export PATH=/hpc/group/bio1/carlos/apps/AMAS/amas:${PATH}

AMAS.py summary -i analyses/phylogenetics/set103/alignments/single/*_ng.fna \
 -f fasta \
 -d dna \
 -c 8 \
 -o analyses/phylogenetics/set103/alignments/sumaries/summary_single_na_ng.txt
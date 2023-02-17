#!/bin/bash

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/concatenate_set73.out
#SBATCH --error=log/concatenate_set73.err
#SBATCH --partition=scavenger

export PATH=/hpc/group/bio1/carlos/apps/AMAS/amas:${PATH}

mkdir -p analyses/phylogenetics/set73/alignments/concat

AMAS.py concat -i analyses/phylogenetics/set73/alignments/single/*ng.fna \
-f fasta -d dna -p analyses/phylogenetics/set73/alignments/concat/Gpart_ng_na \
--part-format raxml --concat-out analyses/phylogenetics/set73/alignments/concat/concat_ng.fna

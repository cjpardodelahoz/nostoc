#!/bin/bash

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/concatenate_set1.out
#SBATCH --error=log/concatenate_set1.err
#SBATCH --partition=common

export PATH=/hpc/group/bio1/carlos/apps/AMAS/amas:${PATH}

AMAS.py concat -i analyses/phylogenetics/set1/alignments/single/*ng.fna \
-f fasta -d dna -p analyses/phylogenetics/set1/alignments/concat/Gpart_ng_na \
--part-format raxml --concat-out analyses/phylogenetics/set1/alignments/concat/concat_ng.fna

AMAS.py concat -i analyses/phylogenetics/set1/alignments/single/*ng.faa \
-f fasta -d aa -p analyses/phylogenetics/set1/alignments/concat/Gpart_ng_aa \
--part-format raxml --concat-out analyses/phylogenetics/set1/alignments/concat/concat_ng.faa

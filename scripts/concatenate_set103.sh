#!/bin/bash

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/concatenate_set103.out
#SBATCH --error=log/concatenate_set103.err
#SBATCH --partition=scavenger

export PATH=/hpc/group/bio1/carlos/apps/AMAS/amas:${PATH}

mkdir -p analyses/phylogenetics/set103/alignments/concat
# Make variable with paths to filtered alignments
aln_dir=analyses/phylogenetics/set103/alignments/single/
aln_paths=$(cat scripts/busco_ids_filtered_set103 | sed 's/$/_ng.fna/' |  sed '1i 16s_aln.fas' | sed '1i trnl_aln.fas' | sed "s|^|${aln_dir}|")

AMAS.py concat -i ${aln_paths} \
-f fasta -d dna -p analyses/phylogenetics/set103/alignments/concat/gene_partition_ng_na \
--part-format raxml --concat-out analyses/phylogenetics/set103/alignments/concat/concat_ng.fna

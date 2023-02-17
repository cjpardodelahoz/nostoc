#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=log/concat_pf_set103.out
#SBATCH --error=log/concat_pf_set103.err
#SBATCH --partition=common

module load IQ-TREE/1.6.12

mkdir -p analyses/phylogenetics/set103/trees/concat

iqtree -nt 16 -s analyses/phylogenetics/set103/alignments/concat/concat_ng.fna \
 -m MF+MERGE -rclusterf 10 -rcluster-max 100 \
 -spp analyses/phylogenetics/set103/alignments/concat/codon_partition_concat_ng_na \
 -pre analyses/phylogenetics/set103/trees/concat/concat_pf_ng_na

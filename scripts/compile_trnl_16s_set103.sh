#!/bin/bash

#SBATCH --mem-per-cpu=4G
#SBATCH -c 1
#SBATCH --error=log/compile_trnl_16s_set103.err
#SBATCH --output=log/compile_trnl_16s_set103.out
#SBATCH --partition=scavenger

bins=$(cat scripts/genome_ids_set103)

# Compile 16s seqs
touch analyses/phylogenetics/set103/seqs/16s.fas
for bin in ${bins} ; do
 if [[ analyses/cyano_genomes_annotation/set103/${bin}/16s.fas  ]] ; then
  cat analyses/cyano_genomes_annotation/set103/${bin}/16s.fas >> \
      analyses/phylogenetics/set103/seqs/16s.fas
 fi
done
# Compile trnl seqs
touch analyses/phylogenetics/set103/seqs/trnl.fas
for bin in ${bins} ; do
 if [[ analyses/cyano_genomes_annotation/set103/${bin}/trnl.fas  ]] ; then
  cat analyses/cyano_genomes_annotation/set103/${bin}/trnl.fas >> \
      analyses/phylogenetics/set103/seqs/trnl.fas
 fi
done

#!/bin/bash

#SBATCH --array=1-773
#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/mafft_set1_aa_%A_%a.out
#SBATCH --error=log/mafft_set1_aa_%A_%a.err
#SBATCH --partition=scavenger

export PATH=/hpc/home/cjp47/mafft-7.475-with-extensions/bin:${PATH}

seq=$(cat scripts/busco_ids_set1 | sed -n ${SLURM_ARRAY_TASK_ID}p)

mafft --maxiterate 1000 --globalpair \
analyses/phylogenetics/set1/seqs/${seq}.faa > \
analyses/phylogenetics/set1/alignments/single/${seq}_aln.faa

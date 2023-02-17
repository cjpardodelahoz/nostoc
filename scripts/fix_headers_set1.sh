#!/bin/bash

#SBATCH --array=1-773
#SBATCH --mem-per-cpu=1G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/fix_headers_set1_%A_%a.out
#SBATCH --error=log/fix_headers_set1_%A_%a.err
#SBATCH --partition=scavenger

seq=$(cat scripts/busco_ids_set1 | sed -n ${SLURM_ARRAY_TASK_ID}p)

sed -i 's/.fa.*/.fa/' analyses/phylogenetics/set1/alignments/single/${seq}_ng.faa
sed -i 's/.fa.*/.fa/' analyses/phylogenetics/set1/alignments/single/${seq}_ng.fna

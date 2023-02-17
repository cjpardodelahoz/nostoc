#!/bin/bash

#SBATCH --array=1-1648
#SBATCH --mem-per-cpu=1G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/fix_headers_set200_%A_%a.out
#SBATCH --error=log/fix_headers_set200_%A_%a.err
#SBATCH --partition=scavenger

seq=$(cat scripts/busco_ids_set200 | sed -n ${SLURM_ARRAY_TASK_ID}p)

sed -i 's/ .*//' analyses/phylogenetics/set200/alignments/single/${seq}_ng.faa
#!/bin/bash

#SBATCH --array=1-127
#SBATCH --mem-per-cpu=4G
#SBATCH -c 8
#SBATCH --error=log/label_contigs_set7_%A_%a.err
#SBATCH --output=log/label_contigs_set7_%A_%a.out
#SBATCH --partition=scavenger

S=$(cat scripts/sample_ids_set7 | sed -n ${SLURM_ARRAY_TASK_ID}p)

scripts/label_binned_contigs.sh analyses/bins/metabat/${S} ${S} analyses/bins/graphbin/binned_csv


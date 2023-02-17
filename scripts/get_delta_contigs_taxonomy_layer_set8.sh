#!/bin/bash

#SBATCH --array=1-127
#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/get_delta_contigs_taxonomy_layer_set8_%A_%a.out
#SBATCH --error=log/get_delta_contigs_taxonomy_layer_set8_%A_%a.err
#SBATCH --partition=scavenger

S=$(cat scripts/sample_ids_set8 | sed -n ${SLURM_ARRAY_TASK_ID}p)

module load R/4.1.1-rhel8

scripts/get_delta_contigs_taxonomy_layer.R analyses/bins/anvio/${S}/graphbin_delta_collection_with_header_splitnames.txt \
 analyses/mmseq_taxonomy/${S}/taxonomy_out.tsv \
 analyses/bins/anvio/${S}


#!/bin/bash

#SBATCH --array=1-127
#SBATCH --mem-per-cpu=32G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=log/mmseqs_taxonomy_set8_%A_%a.out
#SBATCH --error=log/mmseqs_taxonomy_set8_%A_%a.err
#SBATCH --partition=scavenger

source $(conda info --base)/etc/profile.d/conda.sh
conda activate mmseqs2

S=$(cat scripts/sample_ids_set8 | sed -n ${SLURM_ARRAY_TASK_ID}p)

mmseqs taxonomy analyses/mmseq_taxonomy/${S}/contigs_db \
 scripts/uniref90/uniref90_db \
 analyses/mmseq_taxonomy/${S}/taxonomy_out \
 analyses/mmseq_taxonomy/${S}/tmp \
 --threads 16


#!/bin/bash

#SBATCH --array=1-127
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/mmseqs_query_dbs_set8_%A_%a.out
#SBATCH --error=log/mmseqs_query_dbs_set8_%A_%a.err
#SBATCH --partition=scavenger

source $(conda info --base)/etc/profile.d/conda.sh
conda activate mmseqs2

S=$(cat scripts/sample_ids_set8 | sed -n ${SLURM_ARRAY_TASK_ID}p)

mkdir -p analyses/mmseq_taxonomy/${S}

mmseqs createdb analyses/assemblies/${S}/contigs.fasta \
 analyses/mmseq_taxonomy/${S}/contigs_db
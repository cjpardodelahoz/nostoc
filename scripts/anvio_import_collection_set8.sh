#!/bin/bash

#SBATCH --array=1-127
#SBATCH --mem-per-cpu=12G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/anvio_import_graphbin_collection_set8_%A_%a.out
#SBATCH --error=log/anvio_import_graphbin_collection_set8_%A_%a.err
#SBATCH --partition=scavenger

source $(conda info --base)/etc/profile.d/conda.sh
conda activate anvio-7.1

S=$(cat scripts/sample_ids_set8 | sed -n ${SLURM_ARRAY_TASK_ID}p)

anvi-import-collection analyses/bins/graphbin/${S}/graphbin2_anvio_collection.txt \
 -c analyses/bins/anvio/${S}/contigs.db \
 -p analyses/bins/anvio/${S}/profile_1/PROFILE.db \
 --contigs-mode
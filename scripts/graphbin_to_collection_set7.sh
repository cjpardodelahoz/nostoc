#!/bin/bash

#SBATCH --array=1-127
#SBATCH --mem-per-cpu=4G
#SBATCH -c 1
#SBATCH --error=log/graphbin_to_collection_set7_%A_%a.err
#SBATCH --output=log/graphbin_to_collection_set7_%A_%a.out
#SBATCH --partition=scavenger

S=$(cat scripts/sample_ids_set7 | sed -n ${SLURM_ARRAY_TASK_ID}p)


cat analyses/bins/graphbin/${S}/graphbin2_output.csv | sed 's/,/\t/' > \
 analyses/bins/graphbin/${S}/graphbin2_anvio_collection.txt




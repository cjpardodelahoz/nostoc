#!/bin/bash

#SBATCH --array=1-121
#SBATCH --mem-per-cpu=12G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/anvio_export_curated_cyano_bins_%A_%a.out
#SBATCH --error=log/anvio_export_curated_cyano_bins_%A_%a.err
#SBATCH --partition=scavenger

source $(conda info --base)/etc/profile.d/conda.sh
conda activate anvio-7.1

S=$(cat scripts/sample_ids_anvio_edited | sed -n ${SLURM_ARRAY_TASK_ID}p)

anvi-summarize -c analyses/bins/anvio/${S}/contigs.db \
 -p analyses/bins/anvio/${S}/profile_1/PROFILE.db \
 --collection-name graphbin_delta \
 -o analyses/bins/anvio/${S}/summary_out
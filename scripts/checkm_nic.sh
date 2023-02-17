#!/bin/bash

#SBATCH --array=1-120
#SBATCH --mem-per-cpu=10G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/checkm_nic_%A_%a.out
#SBATCH --error=log/checkm_nic_%A_%a.err
#SBATCH --partition=scavenger


source $(conda info --base)/etc/profile.d/conda.sh
conda activate checkm

S=$(cat scripts/sample_ids_nic | sed -n ${SLURM_ARRAY_TASK_ID}p)

checkm lineage_wf --threads 8 -x fa --tab_table \
--file analyses/genome_qc/metabat/nic/${S}/checkm_out \
analyses/bins/metabat/${S} analyses/genome_qc/metabat/nic/${S}

#!/bin/bash

#SBATCH --mem-per-cpu=10G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=log/checkm_set1.out
#SBATCH --error=log/checkm_set1.err
#SBATCH --partition=scavenger

source $(conda info --base)/etc/profile.d/conda.sh
conda activate checkm

checkm lineage_wf --threads 16 -x fa --tab_table \
--file analyses/genome_qc/set1/checkm/checkm_out \
analyses/cyano_genomes/set1 analyses/genome_qc/set1/checkm

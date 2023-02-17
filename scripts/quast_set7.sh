#!/bin/bash

#SBATCH --array=1-129
#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/quast_set7_%A_%a.out
#SBATCH --error=log/quast_set7_%A_%a.err
#SBATCH --partition=scavenger

source $(conda info --base)/etc/profile.d/conda.sh
conda activate quast

bin=$(cat scripts/genome_ids_set7 | sed -n ${SLURM_ARRAY_TASK_ID}p)
S=$(echo ${bin%_*})

quast --threads 8 -o analyses/genome_qc/set7/quast/${bin} --glimmer \
--rna-finding --pe1 analyses/assemblies/${S}/corrected/${S}_R1_paired.fq*cor*.gz \
--pe2 analyses/assemblies/${S}/corrected/${S}_R2_paired.fq*cor*.gz \
--no-check --no-sv analyses/cyano_genomes/set7/${bin}

#!/bin/bash

#SBATCH --array=1-151
#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/prokka_set103_%A_%a.out
#SBATCH --error=log/prokka_set103_%A_%a.err
#SBATCH --partition=scavenger

source $(conda info --base)/etc/profile.d/conda.sh
conda activate prokka

bin=$(cat scripts/genome_ids_set103 | sed -n ${SLURM_ARRAY_TASK_ID}p)

prokka --cpus 8 --metagenome --prefix ${bin} \
--outdir analyses/cyano_genomes_annotation/set103/${bin} --centre dukegcb \
--compliant --force analyses/cyano_genomes/set103/${bin}

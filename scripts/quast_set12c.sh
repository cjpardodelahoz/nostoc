#!/bin/bash

#SBATCH --array=1-151
#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/quast_set12c_%A_%a.out
#SBATCH --error=log/quast_set12c_%A_%a.err
#SBATCH --partition=scavenger

# Quast conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate quast
# Variable with genome filenames
genome=$(cat misc_files/genome_ids_set12c | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Run quast
quast --threads 8 -o analyses/genome_qc/set12c/quast/${genome} --glimmer \
 --rna-finding --no-sv analyses/cyano_genomes/set12c/${genome}

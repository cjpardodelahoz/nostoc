#!/bin/bash

#SBATCH --array=1-151
#SBATCH --mem-per-cpu=10G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=log/busco_set12c_%A_%a.out
#SBATCH --error=log/busco_set12c_%A_%a.err
#SBATCH --partition=scavenger

# Conda environment with BUSCO
source $(conda info --base)/etc/profile.d/conda.sh
conda activate busco
# Variable with genome filenames
genome=$(cat misc_files/genome_ids_set12c | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Run busco
busco -i analyses/cyano_genomes/set12c/${genome} -l nostocales_odb10 -m genome \
 --out analyses/genome_qc/set12c/busco/by_taxon/${genome} \
 --offline --download_path databases/busco_downloads

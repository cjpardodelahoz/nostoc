#!/bin/bash

#SBATCH --array=1-129
#SBATCH --mem-per-cpu=10G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=log/busco_set7_%A_%a.out
#SBATCH --error=log/busco_set7_%A_%a.err
#SBATCH --partition=scavenger

source $(conda info --base)/etc/profile.d/conda.sh
conda activate busco

S=$(cat scripts/genome_ids_set7 | sed -n ${SLURM_ARRAY_TASK_ID}p)

busco -i analyses/cyano_genomes/set7/${S} -l nostocales_odb10 -m genome \
--out analyses/genome_qc/set7/busco/by_taxon/${S} \
--offline --download_path scripts/busco_downloads

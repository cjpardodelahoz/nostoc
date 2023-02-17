#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/gunc_set1.out
#SBATCH --error=log/gunc_set1.err
#SBATCH --partition=common

source $(conda info --base)/etc/profile.d/conda.sh
conda activate gunc

gunc run --input_dir analyses/cyano_genomes/set1/ \
-r scripts/gunc_db/gunc_db_progenomes2.1.dmnd --use_species_level \
--out_dir analyses/genome_qc/set1/gunc \
--temp_dir analyses/genome_qc/set1/gunc \

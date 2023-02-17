#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/gunc_set8.out
#SBATCH --error=log/gunc_set8.err
#SBATCH --partition=common

source $(conda info --base)/etc/profile.d/conda.sh
conda activate gunc

mkdir -p analyses/genome_qc/set8/gunc

gunc run --input_dir analyses/cyano_genomes/set8/ \
-r scripts/gunc_db/gunc_db_progenomes2.1.dmnd --use_species_level \
--out_dir analyses/genome_qc/set8/gunc \
--temp_dir analyses/genome_qc/set8/gunc \

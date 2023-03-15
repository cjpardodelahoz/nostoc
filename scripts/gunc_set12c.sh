#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/gunc_set12c.out
#SBATCH --error=log/gunc_set12c.err
#SBATCH --partition=scavenger

# Conda environment with BUSCO
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gunc
# Out directory
mkdir -p analyses/genome_qc/set12c/gunc
# Run GUNC
gunc run --input_dir analyses/cyano_genomes/set12c/ \
 -r databases/gunc_db/gunc_db_progenomes2.1.dmnd --use_species_level \
 --out_dir analyses/genome_qc/set12c/gunc \
 --temp_dir /work/cjp47 \

#!/bin/bash

#SBATCH --mem-per-cpu=50G
#SBATCH -c 12
#SBATCH --error=log/fastani_set12c.err
#SBATCH --output=log/fastani_set12c.out
#SBATCH --partition=scavenger

module load FastANI/1.31

# Create query list file with path to genome bins
ls -d analyses/cyano_genomes/set12c/* > misc_files/fastani_set12c_ql
# Run FastANI with an all-by-all comparison
fastANI -t 12 --ql misc_files/fastani_set12c_ql --rl misc_files/fastani_set12c_ql \
 --matrix -o analyses/species_delimitation/fastani/set12c/fastani_set12c_ql_out

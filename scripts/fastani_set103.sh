#!/bin/bash

#SBATCH --mem-per-cpu=50G
#SBATCH -c 12
#SBATCH --error=log/fastani_set103.err
#SBATCH --output=log/fastani_set103.out
#SBATCH --partition=scavenger

module load FastANI/1.31

# Create query list file with path to genome bins
ls -d analyses/cyano_genomes/set103/* > scripts/fastani_set103_ql
# Run FastANI with an all-by-all comparison
fastANI -t 12 --ql scripts/fastani_set103_ql --rl scripts/fastani_set103_ql \
 --matrix -o analyses/species_delimitation/fastani/set103/fastani_set103_ql_out

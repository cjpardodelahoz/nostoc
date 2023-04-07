#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 32 # number of threads per process
#SBATCH --output=log/popcogent_set12c.out
#SBATCH --error=log/popcogent_set12c.err
#SBATCH --partition=scavenger

# Load PopCOGenT conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate PopCOGenT
# Run the PopCOGenT script
sh scripts/run_popcogent_set12c.sh
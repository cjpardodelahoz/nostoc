#!/bin/bash

#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/download_uniref90.out
#SBATCH --error=log/download_uniref90.err
#SBATCH --partition=scavenger

source $(conda info --base)/etc/profile.d/conda.sh
conda activate mmseqs2

mkdir -p scripts/uniref90

mmseqs databases UniRef90 scripts/uniref90/uniref90_db scripts/uniref90/tmp \
 --threads 4
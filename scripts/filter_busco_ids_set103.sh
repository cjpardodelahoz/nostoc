#!/bin/bash

#SBATCH --mem-per-cpu=32G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/filter_busco_ids_set103.out
#SBATCH --error=log/filter_busco_ids_set103.err
#SBATCH --partition=scavenger

module load R/4.1.1-rhel8

scripts/filter_busco_ids_set103.R
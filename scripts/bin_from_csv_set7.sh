#!/bin/bash

#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/bin_from_csv_set7.out
#SBATCH --error=log/bin_from_csv_set7.err
#SBATCH --partition=scavenger

module load R/4.1.1-rhel8

scripts/bin_from_csv_set7.R
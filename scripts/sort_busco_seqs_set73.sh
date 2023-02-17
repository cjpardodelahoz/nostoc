#!/bin/bash

#SBATCH --mem-per-cpu=32G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/sort_busco_seqs_set73.out
#SBATCH --error=log/sort_busco_seqs_set73.err
#SBATCH --partition=scavenger

module load R/4.1.1-rhel8

scripts/sort_busco_seqs_set73.R
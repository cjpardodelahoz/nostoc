#!/bin/bash

#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/report_2.out
#SBATCH --error=log/report_2.err
#SBATCH --partition=scavenger

module load R/4.1.1-rhel8

scripts/report_2.R
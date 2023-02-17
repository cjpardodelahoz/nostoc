#!/bin/bash

#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 11 # number of threads per process
#SBATCH --output=log/quast_install.out
#SBATCH --error=log/quast_install.err
#SBATCH --partition=scavenger

mamba create -n quast -c bioconda -c conda-forge quast

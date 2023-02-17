#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=log/wastral_set103.out
#SBATCH --error=log/wastral_set103.err
#SBATCH --partition=scavenger

# Load path to weighted astral hybrid 
export PATH=/hpc/group/bio1/carlos/apps/ASTER-Linux/bin:$PATH
# Run ASTRAL weighted hybrid
astral-hybrid -t 16 \
 -i analyses/phylogenetics/set103/trees/astral/ml_gene.trees \
 -o analyses/phylogenetics/set103/trees/astral/wastral.tree
 
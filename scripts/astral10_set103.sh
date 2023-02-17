#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/astral10_set103.out
#SBATCH --error=log/astral10_set103.err
#SBATCH --partition=scavenger

# Load newick utilities conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate newick
# Load java for astral
module load Java/1.8.0_60
astral_path=/hpc/group/bio1/carlos/apps/Astral

# Collapse nodes with less than 10% UFB. Has to run after running astral_set103.sh
nw_ed analyses/phylogenetics/set103/trees/astral/ml_gene.trees 'i & b<=10' \
 o > analyses/phylogenetics/set103/trees/astral/ml10_gene.trees
# Run ASTRAL
java -jar ${astral_path}/astral.5.15.5.jar \
 -T 8 \
 -i analyses/phylogenetics/set103/trees/astral/ml10_gene.trees \
 -o analyses/phylogenetics/set103/trees/astral/astral10.tree
 
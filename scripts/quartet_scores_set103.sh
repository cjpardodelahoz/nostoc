#!/bin/bash

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=log/quartet_scores_set103.out
#SBATCH --error=log/quartet_scores_set103.err
#SBATCH --partition=scavenger

# Load newick utilities conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate newick
# Path for quartet scores
export PATH=/hpc/group/bio1/carlos/apps/:$PATH
# Directory for quartet score analyses
mkdir -p analyses/phylogenetics/set103/conflict/quartet_scores
# Put all trees into a single file
for locus in $(cat scripts/busco_ids_filtered_set103) ; do
 cat analyses/phylogenetics/set103/trees/single/${locus}*noempty.treefile >> \
 analyses/phylogenetics/set103/conflict/quartet_scores/ml_gene.trees
done
cat analyses/phylogenetics/set103/trees/single/16s*noempty.treefile >> \
 analyses/phylogenetics/set103/conflict/quartet_scores/ml_gene.trees
 cat analyses/phylogenetics/set103/trees/single/trnl*noempty.treefile >> \
 analyses/phylogenetics/set103/conflict/quartet_scores/ml_gene.trees
# Collapse nodes with less than 95% UFB. Has to run after running astral_set103.sh
nw_ed analyses/phylogenetics/set103/conflict/quartet_scores/ml_gene.trees 'i & b<=94' \
 o > analyses/phylogenetics/set103/conflict/quartet_scores/ml95_gene.trees
# Copy the weighted astral tree which will serve as a reference
cp analyses/phylogenetics/set103/trees/astral/wastral.tree \
 analyses/phylogenetics/set103/conflict/quartet_scores/wastral.tree
# Run quartet scores on all trees wihtout collapsing nodes
QuartetScores -t 16 \
 --ref analyses/phylogenetics/set103/conflict/quartet_scores/wastral.tree \
 --eval analyses/phylogenetics/set103/conflict/quartet_scores/ml_gene.trees \
 --output analyses/phylogenetics/set103/conflict/quartet_scores/quartet_scores_raw.tree
# Run quartet scores on trees
QuartetScores -t 16 \
 --ref analyses/phylogenetics/set103/conflict/quartet_scores/wastral.tree \
 --eval analyses/phylogenetics/set103/conflict/quartet_scores/ml95_gene.trees \
 --output analyses/phylogenetics/set103/conflict/quartet_scores/quartet_scores_95.tree

 
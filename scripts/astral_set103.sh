#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/astral_set103.out
#SBATCH --error=log/astral_set103.err
#SBATCH --partition=scavenger

module load Java/1.8.0_60
astral_path=/hpc/group/bio1/carlos/apps/Astral

mkdir -p analyses/phylogenetics/set103/trees/astral

# Put all trees into a single file
for locus in $(cat scripts/busco_ids_filtered_set103) ; do
 cat analyses/phylogenetics/set103/trees/single/${locus}*noempty.treefile >> \
 analyses/phylogenetics/set103/trees/astral/ml_gene.trees
done
cat analyses/phylogenetics/set103/trees/single/16s*noempty.treefile >> \
 analyses/phylogenetics/set103/trees/astral/ml_gene.trees
 cat analyses/phylogenetics/set103/trees/single/trnl*noempty.treefile >> \
 analyses/phylogenetics/set103/trees/astral/ml_gene.trees
# Run ASTRAL
java -jar ${astral_path}/astral.5.15.5.jar \
 -T 8 \
 -i analyses/phylogenetics/set103/trees/astral/ml_gene.trees \
 -o analyses/phylogenetics/set103/trees/astral/astral.tree
 
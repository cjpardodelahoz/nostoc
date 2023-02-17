#!/usr/bin/env Rscript

# Load required packages and functions
library(ape)

# Put all gene trees in a single file (the ones which have no empty ufboot)
system(command = "cat analyses/phylogenetics/set103/trees/single/*noempty.treefile > analyses/phylogenetics/set103/trees/single/all_noempty.trees")
# Load the trees
trees <- read.tree(file = "analyses/phylogenetics/set103/trees/single/all_noempty.trees") 
# Filter to the trees that have all taxa
trees_with_all_taxa <- list()
index <- 1
for (i in 1:length(trees)) {
  tree <- trees[[i]]
  if (length(tree$tip.label) == 151) {
    trees_with_all_taxa[[index]] <- tree
    index <- index + 1
  }
}
# Get majority rule consensus tree from the trees that have all taxa
mr_tree <- ape::consensus(trees_with_all_taxa, p = 0.5, check.labels = T)
# Save the mr tree
dir.create(path = "analyses/phylogenetics/set103/trees/mr")
write.tree(mr_tree, file = "analyses/phylogenetics/set103/trees/mr/mr.tree")
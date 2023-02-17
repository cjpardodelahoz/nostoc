#!/usr/bin/env Rscript

# This is a script to get the the biparts.temp file needed to generate the clade definitions file
# to run Discovista. When the clade_def file is generated with this function, the bipartitions will
# be in the same order as R assigns node labels when reading the tree, which allows connecting the 
# bipartitions with the node labels. Bipartition 1 will have all the taxa in the tree. 
# Bipartition b will correspond to internal node n+b in R, where n is the number of taxa in the tree.

# Usage:
# get_all_biparts.R treefile outfile outgroup_1 outgroup_2 outgroup_n
# 
# treefile      Tree to get the bipartitions.
# outfile       Path to the biparts.temp file 
# outgroupfile  File containing the Label(s) of the outgroup taxa. Each labe in a new line.

# Load required packages and functions
library(tidyverse)
library(treedataverse)
# Function to get a starter file with all bipartitions from a phylogenetic tree tree.
# The output can be used to generate the clade definition file required to run discovista
# tree    A rooted phylogenetic tree in newick format
get_biparts <- function(tree) {
  # Get the number of tips in the tree
  n_tips <- tree$tip.label %>% length()
  # Calculate the maximum number of internal nodes
  max_nodes <- n_tips - 1
  # Create the vector to store the bipartitions
  biparts <- character()
  # Iterate on all internal nodes of the tree
  for (i in 1:max_nodes) {
    node <- n_tips + i
    # Subset the tree to each internal node to get the taxa on one side of each bipartition
    tree_subset <- treeio::tree_subset(tree = tree, node = node, levels_back = 0)
    # Callapse taxa labels into a single character strin separated by commas
    biparts[i] <- tree_subset$tip.label %>%
      paste(collapse = ",")
  }
  return(biparts)
}

# Take variables from command line
args <- commandArgs(trailingOnly = T)
treefile <- args[1]
outfile <- args[2]
outgroupfile <- args[3]
# Load the outgroup file, the tree and root it
outgroup <- scan(file = outgroupfile, what = "character")
tree <- read.tree(file = treefile) %>%
  root.phylo(outgroup = outgroup, resolve.root = T)
# Get the bipartitions
biparts <- get_biparts(tree)
# Write the biparts.temp file
write(biparts, file = outfile)
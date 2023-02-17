#!/usr/bin/env Rscript

# Load required packages
library(ape)

# Read wastral tree and remove brnahc lengths and labels
tree <- read.tree(file = "analyses/phylogenetics/set103/trees/astral/wastral.tree")
tree$edge.length <- NULL
tree$node.label <- NULL
# Reroot
outgroup <- c("Cylindrospermum_stagnale_PCC_7417.fa", "Anabaena_cylindrica_PCC_7122.fa",
              "Aphanizomenon_flos_aquae_NIES_81.fa")
all_taxa <- tree$tip.label
outgroup <- setdiff(all_taxa, outgroup)
tree <- root.phylo(tree, outgroup = outgroup, resolve.root = T)
# Save tree
write.tree(tree, file = "analyses/phylogenetics/set103/trees/astral/wastral_divtime.phy")

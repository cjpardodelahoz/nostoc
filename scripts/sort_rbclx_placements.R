#!/usr/bin/env Rscript
library(tidyverse)
library(treeio)
library(ggtree)

outgroup <- c("Aphanizomenon_flos_aquae_NIES_81.fa",
             "Anabaena_cylindrica_PCC_7122.fa",
             "Cylindrospermum_stagnale_PCC_7417.fa")

tree <- read.tree("analyses/species_delimitation/cooccurrence/trees/placement/RAxML_labelledTree.test") %>%
  root(outgroup = outgroup, resolve.root = T)

MRCA(tree, c("JL34_bin_22.fa", "P2162_bin_10.fa"))

ggtree(tree_subset(tree, 2640, levels_back = 0)) + geom_tiplab()

df <- as_tibble(tree)

plot <- ggtree(tree, branch.length = "none") +
  geom_tiplab()

ggsave(plot, filename = "test.pdf", height = 19)

#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)
library(data.table)
library(ape)
library(treeio)
library(ggtree)

# Define outgroup taxa to root the tree
outgroup <- c("Aphanizomenon_flos_aquae_NIES_81.fa",
              "Anabaena_cylindrica_PCC_7122.fa",
              "Cylindrospermum_stagnale_PCC_7417.fa")
# Load tre with EPA placements
tree <- read.tree("analyses/species_delimitation/cooccurrence/trees/placement/RAxML_labelledTree.test") %>%
  root(outgroup = outgroup, resolve.root = T)
# MRCA node number for focal clades
v <- MRCA(tree, c("JL34_bin_22.fa", "P2162_bin_10.fa"))
# Subset tree to focal groups
tree_v <- tree_subset(tree, node = v, levels_back = 0)
# Get list of queries placed within focal groups
tree_v_df <- as_tibble(tree_v)
v_queries <- tree_v_df %>%
  filter(label %like% "QUERY__") %>%
  pull(label) %>%
  str_remove("QUERY___") %>%
  str_remove("_R_")
# Load edited ABMI sequences
abmi_seqs_all <- read.FASTA("analyses/species_delimitation/cooccurrence/seqs/rbclx_abmi_edited_all.fna")
# Subset ABMI seqs to queries from focal groups
abmi_seqs_v <- abmi_seqs_all[v_queries]
# Save FASTA with seqs within focal groups
write.FASTA(abmi_seqs_v, "analyses/species_delimitation/cooccurrence/seqs/rbclx_abmi_v.fna")




# Plot subset tree
ggtree(tree_v, layout = "circular")

ggtree(tree_subset(tree, 2640, levels_back = 0)) + geom_tiplab()

ggtree(tree) + geom_tiplab()
df <- as_tibble(tree)

plot <- ggtree(tree, branch.length = "none") +
  geom_tiplab()

ggsave(plot, filename = "test.pdf", height = 19)

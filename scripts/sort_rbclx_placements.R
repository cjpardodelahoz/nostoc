#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)
library(data.table)
library(ape)
library(treeio)

# Define outgroup taxa to root the tree
outgroup <- c("Aphanizomenon_flos_aquae_NIES_81.fa",
              "Anabaena_cylindrica_PCC_7122.fa",
              "Cylindrospermum_stagnale_PCC_7417.fa")
# Load tree with EPA placements
tree <- read.tree("analyses/species_delimitation/cooccurrence/trees/placement/RAxML_labelledTree.epa_result") %>%
  root(outgroup = outgroup, resolve.root = T)
# MRCA node number for focal clades
v <- MRCA(tree, c("JL34_bin_22.fa", "P2162_bin_10.fa"))
# Subset tree to focal groups
tree_v <- tree_subset(tree, node = v, levels_back = 0)
# Get list of taxa placed within focal groups
tree_v_df <- as_tibble(tree_v)
taxa_v <- tree_v_df %>%
  filter(!is.na(label)) %>%
  pull(label) %>%
  str_remove("QUERY___")
#queries_v <- tree_v_df %>%
#  filter(label %like% "QUERY__") %>%
#  pull(label) %>%
#  str_remove("QUERY___")
# Load all rbclx sequences
seqs_set103_abmi_global <- read.FASTA("analyses/species_delimitation/cooccurrence/seqs/rbclx_set103_abmi_global.fna")
# Subset rbclx seqs to taxa from focal groups
seqs_set103_abmi_global_v <- seqs_set103_abmi_global[taxa_v]
# Save FASTA with seqs within focal groups
write.FASTA(seqs_set103_abmi_global_v, 
            file = "analyses/species_delimitation/cooccurrence/seqs/rbclx_set103_global_abmi_v.fna")
# Save a newick copy of the tree with the queries placed
write.tree(tree, 
           file = "analyses/species_delimitation/cooccurrence/trees/placement/placement.tree")
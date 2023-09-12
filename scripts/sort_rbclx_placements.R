#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)
library(data.table)
library(ape)
library(treeio)
source("scripts/r_functions.R")

# Split placement tree into species complexes

# Define outgroup taxa to root the tree
outgroup <- c("Aphanizomenon_flos_aquae_NIES_81.fa",
              "Anabaena_cylindrica_PCC_7122.fa",
              "Cylindrospermum_stagnale_PCC_7417.fa",
              "QUERY___AJ632057", "QUERY___AJ293165", "QUERY___Z94888",
              "QUERY___DQ266032", "QUERY___DQ266030", "QUERY___DQ266029",
              "QUERY___DQ185301", "QUERY___DQ185297", "QUERY___DQ185264",
              "QUERY___AJ632066", "QUERY___P9373", "QUERY___P9367")
# Load tree with EPA placements
tree <- read.tree("analyses/species_delimitation/rbclx/clade_assignment/trees/placement/RAxML_labelledTree.epa_result") %>%
  root(outgroup = outgroup, resolve.root = T)
# Save a newick copy of the tree with the queries placed
write.tree(tree, 
           file = "analyses/species_delimitation/rbclx/clade_assignment/trees/placement/placement.tree")
# MRCA node number for focal clades
# Get species complex nodes
node_3_1 <- MRCA(tree, c("P2081_bin_11.fa", "P8768_bin_2.fa"))
node_3_2 <- MRCA(tree, c("P6447_bin_3.fa",  "Nostoc_sp_Lobaria_pulmonaria_5183.fa"))
node_3_3 <- MRCA(tree, c("NMS4_bin_2.fa", "QUERY___P6394"))
node_3_4 <- MRCA(tree, c("P6636_bin_11.fa", "P2213_bin_24.fa"))
node_3_5 <- MRCA(tree, c("S67_bin_3.fa", "P6602_bin_13.fa"))
node_3_6 <- MRCA(tree, c("P2162_bin_10.fa", "P8202_bin_9.fa"))
node_3_7 <- MRCA(tree, c("P6620_bin_4.fa", "P3068_bin_7.fa"))
node_3_8 <- MRCA(tree, c("Nostoc_KVS11.fa", "Nmoss2.fa"))
node_3_9 <- MRCA(tree, c("P943_bin_5.fa", "QUERY___P8820"))
node_3_10 <- MRCA(tree, c("P2039_bin_23.fa", "S32_bin_15.fa"))
node_3_11 <- MRCA(tree, c("NMS1_bin_5.fa", "P6963_bin_9.fa"))
node_3_12 <- MRCA(tree, c("P9820_bin_6.fa", "Nostoc_commune_NIES_4072.fa"))
node_2_1 <- MRCA(tree, c("P12545_bin_5.fa", "P12564_bin_10.fa"))
node_2_2 <- MRCA(tree, c("P12537_bin_5.fa", "P12591_bin_6.fa"))
node_2_3 <- MRCA(tree, c("P12588_bin_4.fa", "QUERY___EF102301"))
node_2_4 <- MRCA(tree, c("P10247_bin_16.fa", "NOS_bin_1.fa"))
node_subclade1 <- MRCA(tree, c("P12642_bin_1.fa", "QUERY___DQ185299"))
# Subset tree to focal groups
tree_3_1 <- tree_subset(tree, node = node_3_1, levels_back = 0)
tree_3_2 <- tree_subset(tree, node = node_3_2, levels_back = 0)
tree_3_3 <- tree_subset(tree, node = node_3_3, levels_back = 0)
tree_3_4 <- tree_subset(tree, node = node_3_4, levels_back = 0)
tree_3_5 <- tree_subset(tree, node = node_3_5, levels_back = 0)
tree_3_6 <- tree_subset(tree, node = node_3_6, levels_back = 0)
tree_3_7 <- tree_subset(tree, node = node_3_7, levels_back = 0)
tree_3_8 <- tree_subset(tree, node = node_3_8, levels_back = 0)
tree_3_9 <- tree_subset(tree, node = node_3_9, levels_back = 0)
tree_3_10 <- tree_subset(tree, node = node_3_10, levels_back = 0)
tree_3_11 <- tree_subset(tree, node = node_3_11, levels_back = 0)
tree_3_12 <- tree_subset(tree, node = node_3_12, levels_back = 0)
tree_2_1 <- tree_subset(tree, node = node_2_1, levels_back = 0)
tree_2_2 <- tree_subset(tree, node = node_2_2, levels_back = 0)
tree_2_3 <- tree_subset(tree, node = node_2_3, levels_back = 0)
tree_2_4 <- tree_subset(tree, node = node_2_4, levels_back = 0)
tree_subclade1 <- tree_subset(tree, node = node_subclade1, levels_back = 0)

# Extract tree with unassigned queries

# Make table to indicate where queries come from (ABMI or public)
ref_taxa <- scan(file = "misc_files/genome_ids_set103", what = "character")
query_source_df <- tree$tip.label %>%
  as_tibble() %>%
  rename(taxon = value) %>%
  mutate(source = case_when(taxon %in% ref_taxa ~ "reference_tree",
                            str_detect(taxon, "QUERY___P[1-9]*") ~ "abmi",
                            .default = "public"))
# All Nostoc queries
all_queries <- query_source_df %>% filter(source == "public" | source == "abmi") %>%
  filter(taxon %nin% outgroup) %>%
  pull(taxon)
# queries placed within predefined focal groups
assigned_queries <- c(tree_3_1$tip.label, tree_3_2$tip.label, tree_3_3$tip.label,
                      tree_3_4$tip.label, tree_3_5$tip.label, tree_3_6$tip.label,
                      tree_3_7$tip.label, tree_3_8$tip.label, tree_3_9$tip.label,
                      tree_3_10$tip.label, tree_3_11$tip.label, tree_3_12$tip.label,
                      tree_2_1$tip.label, tree_2_2$tip.label, tree_2_3$tip.label,
                      tree_2_4$tip.label, tree_subclade1$tip.label) %>%
  as_tibble() %>%
  filter(value %nin% ref_taxa) %>%
  pull(value)
# queries that fell outside of focal groups
unassigned_queries <- setdiff(all_queries, assigned_queries)
write(unassigned_queries, file = "analyses/species_delimitation/rbclx/clade_assignment/trees/placement/unassigned_queries.txt")
# Tree showing only unassigned queries
tree_unassigned <- keep.tip(tree, c(ref_taxa, unassigned_queries))
write.tree(tree_unassigned, file = "analyses/species_delimitation/rbclx/clade_assignment/trees/placement/unassigned.tree")

# Sort sequences of assigned queries

# Sort query labels from focal groups
tree_3_1_labels <- tree_3_1$tip.label %>% str_remove("QUERY___")
tree_3_2_labels <- tree_3_2$tip.label %>% str_remove("QUERY___")
tree_3_3_labels <- tree_3_3$tip.label %>% str_remove("QUERY___")
tree_3_4_labels <- tree_3_4$tip.label %>% str_remove("QUERY___")
tree_3_5_labels <- tree_3_5$tip.label %>% str_remove("QUERY___")
tree_3_6_labels <- tree_3_6$tip.label %>% str_remove("QUERY___")
tree_3_7_labels <- tree_3_7$tip.label %>% str_remove("QUERY___")
tree_3_8_labels <- tree_3_8$tip.label %>% str_remove("QUERY___")
tree_3_9_labels <- tree_3_9$tip.label %>% str_remove("QUERY___")
tree_3_10_labels <- tree_3_10$tip.label %>% str_remove("QUERY___")
tree_3_11_labels <- tree_3_11$tip.label %>% str_remove("QUERY___")
tree_3_12_labels <- tree_3_12$tip.label %>% str_remove("QUERY___")
tree_2_1_labels <- tree_2_1$tip.label %>% str_remove("QUERY___")
tree_2_2_labels <- tree_2_2$tip.label %>% str_remove("QUERY___")
tree_2_3_labels <- tree_2_3$tip.label %>% str_remove("QUERY___")
tree_2_4_labels <- tree_2_4$tip.label %>% str_remove("QUERY___")
tree_subclade1_labels <- tree_subclade1$tip.label %>% str_remove("QUERY___")
# Load all rbclx sequences
seqs_set103_abmi_public <- read.FASTA("analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_set103_abmi_public.fna")
# Subset rbclx seqs to taxa from focal groups
seqs_3_1 <- seqs_set103_abmi_public[tree_3_1_labels]
seqs_3_2 <- seqs_set103_abmi_public[tree_3_2_labels]
seqs_3_3 <- seqs_set103_abmi_public[tree_3_3_labels]
seqs_3_4 <- seqs_set103_abmi_public[tree_3_4_labels]
seqs_3_5 <- seqs_set103_abmi_public[tree_3_5_labels]
seqs_3_6 <- seqs_set103_abmi_public[tree_3_6_labels]
seqs_3_7 <- seqs_set103_abmi_public[tree_3_7_labels]
seqs_3_8 <- seqs_set103_abmi_public[tree_3_8_labels]
seqs_3_9 <- seqs_set103_abmi_public[tree_3_9_labels]
seqs_3_10 <- seqs_set103_abmi_public[tree_3_10_labels]
seqs_3_11 <- seqs_set103_abmi_public[tree_3_11_labels]
seqs_3_12 <- seqs_set103_abmi_public[tree_3_12_labels]
seqs_2_1 <- seqs_set103_abmi_public[tree_2_1_labels]
seqs_2_2 <- seqs_set103_abmi_public[tree_2_2_labels]
seqs_2_3 <- seqs_set103_abmi_public[tree_2_3_labels]
seqs_2_4 <- seqs_set103_abmi_public[tree_2_4_labels]
seqs_subclade1 <- seqs_set103_abmi_public[tree_subclade1_labels]
# Save FASTA with seqs within focal groups
write.FASTA(seqs_3_1, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_3_1.fna")
write.FASTA(seqs_3_2, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_3_2.fna")
write.FASTA(seqs_3_3, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_3_3.fna")
write.FASTA(seqs_3_4, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_3_4.fna")
write.FASTA(seqs_3_5, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_3_5.fna")
write.FASTA(seqs_3_6, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_3_6.fna")
write.FASTA(seqs_3_7, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_3_7.fna")
write.FASTA(seqs_3_8, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_3_8.fna")
write.FASTA(seqs_3_9, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_3_9.fna")
write.FASTA(seqs_3_10, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_3_10.fna")
write.FASTA(seqs_3_11, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_3_11.fna")
write.FASTA(seqs_3_12, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_3_12.fna")
write.FASTA(seqs_2_1, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_2_1.fna")
write.FASTA(seqs_2_2, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_2_2.fna")
write.FASTA(seqs_2_3, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_2_3.fna")
write.FASTA(seqs_2_4, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_2_4.fna")
write.FASTA(seqs_3_12, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_3_12.fna")
write.FASTA(seqs_subclade1, file = "analyses/species_delimitation/rbclx/clade_assignment/seqs/rbclx_subclade1.fna")

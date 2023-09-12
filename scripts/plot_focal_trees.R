#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)
library(data.table)
library(ape)
library(phytools)
library(treeio)
library(ggtree)
source("scripts/r_functions.R")


#### PLOT TREES WITH PUBLIC METADATA ####

# Load ANI and popcogent clusters
ani95_clusters <- read_csv(file = "analyses/species_delimitation/fastani/set12c/ani_95_clusters.csv") %>%
  mutate(ani_95_cluster = str_replace(ani_95_cluster, "c", "ANI-"))
popcogent_clusters <- read_delim(file = "analyses/species_delimitation/popcogent/set12c/set12c_0.000355362.txt.cluster.tab.txt",
                                 delim = "\t") %>%
  select(Strain, Main_cluster) %>%
  mutate(Main_cluster =
           Main_cluster + 1) %>%
  mutate(Strain = 
           str_replace(Strain, "_chromosome", ".fa")) %>%
  mutate(Main_cluster = paste("PopCO-", Main_cluster, sep = ""))
# Load metadata for public seqs and update the GIs with their corresponding
# accessions
unmatched_gis <- scan(file = "analyses/species_delimitation/rbclx/public/unmatched_gis.txt", what = "character")
unmatched_accessions <- scan(file = "analyses/species_delimitation/rbclx/public/unmatched_accessions.txt", what = "character")
gi_accession_key <- tibble(gi = unmatched_gis, accession = unmatched_accessions)
public_rbclx_metadata <- read_csv(file = "analyses/species_delimitation/rbclx/public/public_rbclx_metadata.csv") %>% 
  left_join(gi_accession_key, by = c("rbclx_accession" = "gi")) %>%
  mutate(rbclx_accession = if_else(is.na(accession),
                                   rbclx_accession, accession)) %>%
  select(-accession) 
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
# Make table to indicate where queries come from (ABMI or public)
ref_taxa <- scan(file = "misc_files/genome_ids_set103", what = "character")
query_source_df <- tree$tip.label %>%
  as_tibble() %>%
  rename(taxon = value) %>%
  mutate(source = case_when(taxon %in% ref_taxa ~ "reference_tree",
                            str_detect(taxon, "QUERY___P[1-9]*") ~ "abmi",
                            .default = "public")) %>%
  mutate(taxon = str_remove(taxon, "QUERY___"))
# Reference and public taxa
ref_public_taxa <- query_source_df %>%
  filter(source == "public" | source == "reference_tree") %>%
  pull(taxon)
# Load table with ref taxon name key
ref_key_df <- read_csv("document/tables/voucher_v1.csv") %>%
  select(genome_id, taxon_name, country) %>%
  mutate(country = str_replace(country, "United Kingdom", "UK")) %>%
  mutate(country = str_replace(country, "United States of America", "USA")) %>%
  mutate(country = str_replace(country, "United States of America", "Panama")) %>%
  mutate(taxon_name =
           ifelse(taxon_name == genome_id,
                  taxon_name,
                  paste(genome_id, taxon_name, sep = "_"))) %>%
  mutate(genome_id = 
           paste(genome_id, ".fa", sep = "")) %>%
  mutate(taxon_name = 
           ifelse(!is.na(country),
                  paste(taxon_name, country, sep = "_"),
                  taxon_name)) %>%
  select(genome_id, taxon_name)
# Get df key to replace tip labes
tiplabel_key <- query_source_df %>%
  filter(source == "public" | source == "reference_tree") %>%
  left_join(public_rbclx_metadata, by = c("taxon" = "rbclx_accession")) %>%
  left_join(ani95_clusters, by = c("taxon" = "genome")) %>%
  left_join(popcogent_clusters, by = c("taxon" = "Strain")) %>%
  left_join(ref_key_df, by = c("taxon" = "genome_id")) %>%
  mutate(new_label = if_else(source == "public",
                             paste(taxon, dna_source, region, "Phylogroup", nostoc_phylogroup),
                             paste(taxon_name, ani_95_cluster, Main_cluster))) %>%
  select(taxon, new_label)
# Function to prepare data to plot rbclx trees with public metadata
prep_public_focal_trees <- function(treefile, outgroup_taxa, font_size,
                                    col_offset) {
  tree <- read.tree(file = treefile) %>%
    root(outgroup = outgroup_taxa, resolve.root = T)
  tree_taxa <- intersect(tree$tip.label, ref_public_taxa)
  tree_tiplabel_key <- tiplabel_key %>%
    filter(taxon %in% tree_taxa)
  tree <- keep.tip(tree, tree_taxa) %>%
    rename_taxa(tree_tiplabel_key)
  toy_df <- tree_tiplabel_key %>% distinct() %>%
    column_to_rownames(var = "new_label") %>%
    mutate(taxon = "not_a_variable")
  tree_plot <- ggtree(tree) +
    geom_tiplab(align = T, size = font_size)
  tree_plot <- gheatmap(tree_plot, toy_df, 
                        width = 0.1, offset = col_offset, colnames = F) +
    theme(legend.position = "none")
}
# Clade 3-1
outgroup_3_1 <- c("P2081_bin_11.fa", "MH770776")
clade_3_1_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_1.treefile",
  outgroup_taxa = outgroup_3_1, font_size = 3, col_offset = 0.1) +
  geom_treescale(width = 0.005, y = 50, x = 0.003, offset = 1.5)
ggsave(clade_3_1_public_plot, file = "document/plots/clade_3_1_public.pdf",
       units = "cm", height = 40, width = 40)
# Clade 3-2
outgroup_3_2 <- c("NMS9_bin_7.fa", "S44_bin_8.fa", "KX923109")
clade_3_2_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_2.treefile",
  outgroup_taxa = outgroup_3_2, font_size = 3, col_offset = 0.1) +
  geom_treescale(width = 0.005, y = 20, x = 0.003, offset = 0.5)
ggsave(clade_3_2_public_plot, file = "document/plots/clade_3_2_public.pdf",
       units = "cm", height = 20, width = 30)
 # Clade 3-4
outgroup_3_4 <- c("KX923097")
clade_3_4_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_4.treefile",
  outgroup_taxa = outgroup_3_4, font_size = 3, col_offset = 0.1) +
  geom_treescale(width = 0.005, y = 30, x = 0.003, offset = 1.5)
ggsave(clade_3_4_public_plot, file = "document/plots/clade_3_4_public.pdf")
# Clade 3-5
outgroup_3_5 <- c("DQ185281", "KX922920", "KX922926", "KX922928", "S67_bin_3.fa",
                  "KX922927","KX922930", "KX922971", "KX923043", 
                  "P3034_bin_5.fa", "S66_bin_1.fa")
clade_3_5_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_5_with_spacer.treefile",
  outgroup_taxa = outgroup_3_5, font_size = 3, col_offset = 0.08) +
  geom_treescale(width = 0.005, y = 30, x = 0.003, offset = 1.5)
ggsave(clade_3_5_public_plot, file = "document/plots/clade_3_5_public.pdf")
# Clade 3-6
tmp_tree <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_6.treefile")
#root_node_3_6 <- MRCA(tmp_tree, c("P8202_bin_9.fa", "P8569_bin_1.fa")) 
root_node_3_6 <- MRCA(tmp_tree, c("P8202_bin_9.fa", "P8569_bin_1.fa")) 
outgroup_3_6 <- tree_subset(tmp_tree, node = root_node_3_6, levels_back = 0)$tip.label
clade_3_6_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_6.treefile",
  outgroup_taxa = outgroup_3_6, font_size = 3, col_offset = 0.08) +
  geom_treescale(width = 0.005, y = 30, x = 0.003, offset = 1.5)
ggsave(clade_3_6_public_plot, file = "document/plots/clade_3_6_public.pdf",
       units = "cm", height = 40, width = 40)
# Clade 3-7
tmp_tree <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_7.treefile")
root_node_3_7 <- MRCA(tmp_tree, c("Z94893", "P3068_bin_7.fa")) 
outgroup_3_7 <- tree_subset(tmp_tree, node = root_node_3_7, levels_back = 0)$tip.label
clade_3_7_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_7.treefile",
  outgroup_taxa = outgroup_3_7, font_size = 3, col_offset = 0.08) +
  geom_treescale(width = 0.005, y = 20, x = 0.003, offset = 1.0)
ggsave(clade_3_7_public_plot, file = "document/plots/clade_3_7_public.pdf")
# Clade 3-8
outgroup_3_8 <- c("Nostoc_KVJ2.fa", "Nostoc_KVS11.fa")
clade_3_8_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_8.treefile",
  outgroup_taxa = outgroup_3_8, font_size = 3, col_offset = 0.08) +
  geom_treescale(width = 0.005, y = 5, x = 0.003, offset = 0.2)
ggsave(clade_3_8_public_plot, file = "document/plots/clade_3_8_public.pdf")
# Clade 3-9
outgroup_3_9 <- c("KX922962", "KX922961")
clade_3_9_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_9.treefile",
  outgroup_taxa = outgroup_3_9, font_size = 3, col_offset = 0.08) +
  geom_treescale(width = 0.005, y = 5, x = 0.003, offset = .2)
ggsave(clade_3_9_public_plot, file = "document/plots/clade_3_9_public.pdf")
# Clade 3-10
tmp_tree <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_10.treefile")
root_node_3_10 <- MRCA(tmp_tree, c("P2037_bin_28.fa", "P2039_bin_23.fa")) 
outgroup_3_10 <- tree_subset(tmp_tree, node = root_node_3_10, levels_back = 0)$tip.label
clade_3_10_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_10.treefile",
  outgroup_taxa = outgroup_3_10, font_size = 3, col_offset = 0.08) +
  geom_treescale(width = 0.005, y = 20, x = 0.003, offset = 1)
ggsave(clade_3_10_public_plot, file = "document/plots/clade_3_10_public.pdf",
       units = "cm", height = 40, width = 40)
# Clade 3-11
outgroup_3_11 <- c("NMS1_bin_5.fa", "NMS2_bin_6.fa")
clade_3_11_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_11.treefile",
  outgroup_taxa = outgroup_3_11, font_size = 3, col_offset = 0.09) +
  geom_treescale(width = 0.005, y = 10, x = 0.003, offset = .2)
ggsave(clade_3_11_public_plot, file = "document/plots/clade_3_11_public.pdf")
# Clade 3-12
outgroup_3_12 <- c("P9820_bin_6.fa")
clade_3_12_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_12.treefile",
  outgroup_taxa = outgroup_3_12, font_size = 3, col_offset = 0.09) +
  geom_treescale(width = 0.005, y = 25, x = 0.003, offset = .2)
ggsave(clade_3_12_public_plot, file = "document/plots/clade_3_12_public.pdf")
# Clade 2-1
outgroup_2_1 <- c("EF102347")
clade_2_1_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_2_1.treefile",
  outgroup_taxa = outgroup_2_1, font_size = 3, col_offset = 0.09) +
  geom_treescale(width = 0.005, y = 30, x = 0.003, offset = .5)
ggsave(clade_2_1_public_plot, file = "document/plots/clade_2_1_public.pdf")
# Clade 2-2
outgroup_2_2 <- c("S8_bin_2.fa", "S9_bin_6.fa", "KX923058", "EU877488", "P12591_bin_6.fa")
clade_2_2_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_2_2.treefile",
  outgroup_taxa = outgroup_2_2, font_size = 3, col_offset = 0.09) +
  geom_treescale(width = 0.005, y = 5, x = 0.003, offset = .2)
ggsave(clade_2_2_public_plot, file = "document/plots/clade_2_2_public.pdf")
# Clade 2-4
tmp_tree <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_2_4.treefile")
root_node_2_4 <- MRCA(tmp_tree, c("P8690_bin_1.fa", "JL33_bin_16.fa", "P10247_bin_16.fa")) 
outgroup_2_4 <- tree_subset(tmp_tree, node = root_node_2_4, levels_back = 0)$tip.label
clade_2_4_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_2_4.treefile",
  outgroup_taxa = outgroup_2_4, font_size = 3, col_offset = 0.09) +
  geom_treescale(width = 0.005, y = 50, x = 0.003, offset = 1)
ggsave(clade_2_4_public_plot, file = "document/plots/clade_2_4_public.pdf",
       units = "cm", height = 30, width = 30)


#### GENERATE TABLES WITH REVISED RBCLX CLASSIFICATION ####

# Section 3-1

# Load and root section tree the same way as we plotted before
outgroup_3_1 <- c("P2081_bin_11.fa", "MH770776")
tree_3_1 <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_1.treefile") %>%
  root(outgroup = outgroup_3_1, resolve.root = T)
# Classify phylogroup XXXVII taxa
phylogroup_xxxvii_node <- MRCA(tree_3_1, c("P2081_bin_11.fa", "MH770776"))
phylogroup_xxxvii_taxa <- tree_subset(tree_3_1, node = phylogroup_xxxvii_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XXXVII") %>%
  rename(dna_id = value)
# Classify phylogroup XXXIV taxa
phylogroup_xxxiv_node <- MRCA(tree_3_1, c("MH770581", "MH770583"))
phylogroup_xxxiv_taxa <- tree_subset(tree_3_1, node = phylogroup_xxxiv_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XXXIV") %>%
  rename(dna_id = value)
# Classify phylogroup XXXVIII taxa
phylogroup_xxxviii_node <- MRCA(tree_3_1, c("MH770775", "DQ185283"))
phylogroup_xxxviii_taxa <- tree_subset(tree_3_1, node = phylogroup_xxxviii_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XXXVIII") %>%
  rename(dna_id = value)
# Classify phylogroup XXXIII taxa
phylogroup_xxxiii_node <- MRCA(tree_3_1, c("P8575_bin_3.fa", "JL31_bin_38.fa"))
phylogroup_xxxiii_taxa <- tree_subset(tree_3_1, node = phylogroup_xxxiii_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XXXIII") %>%
  rename(dna_id = value)
# Classify phylogroup XXIX taxa
phylogroup_xxix_node <- MRCA(tree_3_1, c("MH770767", "EU877530"))
phylogroup_xxix_taxa <- tree_subset(tree_3_1, node = phylogroup_xxix_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XXIX") %>%
  rename(dna_id = value)
# Classify phylogroup XXVII taxa
phylogroup_xxvii_node <- MRCA(tree_3_1, c("MH770576", "DQ185308"))
phylogroup_xxvii_taxa <- tree_subset(tree_3_1, node = phylogroup_xxvii_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XXVII") %>%
  rename(dna_id = value)
# Classify phylogroup XXV taxa
phylogroup_xxv_node <- MRCA(tree_3_1, c("MH770616", "MH770769"))
phylogroup_xxv_taxa <- tree_subset(tree_3_1, node = phylogroup_xxv_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XXV") %>%
  rename(dna_id = value)
# Classify species complex 3-1a
complex_3_1_node <- MRCA(tree_3_1, c("KX922914", "MH770769"))
complex_3_1_taxa <- tree_subset(tree_3_1, node = complex_3_1_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(species_complex = "3.1a") %>%
  rename(dna_id = value)
# Classify section 3-1
section_3_1_taxa <- tree_3_1$tip.label %>%
  as_tibble() %>%
  mutate(section = "3.1") %>%
  rename(dna_id = value)
# Join all assignments
clade_assignments_3_1 <- bind_rows(phylogroup_xxxvii_taxa, phylogroup_xxxiv_taxa, phylogroup_xxxviii_taxa,
                                   phylogroup_xxxiii_taxa, phylogroup_xxix_taxa, phylogroup_xxvii_taxa, 
                                   phylogroup_xxv_taxa) %>%
  right_join(section_3_1_taxa, by = "dna_id") %>%
  left_join(complex_3_1_taxa, by = "dna_id")

# Section 3-2

# Load and root section tree the same way as we plotted before
outgroup_3_2 <- c("NMS9_bin_7.fa", "S44_bin_8.fa", "KX923109")
tree_3_2 <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_2.treefile") %>%
  root(outgroup = outgroup_3_2, resolve.root = T)
# Classify phylogroup NEW1 taxa
phylogroup_new1_node <- MRCA(tree_3_2, c("P11060_bin_1.fa", "EU877528"))
phylogroup_new1_taxa <- tree_subset(tree_3_2, node = phylogroup_new1_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "NEW1") %>%
  rename(dna_id = value)
# Classify phylogroup NEW2 taxa
phylogroup_new2_node <- MRCA(tree_3_2, c("S44_bin_8.fa", "KX923109"))
phylogroup_new2_taxa <- tree_subset(tree_3_2, node = phylogroup_new2_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "NEW2") %>%
  rename(dna_id = value)
# Classify section 3-2
section_3_2_taxa <- tree_3_2$tip.label %>%
  as_tibble() %>%
  mutate(section = "3.2") %>%
  rename(dna_id = value) %>%
  add_column(species_complex = "3.2a")
# Join all assignments
clade_assignments_3_2 <- bind_rows(phylogroup_new1_taxa, phylogroup_new2_taxa) %>%
  right_join(section_3_2_taxa, by = "dna_id")

# Section 3-3

# Only one genome in the section, and sequences that were placed with can be 94%
# similar so, will consider them aff_3.3
clade_assignments_3_3 <- read.FASTA(file = "analyses/species_delimitation/rbclx/clade_assignment/alignments/rbclx_3_3_aln.fna") %>%
  names() %>%
  as_tibble() %>%
  rename(dna_id = value) %>%
  mutate(section = ifelse(dna_id == "NMS4_bin_2.fa",
                          "3.3", "aff_3.3")) %>%
  add_column(species_complex = NA) %>%
  add_column(phylogroup = NA)

# Section 3-4

# This is all phylogroup VI
clade_assignments_3_4 <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_4.treefile")$tip.label %>%
  as_tibble() %>%
  rename(dna_id = value) %>%
  mutate(section = "3.4") %>%
  add_column(species_complex = NA) %>%
  add_column(phylogroup = "VI")

# Section 3-5
  
# Load and root section tree the same way as we plotted before
outgroup_3_5 <- c("DQ185281", "KX922920", "KX922926", "KX922928", "S67_bin_3.fa",
                  "KX922927","KX922930", "KX922971", "KX923043", 
                  "P3034_bin_5.fa", "S66_bin_1.fa")
tree_3_5 <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_5_with_spacer.treefile") %>%
  root(outgroup = outgroup_3_5, resolve.root = T)
# Classify phylogroup VIIb taxa
phylogroup_viib_node <- MRCA(tree_3_5, c("S66_bin_1.fa", "KX922926"))
phylogroup_viib_taxa <- tree_subset(tree_3_5, node = phylogroup_viib_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "VIIb") %>%
  rename(dna_id = value)
# Classify phylogroup VIId taxa
phylogroup_viid_node <- MRCA(tree_3_5, c("P10089", "KX923001"))
phylogroup_viid_taxa <- tree_subset(tree_3_5, node = phylogroup_viid_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "VIId") %>%
  rename(dna_id = value)
# Classify phylogroup VIIa taxa
phylogroup_viia_node <- MRCA(tree_3_5, c("P6455", "KX923000"))
phylogroup_viia_taxa <- tree_subset(tree_3_5, node = phylogroup_viia_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "VIIa") %>%
  rename(dna_id = value)
# Classify phylogroup VIIc taxa REVISE THIS BECAUSE VIIc DOESNT WORK
phylogroup_viic_node <- MRCA(tree_3_5, c("P6465_bin_9.fa", "S43_bin_7.fa"))
phylogroup_viic_taxa <- tree_subset(tree_3_5, node = phylogroup_viic_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "VIIc") %>%
  rename(dna_id = value)
# Classify phylogroup Xa taxa
phylogroup_xa_node <- MRCA(tree_3_5, c("KX922934", "KX923014"))
phylogroup_xa_taxa <- tree_subset(tree_3_5, node = phylogroup_xa_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "Xa") %>%
  rename(dna_id = value)
# Classify phylogroup Xb taxa
phylogroup_xb_node <- MRCA(tree_3_5, c("KX922933", "KX922935"))
phylogroup_xb_taxa <- tree_subset(tree_3_5, node = phylogroup_xb_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "Xb") %>%
  rename(dna_id = value)
# Classify phylogroup XV taxa
phylogroup_xv_node <- MRCA(tree_3_5, c("KX923032", "KX922990"))
phylogroup_xv_taxa <- tree_subset(tree_3_5, node = phylogroup_xv_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XV") %>%
  rename(dna_id = value)
# Classify phylogroup XXVIa taxa
phylogroup_xxvia_node <- MRCA(tree_3_5, c("MH770562", "MH770561"))
phylogroup_xxvia_taxa <- tree_subset(tree_3_5, node = phylogroup_xxvia_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XXVIa") %>%
  rename(dna_id = value)
# Classify species complex 3-5a
complex_3_5_node <- MRCA(tree_3_5, c("KX923001", "MH770561"))
complex_3_5_taxa <- tree_subset(tree_3_5, node = complex_3_5_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(species_complex = "3.5a") %>%
  rename(dna_id = value)
# Classify section 3-5
section_3_5_taxa <- tree_3_5$tip.label %>%
  as_tibble() %>%
  mutate(section = "3.5") %>%
  rename(dna_id = value)
# Join all assignments
clade_assignments_3_5 <- bind_rows(phylogroup_viia_taxa, phylogroup_viib_taxa, phylogroup_viic_taxa, 
                                   phylogroup_viid_taxa, phylogroup_xa_taxa, phylogroup_xb_taxa, 
                                   phylogroup_xv_taxa, phylogroup_xxvia_taxa) %>%
  right_join(section_3_5_taxa, by = "dna_id") %>%
  left_join(complex_3_5_taxa, by = "dna_id")


# Section 3-6

# Load and root section tree the same way as we plotted before
tmp_tree <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_6.treefile")
root_node_3_6 <- MRCA(tmp_tree, c("P8202_bin_9.fa", "P8569_bin_1.fa")) 
outgroup_3_6 <- tree_subset(tmp_tree, node = root_node_3_6, levels_back = 0)$tip.label
tree_3_6 <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_6.treefile") %>%
  root(outgroup = outgroup_3_6, resolve.root = T)
# Classify phylogroup V taxa
phylogroup_v_node <- MRCA(tree_3_6, c("P6572", "KX922985"))
phylogroup_v_taxa <- tree_subset(tree_3_6, node = phylogroup_v_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "V") %>%
  rename(dna_id = value)
# Classify phylogroup XLII taxa
phylogroup_xlii_node <- MRCA(tree_3_6, c("MH770701", "P8202_bin_9.fa"))
phylogroup_xlii_taxa <- tree_subset(tree_3_6, node = phylogroup_xlii_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XLII") %>%
  rename(dna_id = value)
# Classify phylogroup XLII taxa
phylogroup_xlii_node <- MRCA(tree_3_6, c("MH770701", "P8202_bin_9.fa"))
phylogroup_xlii_taxa <- tree_subset(tree_3_6, node = phylogroup_xlii_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XLII") %>%
  rename(dna_id = value)
# Classify species complex 3.6a taxa
complex_3_6_node <- MRCA(tree_3_6, c("KX923047", "MH770745"))
complex_3_6_taxa <- tree_subset(tree_3_6, node = complex_3_6_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(species_complex = "3.6a") %>%
  rename(dna_id = value)
# Classify section 3-6
section_3_6_taxa <- tree_3_6$tip.label %>%
  as_tibble() %>%
  mutate(section = "3.6") %>%
  rename(dna_id = value)
# Join all assignments
clade_assignments_3_6 <- bind_rows(phylogroup_v_taxa, phylogroup_xlii_taxa) %>%
  right_join(section_3_6_taxa, by = "dna_id") %>%
  left_join(complex_3_6_taxa, by = "dna_id")

# Section 3-7

# Load and root section tree the same way as we plotted before
tmp_tree <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_7.treefile")
root_node_3_7 <- MRCA(tmp_tree, c("Z94893", "P3068_bin_7.fa")) 
outgroup_3_7 <- tree_subset(tmp_tree, node = root_node_3_7, levels_back = 0)$tip.label
tree_3_7 <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_7.treefile") %>%
  root(outgroup = outgroup_3_7, resolve.root = T)
# Classify phylogroup XIa taxa
phylogroup_xia_node <- MRCA(tree_3_7, c("P5023_bin_4.fa", "KX923052"))
phylogroup_xia_taxa <- tree_subset(tree_3_7, node = phylogroup_xia_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XIa") %>%
  rename(dna_id = value)
# Classify phylogroup XIb taxa
phylogroup_xib_node <- MRCA(tree_3_7, c("P3068_bin_7.fa", "KX922955"))
phylogroup_xib_taxa <- tree_subset(tree_3_7, node = phylogroup_xib_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XIb") %>%
  rename(dna_id = value)
# Classify section 3-7
section_3_7_taxa <- tree_3_7$tip.label %>%
  as_tibble() %>%
  mutate(section = "3.7") %>%
  rename(dna_id = value)
# Join all assignments
clade_assignments_3_7 <- bind_rows(phylogroup_xia_taxa, phylogroup_xib_taxa) %>%
  right_join(section_3_7_taxa, by = "dna_id") %>%
  add_column(species_complex = NA)

# Section 3-8

# Load and root section tree the same way as we plotted before
outgroup_3_8 <- c("Nostoc_KVJ2.fa", "Nostoc_KVS11.fa")
tree_3_8 <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_8.treefile") %>%
  root(outgroup = outgroup_3_8, resolve.root = T)
# Classify phylogroup NEW3 taxa
phylogroup_new3_node <- MRCA(tree_3_8,  c("Nostoc_KVJ2.fa", "Nostoc_KVS11.fa"))
phylogroup_new3_taxa <- tree_subset(tree_3_8, node = phylogroup_new3_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "NEW3") %>%
  rename(dna_id = value)
# Classify phylogroup NEW4 taxa
phylogroup_new4_taxa <- drop.tip(tree_3_8, outgroup_3_8)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "NEW4") %>%
  rename(dna_id = value)
# Join all assignments
clade_assignments_3_8 <- bind_rows(phylogroup_new3_taxa, phylogroup_new4_taxa) %>%
  add_column(section = "3.8") %>%
  add_column(species_complex = NA)

# Section 3-9

# Load and root section tree the same way as we plotted before
outgroup_3_9 <- c("KX922962", "KX922961")
tree_3_9 <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_9.treefile") %>%
  root(outgroup = outgroup_3_9, resolve.root = T)
# Classify phylogroup XXI taxa
phylogroup_xxi_node <- MRCA(tree_3_9,  c("MH770584", "P8820"))
phylogroup_xxi_taxa <- tree_subset(tree_3_9, node = phylogroup_xxi_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XXI") %>%
  rename(dna_id = value) %>%
  add_column(section = "3.9")
# Classify phylogroup XXI taxa
phylogroup_xvii_node <- MRCA(tree_3_9,  c("KX922962", "KX922961"))
phylogroup_xvii_taxa <- tree_subset(tree_3_9, node = phylogroup_xvii_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XVII") %>%
  rename(dna_id = value) %>%
  add_column(section = "aff_3.9")
# Join all assignments
clade_assignments_3_9 <- bind_rows(phylogroup_xxi_taxa, phylogroup_xvii_taxa) %>%
  add_column(species_complex = NA)

# Section 3-10

# Load and root section tree the same way as we plotted before
tmp_tree <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_10.treefile")
root_node_3_10 <- MRCA(tmp_tree, c("P2037_bin_28.fa", "P2039_bin_23.fa")) 
outgroup_3_10 <- tree_subset(tmp_tree, node = root_node_3_10, levels_back = 0)$tip.label
tree_3_10 <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_10.treefile") %>%
  root(outgroup = outgroup_3_10, resolve.root = T)
# Classify phylogroup XXII taxa
phylogroup_xxii_node <- MRCA(tree_3_10,  c("P2037_bin_28.fa", "MH770512"))
phylogroup_xxii_taxa <- tree_subset(tree_3_10, node = phylogroup_xxii_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XXII") %>%
  rename(dna_id = value)
# Classify phylogroup XXIII taxa
phylogroup_xxiii_node <- MRCA(tree_3_10,  c("P2039_bin_23.fa", "MH770527"))
phylogroup_xxiii_taxa <- tree_subset(tree_3_10, node = phylogroup_xxiii_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XXIII") %>%
  rename(dna_id = value)
# Classify phylogroup XIII.XLIII taxa
phylogroup_xiii_xliii_node <- MRCA(tree_3_10,  c("MK720761", "KX922938"))
phylogroup_xiii_xliii_taxa <- tree_subset(tree_3_10, node = phylogroup_xiii_xliii_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XIII.XLIII") %>%
  rename(dna_id = value)
# Classify phylogroup XIII.XLIII taxa
phylogroup_xiii_xliii_node <- MRCA(tree_3_10,  c("MK720761", "KX922938"))
phylogroup_xiii_xliii_taxa <- tree_subset(tree_3_10, node = phylogroup_xiii_xliii_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XIII.XLIII") %>%
  rename(dna_id = value)
# Classify phylogroup VIIIa taxa
phylogroup_viiia_node <- MRCA(tree_3_10,  c("MH770511", "KX923074"))
phylogroup_viiia_taxa <- tree_subset(tree_3_10, node = phylogroup_viiia_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "VIIIa") %>%
  rename(dna_id = value)
# Classify phylogroup XX taxa
phylogroup_xx_node <- MRCA(tree_3_10,  c("KX923038", "MH770720"))
phylogroup_xx_taxa <- tree_subset(tree_3_10, node = phylogroup_xx_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XX") %>%
  rename(dna_id = value)
# Classify phylogroup XVI.XVIII taxa
phylogroup_xvi_xviii_node <- MRCA(tree_3_10,  c("P1229_bin_25.fa", "MH770712"))
phylogroup_xvi_xviii_taxa <- tree_subset(tree_3_10, node = phylogroup_xvi_xviii_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XVI.XVIII") %>%
  rename(dna_id = value)
# Classify phylogroup XL taxa
phylogroup_xl_node <- MRCA(tree_3_10,  c("MH770644", "MH770659"))
phylogroup_xl_taxa <- tree_subset(tree_3_10, node = phylogroup_xl_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XL") %>%
  rename(dna_id = value)
# Classify phylogroup XIX taxa
phylogroup_xix_node <- MRCA(tree_3_10,  c("KX922887", "P1574_bin_1.fa"))
phylogroup_xix_taxa <- tree_subset(tree_3_10, node = phylogroup_xix_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "XIX") %>%
  rename(dna_id = value)
# Classify species complex 3-10a
complex_3_10_node <- MRCA(tree_3_10, c("P1574_bin_1.fa", "MH770515"))
complex_3_10_taxa <- tree_subset(tree_3_10, node = complex_3_10_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(species_complex = "3.10a") %>%
  rename(dna_id = value)
# Classify section 3-10
section_3_10_taxa <- tree_3_10$tip.label %>%
  as_tibble() %>%
  mutate(section = "3.10") %>%
  rename(dna_id = value)
# Join all assignments
clade_assignments_3_10 <- bind_rows(phylogroup_xxii_taxa, phylogroup_xxiii_taxa, phylogroup_xiii_xliii_taxa, 
                                   phylogroup_viiia_taxa, phylogroup_xx_taxa, phylogroup_xvi_xviii_taxa, 
                                   phylogroup_xl_taxa, phylogroup_xix_taxa) %>%
  right_join(section_3_10_taxa, by = "dna_id") %>%
  left_join(complex_3_10_taxa, by = "dna_id")

# Section 3-11

# Load and root section tree the same way as we plotted before
outgroup_3_11 <- c("NMS1_bin_5.fa", "NMS2_bin_6.fa")
tree_3_11 <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_11.treefile") %>%
  root(outgroup = outgroup_3_11, resolve.root = T)
# Classify phylogroup I? taxa
phylogroup_i_node <- MRCA(tree_3_11,  c("NMS1_bin_5.fa", "NMS2_bin_6.fa"))
phylogroup_i_taxa <- tree_subset(tree_3_11, node = phylogroup_i_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "I") %>%
  rename(dna_id = value)
# Classify phylogroup VIIIb taxa
phylogroup_viiib_node <- MRCA(tree_3_11,  c("KX922897", "KX923088"))
phylogroup_viiib_taxa <- tree_subset(tree_3_11, node = phylogroup_viiib_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "VIIIb") %>%
  rename(dna_id = value)
# Classify phylogroup NEW5 taxa
phylogroup_new5_node <- MRCA(tree_3_11,  c("P8256_bin_3.fa", "MH770580"))
phylogroup_new5_taxa <- tree_subset(tree_3_11, node = phylogroup_new5_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "VIIIb") %>%
  rename(dna_id = value)
# Classify species complex 3-11a
complex_3_11_node <- MRCA(tree_3_11, c("KX922897", "P6963_bin_9.fa"))
complex_3_11_taxa <- tree_subset(tree_3_11, node = complex_3_11_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(species_complex = "3.11a") %>%
  rename(dna_id = value)
# Classify section 3-11
section_3_11_taxa <- tree_3_11$tip.label %>%
  as_tibble() %>%
  mutate(section = "3.11") %>%
  rename(dna_id = value)
# Join all assignments
clade_assignments_3_11 <- bind_rows(phylogroup_i_taxa, phylogroup_viiib_taxa, phylogroup_new5_taxa) %>%
  right_join(section_3_11_taxa, by = "dna_id") %>%
  left_join(complex_3_11_taxa, by = "dna_id")

# Section 3-12

# Load and root section tree the same way as we plotted before
outgroup_3_12 <- c("P9820_bin_6.fa", "P9820")
tree_3_12 <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_12.treefile") %>%
  root(outgroup = outgroup_3_12, resolve.root = T)
# Classify phylogroup NEW6 taxa
phylogroup_new6_taxa <- c("P9820_bin_6.fa", "P9820") %>%
  as_tibble() %>%
  mutate(phylogroup = "NEW6") %>%
  rename(dna_id = value)
# Classify phylogroup XIV taxa
phylogroup_xiv_taxa <- c("KX923070", "KX923071", "KX923069") %>%
  as_tibble() %>%
  mutate(phylogroup = "XIV") %>%
  rename(dna_id = value)
# Classify species complex 3-11a
complex_3_12_node <- MRCA(tree_3_12, c("AJ632063", "KX923069"))
complex_3_12_taxa <- tree_subset(tree_3_12, node = complex_3_12_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(species_complex = "3.12a") %>%
  rename(dna_id = value)
# Classify section 3-12
section_3_12_taxa <- tree_3_12$tip.label %>%
  as_tibble() %>%
  mutate(section = "3.12") %>%
  rename(dna_id = value)
# Join all assignments
clade_assignments_3_12 <- bind_rows(phylogroup_new5_taxa, phylogroup_xiv_taxa) %>%
  right_join(section_3_12_taxa, by = "dna_id") %>%
  left_join(complex_3_12_taxa, by = "dna_id")

# Section 2-1

# Load and root section tree the same way as we plotted before
outgroup_2_1 <- c("EF102347")
tree_2_1 <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_2_1.treefile") %>%
  root(outgroup = outgroup_2_1, resolve.root = T)
# Classify section 2-1 taxa
section_2_1_node <- MRCA(tree_2_1,  c("DQ266020", "P12559_bin_8.fa"))
section_2_1_taxa <- tree_subset(tree_2_1, node = section_2_1_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(section = "2.1") %>%
  rename(dna_id = value) %>%
  add_column(phylogroup = NA) %>%
  add_column(species_complex = "2.1a")
# Classify aff section 2-1 taxa
aff_section_2_1_node <- MRCA(tree_2_1,  c("DQ185294", "DQ266002"))
aff_section_2_1_taxa <- tree_subset(tree_2_1, node = aff_section_2_1_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  add_row(value = "EF102347") %>%
  mutate(section = "aff_2.1") %>%
  rename(dna_id = value) %>%
  add_column(phylogroup = NA) %>%
  add_column(species_complex = NA)
# Join all assignments
clade_assignments_2_1 <- bind_rows(section_2_1_taxa, aff_section_2_1_taxa)

# Section 2-2

# Load and root section tree the same way as we plotted before
outgroup_2_2 <- c("S8_bin_2.fa", "S9_bin_6.fa", "KX923058", "EU877488", "P12591_bin_6.fa")
tree_2_2 <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_2_2.treefile") %>%
  root(outgroup = outgroup_2_2, resolve.root = T) 
# Classify phylogroup IX
phylogroup_ix_taxa <- c("KX923082", "KX923081") %>%
  as_tibble() %>%
  rename(dna_id = value) %>%
  add_column(phylogroup = "IX")
# Classify species complex 2.2a
complex_2_2_taxa <- tree_2_2$tip.label %>%
  as_tibble() %>%
  mutate(section = "2.2") %>%
  rename(dna_id = value) %>%
  add_column(species_complex = "2.2a")
# Join all assignments
clade_assignments_2_2 <- left_join(complex_2_2_taxa, phylogroup_ix_taxa, 
                                   by = "dna_id") 

# Section 2-3

# Only one genome in the section, and sequences that were placed with can be 94%
# similar so, will consider them aff_2.3
clade_assignments_2_3 <- read.FASTA(file = "analyses/species_delimitation/rbclx/clade_assignment/alignments/rbclx_2_3_aln.fna") %>%
  names() %>%
  as_tibble() %>%
  rename(dna_id = value) %>%
  mutate(section = ifelse(dna_id == "P12588_bin_4.fa",
                          "2.3", "aff_2.3")) %>%
  add_column(species_complex = NA) %>%
  add_column(phylogroup = NA)

# Section 2-4

# Load and root section tree the same way as we plotted before
tmp_tree <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_2_4.treefile")
root_node_2_4 <- MRCA(tmp_tree, c("P8690_bin_1.fa", "JL33_bin_16.fa", "P10247_bin_16.fa")) 
outgroup_2_4 <- tree_subset(tmp_tree, node = root_node_2_4, levels_back = 0)$tip.label
tree_2_4 <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_2_4.treefile") %>%
  root(outgroup = outgroup_2_4, resolve.root = T)
# Classify phylogroup III taxa
phylogroup_iii_node <- MRCA(tree_2_4,  c("P10247_bin_16.fa", "P8690_bin_1.fa"))
phylogroup_iii_taxa <- tree_subset(tree_2_4, node = phylogroup_iii_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "III") %>%
  rename(dna_id = value)
# Classify phylogroup IVa taxa
phylogroup_iva_node <- MRCA(tree_2_4,  c("MH757069", "NOS_bin_1.fa"))
phylogroup_iva_taxa <- tree_subset(tree_2_4, node = phylogroup_iva_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "IVa") %>%
  rename(dna_id = value)
# Classify phylogroup IVb taxa
phylogroup_ivb_node <- MRCA(tree_2_4,  c("P10894_bin_1.fa", "MH757080"))
phylogroup_ivb_taxa <- tree_subset(tree_2_4, node = phylogroup_ivb_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "IVb") %>%
  rename(dna_id = value)
# Classify phylogroup IVc taxa
phylogroup_ivc_node <- MRCA(tree_2_4,  c("P9728_bin_6.fa", "NMS7_bin_22.fa"))
phylogroup_ivc_taxa <- tree_subset(tree_2_4, node = phylogroup_ivc_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(phylogroup = "IVc") %>%
  rename(dna_id = value)
# Classify section and complex 2-4
section_2_4_taxa <- tree_2_4$tip.label %>%
  as_tibble() %>%
  add_column(section = "2.4") %>%
  add_column(species_complex = "2.4") %>%
  rename(dna_id = value)
# Join all assignments
clade_assignments_2_4 <- bind_rows(phylogroup_iii_taxa, phylogroup_iva_taxa, phylogroup_ivb_taxa,
                                   phylogroup_ivc_taxa) %>%
  right_join(section_2_4_taxa, by = "dna_id")

# Subclade 1

# No phylogroups delimited within subclade 1, so will consider only the subclade
# level
clade_assignments_subclade_1 <- read.FASTA(file = "analyses/species_delimitation/rbclx/clade_assignment/alignments/rbclx_subclade1_aln_edited.fna") %>%
  names() %>%
  as_tibble() %>%
  rename(dna_id = value) %>%
  add_column(subclade = "subclade_1") %>%
  add_column(section = NA) %>%
  add_column(species_complex = NA) %>%
  add_column(phylogroup = NA)

# Tables with all assignments

# Bring in unassigned queries
unassigned_queries <- scan(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/placement/unassigned_queries.txt", 
                           what = "character") %>%
  str_remove("QUERY___")
# Join all asignments
clade_assignments_all <- bind_rows(clade_assignments_2_1, clade_assignments_2_2,
                                   clade_assignments_2_3, clade_assignments_2_4, clade_assignments_3_1,
                                   clade_assignments_3_2, clade_assignments_3_3, clade_assignments_3_4,
                                   clade_assignments_3_5, clade_assignments_3_6, clade_assignments_3_7,
                                   clade_assignments_3_8, clade_assignments_3_9, clade_assignments_3_10,
                                   clade_assignments_3_11, clade_assignments_3_12, clade_assignments_subclade_1) %>%
  mutate(subclade = ifelse(is.na(subclade),
                           case_when(str_detect(section, "2\\.") ~ "subclade_2",
                                     str_detect(section, "3\\.") ~ "subclade_3"),
                           subclade)) %>%
  add_row(dna_id = unassigned_queries, subclade = "unassigned")
# Public clade assignments
clade_assignments_public <- public_rbclx_metadata %>%
  distinct(rbclx_accession, .keep_all = T) %>%
  left_join(clade_assignments_all, by = c("rbclx_accession" = "dna_id"))
# Write tables
write.csv(clade_assignments_all, file = "analyses/species_delimitation/rbclx/clade_assignment/clade_assignments_all.csv",
          row.names = F)

z
#### PLOT TREES WITH SPATIAL DATA ####

# Function to prepare ABMI site data for plotting
prep_abmi_data <- function(clade_assignments, source_data, abmi_metadata) {
  # Get clade assignments for ABMI taxa in focal tree
  clade_assignments_abmi <- clade_assignments %>%
    left_join(source_data, by = c("dna_id" = "taxon")) %>%
    filter(source == "abmi")
  # Get list of ABMI taxa in focal tree
  abmi_taxa <- clade_assignments_abmi %>%
    pull(dna_id)
  # Get counts of OTU per site
  # An TOU could be a phylogroup, or if no phylogroup was assigned, then all
  # floating taxa within a species complex were considered the same OTU
  cooccurrence_df <- clade_assignments_abmi %>%
    left_join(abmi_metadata, by = c("dna_id" = "duke_dna")) %>%
    mutate(best_otu = ifelse(!is.na(phylogroup),
                             phylogroup, species_complex)) %>%
    filter(!is.na(site)) %>% # FIGURE OUT WHY P6314, P6535 and P6314 are not in DB
    select(site, best_otu) %>%
    distinct() %>% # KEEPS ONLY ONE DUPLICATE BUT BOTH RECORDS HAVE THE SAME SITE SO IT'S SAFE
    group_by(site) %>%
    summarise(n_otu = length(best_otu))
  # Join cooccurrence df with cyano db and natural regions from site data
  site_data <- left_join(clade_assignments_abmi, abmi_metadata, 
                             by = c("dna_id" = "duke_dna")) %>%
    distinct(dna_id, .keep_all = T) %>% # KEEPS ONLY ONE DUPLICATE BUT BOTH RECORDS HAVE THE SAME SITE SO IT'S SAFE
    filter(!is.na(site)) %>% # FIGURE OUT WHY P6314, P6535 and P6314 are not in DB
    left_join(cooccurrence_df, by = "site")
  # Results into table
  result <- list(site_data = site_data, abmi_taxa = abmi_taxa, 
                 clade_assignments_abmi = clade_assignments_abmi)
}
# Load ABMI data
cyano_db <- read.delim("misc_files/ABMI_cyanolichen_db_T-BAS_v6.csv",
                       header = T, sep = ",") %>%
  select(duke_dna, site)
# Prep site and list of taxa for abmi plots
abmi_data_3_1 <- prep_abmi_data(clade_assignments_3_1, query_source_df, cyano_db)
abmi_data_3_5 <- prep_abmi_data(clade_assignments_3_5, query_source_df, cyano_db)
abmi_data_3_6 <- prep_abmi_data(clade_assignments_3_6, query_source_df, cyano_db)
abmi_data_2_4 <- prep_abmi_data(clade_assignments_2_4, query_source_df, cyano_db)
# Get list of ABMI taxa in focal trees
abmi_taxa_3_1 <- abmi_data_3_1$abmi_taxa
abmi_taxa_3_5 <- abmi_data_3_5$abmi_taxa
abmi_taxa_3_6 <- abmi_data_3_6$abmi_taxa
abmi_taxa_2_4 <- abmi_data_2_4$abmi_taxa
# Get site data for focal trees
site_data_3_1 <- abmi_data_3_1$site_data
site_data_3_5 <- abmi_data_3_5$site_data
site_data_3_6 <- abmi_data_3_6$site_data
site_data_2_4 <- abmi_data_2_4$site_data
# Get clade assignments for focal trees
clade_assignments_2_4_abmi <- abmi_data_2_4$clade_assignments_abmi
# Remove reference and public taxa from focal trees
tree_3_1_abmi <- keep.tip(tree_3_1, abmi_taxa_3_1)
tree_3_5_abmi <- keep.tip(tree_3_5, abmi_taxa_3_5)
tree_3_6_abmi <- keep.tip(unroot(tree_3_6), abmi_taxa_3_6)
tree_2_4_abmi <- keep.tip(unroot(tree_2_4), abmi_taxa_2_4)
# Join site data to trees for plotting
tree_3_1_abmi_meta <- as_tibble(tree_3_1_abmi) %>%
  left_join(site_data_3_1, by = c("label" = "dna_id")) %>%
  as.treedata()
tree_3_5_abmi_meta <- as_tibble(tree_3_5_abmi) %>%
  left_join(site_data_3_5, by = c("label" = "dna_id")) %>%
  as.treedata()
tree_3_6_abmi_meta <- as_tibble(tree_3_6_abmi) %>%
  left_join(site_data_3_6, by = c("label" = "dna_id")) %>%
  as.treedata()
tree_2_4_abmi_meta <- as_tibble(tree_2_4_abmi) %>%
  left_join(site_data_2_4, by = c("label" = "dna_id")) %>%
  as.treedata()
# Get node numbers for phylogroups
abmi_v_node <- MRCA(tree_3_6_abmi, 
                    filter(clade_assignments_3_6_abmi, phylogroup == "V") %>%
                      pull(dna_id) %>% 
                      unique()
                    )
abmi_xlii_node <- MRCA(tree_3_6_abmi, 
                    filter(clade_assignments_3_6_abmi, phylogroup == "XLII") %>%
                      pull(dna_id) %>% 
                      unique()
                    )
abmi_iva_node <- MRCA(tree_2_4_abmi, 
                       filter(clade_assignments_2_4_abmi, phylogroup == "IVa") %>%
                         pull(dna_id) %>% 
                         unique()
)
abmi_ivb_node <- MRCA(tree_2_4_abmi, 
                      filter(clade_assignments_2_4_abmi, phylogroup == "IVb") %>%
                        pull(dna_id) %>% 
                        unique()
)
abmi_ivc_node <- MRCA(tree_2_4_abmi, 
                      filter(clade_assignments_2_4_abmi, phylogroup == "IVc") %>%
                        pull(dna_id) %>% 
                        unique()
)
abmi_iii_node <- MRCA(tree_2_4_abmi, 
                      filter(clade_assignments_2_4_abmi, phylogroup == "III") %>%
                        pull(dna_id) %>% 
                        unique()
)
# Colors for cooccurrence
colors_3 <- c("1" = "#F29AAA", "2" = "#1F78B4", "3" = "#B2DF8A")
colors_4 <- c("1" = "#F29AAA", "2" = "#1F78B4", "3" = "#B2DF8A", "4" = "#AD6E1A")
# Plot trees with cooccurrence data
tree_3_6_cooccurrence_plot <- ggtree(tree_3_6_abmi_meta) +
  geom_tippoint(aes(color = factor(n_otu), alpha = 0.1)) +
  geom_cladelab(node = abmi_v_node, label = "V") +
  geom_cladelab(node = abmi_xlii_node, label = "XLII") +
  scale_color_manual(name = "No. of co-occurring OTUs\nin collection site", values = colors_3) +
  xlim(0, 0.037) +
  ylim(0, 460) +
  geom_treescale(width = 0.005, y = 50, x = 0.003, offset = 3.5) +
  theme(legend.position = c(0.05, 0.7),
        legend.justification = "left", 
        legend.direction = "vertical")
tree_2_4_cooccurrence_plot <- ggtree(tree_2_4_abmi_meta) +
  geom_tippoint(aes(color = factor(n_otu), alpha = 0.1), na.rm = T) +
  geom_cladelab(node = abmi_iva_node, label = "IVa") +
  geom_cladelab(node = abmi_ivb_node, label = "IVb") +
  geom_cladelab(node = abmi_ivc_node, label = "IVc") +
  geom_cladelab(node = abmi_iii_node, label = "III") +
  scale_color_manual(name = "No. of co-occurring OTUs\nin collection site", values = colors_4) +
  geom_treescale(width = 0.005, y = 50, x = 0.003, offset = 3.5)

ggtree(tree_3_1_abmi_meta) +
  geom_tippoint(aes(color = factor(n_otu), alpha = 0.1), na.rm = T)
ggtree(tree_3_5_abmi_meta) +
  geom_tippoint(aes(color = factor(n_otu), alpha = 0.1), na.rm = T)
# Save tree plots
ggsave(tree_3_6_cooccurrence_plot, 
       filename = "document/plots/tree_3_6_cooccurrence_plot.pdf")
ggsave(tree_2_4_cooccurrence_plot, 
       filename = "document/plots/tree_2_4_cooccurrence_plot.pdf") 
  
# Site data for plotting maps in ArcGIS

# Add count of number of specimens and lat lons
# 3-1
site_data_3_1_plot_sum <- site_data_3_1 %>%
  group_by(site) %>%
  summarise(n_specimens = length(dna_id), n_otu = unique(n_otu)) %>%
  left_join(all_site_data, by = "site") %>%
  select(site, n_specimens, n_otu, lat, long) %>%
  filter(!is.na(lat))
# 3-5
site_data_3_5_plot_sum <- site_data_3_5 %>%
  group_by(site) %>%
  summarise(n_specimens = length(dna_id), n_otu = unique(n_otu)) %>%
  left_join(all_site_data, by = "site") %>%
  select(site, n_specimens, n_otu, lat, long) %>%
  filter(!is.na(lat))
# 3-6
site_data_3_6_plot_sum <- site_data_3_6 %>%
  group_by(site) %>%
  summarise(n_specimens = length(dna_id), n_otu = unique(n_otu)) %>%
  left_join(all_site_data, by = "site") %>%
  select(site, n_specimens, n_otu, lat, long) %>%
  filter(!is.na(lat))
# 2-4
site_data_2_4_plot_sum <- site_data_2_4 %>%
  group_by(site) %>%
  summarise(n_specimens = length(dna_id), n_otu = unique(n_otu)) %>%
  left_join(all_site_data, by = "site") %>%
  select(site, n_specimens, n_otu, lat, long) %>%
  filter(!is.na(lat))
# Save site data as csvs
dir.create("analyses/species_delimitation/rbclx/cooccurrence_maps")
write_csv(site_data_3_1_plot_sum, file = "analyses/species_delimitation/rbclx/cooccurrence_maps/site_data_3_1_plot_sum.csv")
write_csv(site_data_3_5_plot_sum, file = "analyses/species_delimitation/rbclx/cooccurrence_maps/site_data_3_5_plot_sum.csv")
write_csv(site_data_3_6_plot_sum, file = "analyses/species_delimitation/rbclx/cooccurrence_maps/site_data_3_6_plot_sum.csv")
write_csv(site_data_2_4_plot_sum, file = "analyses/species_delimitation/rbclx/cooccurrence_maps/site_data_2_4_plot_sum.csv")

#### PLOT TREES WITH ALIGNMENT ####

# Load alignment with unique haplotypes. I generated the fasta from the phylip
# produced by iqtree. I removed the 3' end of the rbcL that is only present in the
# taxa from set103 and removed those taxa from the alignment
haplotype_seqs_v <- read.FASTA("analyses/species_delimitation/cooccurrence/trees/focal/rbclx_set103_global_abmi_v.uniqueseq.fna")
# Subset focal tree to haplotyes
haplotype_labels_v <- names(haplotype_seqs_v)
# Subset focal tree to haplotype taxa
tree_v_haplotypes <- keep.tip(tree_v, tip = haplotype_labels_v)
# Plot focal tree with alignment
tree_v_haplotypes_plot <- ggtree(tree_v_haplotypes) +
  geom_treescale(width = 0.005, y = 70, x = 0.003, offset = 2) +
  theme(legend.position = "none")
tree_v_haplotypes_msa <- msaplot(tree_v_haplotypes_plot, haplotype_seqs_v)
# Save plot
ggsave(tree_v_haplotypes_msa, filename = "document/plots/tree_haplotypes_v_msa.pdf")

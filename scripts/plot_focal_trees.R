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
# Get df key to replace tip labes
tiplabel_key <- query_source_df %>%
  filter(source == "public" | source == "reference_tree") %>%
  left_join(public_rbclx_metadata, by = c("taxon" = "rbclx_accession")) %>%
  left_join(ani95_clusters, by = c("taxon" = "genome")) %>%
  left_join(popcogent_clusters, by = c("taxon" = "Strain")) %>%
  mutate(new_label = if_else(source == "public",
                             paste(taxon, dna_source, region, "Phylogroup", nostoc_phylogroup),
                             paste(taxon, ani_95_cluster, Main_cluster))) %>%
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
  outgroup_taxa = outgroup_3_1, font_size = 3, col_offset = 0.1)
ggsave(clade_3_1_public_plot, file = "document/plots/clade_3_1_public.pdf",
       units = "cm", height = 40, width = 40)
# Clade 3-4
outgroup_3_4 <- c("KX923097")
clade_3_4_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_4.treefile",
  outgroup_taxa = outgroup_3_4, font_size = 3, col_offset = 0.1)
ggsave(clade_3_4_public_plot, file = "document/plots/clade_3_4_public.pdf")
# Clade 3-5
outgroup_3_5 <- c("DQ185281", "KX922920", "KX922926", "KX922928", "S67_bin_3.fa",
                  "KX922927","KX922930", "KX922971", "KX923043", 
                  "P3034_bin_5.fa", "S66_bin_1.fa")
clade_3_5_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_5_with_spacer.treefile",
  outgroup_taxa = outgroup_3_5, font_size = 3, col_offset = 0.08)
ggsave(clade_3_5_public_plot, file = "document/plots/clade_3_5_public.pdf")
# Clade 3-6
tmp_tree <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_6.treefile")
#root_node_3_6 <- MRCA(tmp_tree, c("P8202_bin_9.fa", "P8569_bin_1.fa")) 
root_node_3_6 <- MRCA(tmp_tree, c("KX922980", "P8569_bin_1.fa")) 
outgroup_3_6 <- tree_subset(tmp_tree, node = root_node_3_6, levels_back = 0)$tip.label
clade_3_6_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_6.treefile",
  outgroup_taxa = outgroup_3_6, font_size = 3, col_offset = 0.08)
ggsave(clade_3_6_public_plot, file = "document/plots/clade_3_6_public.pdf",
       units = "cm", height = 40, width = 40)
# Clade 3-7
tmp_tree <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_7.treefile")
root_node_3_7 <- MRCA(tmp_tree, c("Z94893", "P3068_bin_7.fa")) 
outgroup_3_7 <- tree_subset(tmp_tree, node = root_node_3_7, levels_back = 0)$tip.label
clade_3_7_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_7.treefile",
  outgroup_taxa = outgroup_3_7, font_size = 3, col_offset = 0.08)
ggsave(clade_3_7_public_plot, file = "document/plots/clade_3_7_public.pdf")
# Clade 3-8
outgroup_3_8 <- c("Nostoc_KVJ2.fa", "Nostoc_KVS11.fa")
clade_3_8_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_8.treefile",
  outgroup_taxa = outgroup_3_8, font_size = 3, col_offset = 0.08)
ggsave(clade_3_8_public_plot, file = "document/plots/clade_3_8_public.pdf")
# Clade 3-9
outgroup_3_9 <- c("KX922962", "KX922961")
clade_3_9_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_9.treefile",
  outgroup_taxa = outgroup_3_9, font_size = 3, col_offset = 0.08)
ggsave(clade_3_9_public_plot, file = "document/plots/clade_3_9_public.pdf")
# Clade 3-10
tmp_tree <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_10.treefile")
root_node_3_10 <- MRCA(tmp_tree, c("P2037_bin_28.fa", "P2039_bin_23.fa")) 
outgroup_3_10 <- tree_subset(tmp_tree, node = root_node_3_10, levels_back = 0)$tip.label
clade_3_10_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_10.treefile",
  outgroup_taxa = outgroup_3_10, font_size = 3, col_offset = 0.08)
ggsave(clade_3_10_public_plot, file = "document/plots/clade_3_10_public.pdf",
       units = "cm", height = 40, width = 40)
# Clade 3-11
outgroup_3_11 <- c("NMS1_bin_5.fa", "NMS2_bin_6.fa")
clade_3_11_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_11.treefile",
  outgroup_taxa = outgroup_3_11, font_size = 3, col_offset = 0.09)
ggsave(clade_3_11_public_plot, file = "document/plots/clade_3_11_public.pdf")
# Clade 3-12
outgroup_3_12 <- c("P9820_bin_6.fa")
clade_3_12_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_3_12.treefile",
  outgroup_taxa = outgroup_3_12, font_size = 3, col_offset = 0.09)
ggsave(clade_3_12_public_plot, file = "document/plots/clade_3_12_public.pdf")
# Clade 2-1
outgroup_2_1 <- c("EF102347")
clade_2_1_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_2_1.treefile",
  outgroup_taxa = outgroup_2_1, font_size = 3, col_offset = 0.09)
ggsave(clade_2_1_public_plot, file = "document/plots/clade_2_1_public.pdf")
# Clade 2-2
outgroup_2_2 <- c("S8_bin_2.fa", "S9_bin_6.fa", "KX923058", "EU877488", "P12591_bin_6.fa")
clade_2_2_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_2_2.treefile",
  outgroup_taxa = outgroup_2_2, font_size = 3, col_offset = 0.09)
ggsave(clade_2_2_public_plot, file = "document/plots/clade_2_2_public.pdf")
# Clade 2-3
# Clade 2-4
tmp_tree <- read.tree(file = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_2_4.treefile")
root_node_2_4 <- MRCA(tmp_tree, c("P8690_bin_1.fa", "JL33_bin_16.fa", "P10247_bin_16.fa")) 
outgroup_2_4 <- tree_subset(tmp_tree, node = root_node_2_4, levels_back = 0)$tip.label
clade_2_4_public_plot <- prep_public_focal_trees(
  treefile = "analyses/species_delimitation/rbclx/clade_assignment/trees/focal/rbclx_2_4.treefile",
  outgroup_taxa = outgroup_2_4, font_size = 3, col_offset = 0.09)
ggsave(clade_2_4_public_plot, file = "document/plots/clade_2_4_public.pdf",
       units = "cm", height = 30, width = 30)



  clade_3_5_public_plot <- ggtree(clade_3_5_public_bundle$tree) +
  geom_tiplab(align = T, size = 3)
clade_3_5_public_plot <- gheatmap(clade_3_5_public_plot, clade_3_5_public_bundle$df, 
                                  width = 0.1, offset = 0.08, colnames = F) +
  theme(legend.position = "none")


#### PLOT TREES WITH SPATIAL DATA ####

# Load backbone tree
backbone_tree <- read.tree("analyses/species_delimitation/cooccurrence/trees/placement/backbone.tree")
# Get labels of reference taxa to remove from the focal ml trees 
reference_node_v <- MRCA(backbone_tree, c("JL34_bin_22.fa", "P2162_bin_10.fa"))
reference_tree_v <- tree_subset(backbone_tree, node = reference_node_v,
                                levels_back = 0)
reference_taxa_v <- reference_tree_v$tip.label
# Load focal ml trees
tree_v <- read.tree("analyses/species_delimitation/cooccurrence/trees/focal/rbclx_set103_global_abmi_v.treefile")
# Remove reference taxa from focal trees
tree_v <- drop.tip(tree_v, reference_taxa_v)

# Classify taxa from focal groups according to clades

# Load ABMI data
cyano_db <- read.delim("misc_files/ABMI_cyanolichen_db_T-BAS_v6.csv",
                       header = T, sep = ",") %>%
  select(duke_dna, site)
# Load site metadata
all_site_data <- read.delim("misc_files/cyano_db_site_data.csv", header = T, sep = ",")
# Split focal group V into clades V, XLII, and XXXIX
clade_v_anchor <- c("P2176", "P2156")
clade_v_node <- MRCA(tree_v, clade_v_anchor)
clade_v_taxa <- tree_subset(tree_v, node = clade_v_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(clade = "V") %>%
  rename(duke_dna = value)
clade_xlii_anchor <- c("P8203", "P10943")
clade_xlii_node <- MRCA(tree_v, clade_xlii_anchor)
clade_xlii_taxa <- tree_subset(tree_v, node = clade_xlii_node, levels_back = 0)$tip.label %>%
  as_tibble() %>%
  mutate(clade = "XLII") %>%
  rename(duke_dna = value)
clade_xxxix_taxa <- tree_v$tip.label %>%
  setdiff(c(clade_v_taxa$duke_dna, clade_xlii_taxa$duke_dna)) %>%
  as_tibble() %>%
  mutate(clade = "XXXIX") %>%
  rename(duke_dna = value)
clade_assignments_v <- bind_rows(clade_v_taxa, clade_xlii_taxa, clade_xxxix_taxa)
# Get counts of clades per site
cooccurrence_df <- left_join(cyano_db, clade_assignments_v, by = "duke_dna") %>%
  filter(!is.na(clade)) %>%
  select(2:3) %>%
  distinct() %>%
  group_by(site) %>%
  summarise(n_clades = length(clade))
# Join cooccurrence df with cyano db and natural regions from site data
# Remove row 67 which has a wrong dupplicate
# FIGURE OUT DUPLICATES
site_data_v <- left_join(clade_assignments_v, cyano_db, by = "duke_dna") %>%
  filter(!row_number() %in% 67) %>%
  left_join(cooccurrence_df, by = "site") %>%
  left_join(select(all_site_data, site, nat_reg, nat_subreg), by = "site") %>%
  mutate(nat_reg = 
           na_if(nat_reg, "#VALUE!")) %>%
  mutate(nat_subreg = 
           na_if(nat_subreg, "#VALUE!")) %>%
  distinct()
# Join site data to trees for plotting
tree_v_meta <- as_tibble(tree_v) %>%
  left_join(site_data_v, by = c("label" = "duke_dna")) %>%
  as.treedata()
# Plot tree with cooccurrence data
tree_v_cooccurrence_plot <- ggtree(tree_v_meta) +
  geom_tippoint(aes(color = factor(n_clades), alpha = 0.3), na.rm = T) +
  geom_cladelab(node = clade_v_node, label = "V") +
  geom_cladelab(node = clade_xlii_node, label = "XLII") +
  geom_cladelab(node = 456, label = "XXXIX") + 
  scale_color_discrete(name = "No. of co-occurring clades") +
  xlim(0, 0.037) +
  ylim(0, 460) +
  geom_treescale(width = 0.005, y = 50, x = 0.003, offset = 3.5) +
  theme(legend.position = c(0.05, 0.7),
        legend.justification = "left", 
        legend.direction = "vertical")
# Save tree plot
ggsave(tree_v_cooccurrence_plot, 
       filename = "document/plots/cooccurrence_tree_v.pdf")
  


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

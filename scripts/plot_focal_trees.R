#!/usr/bin/env Rscript

# Load required packages
library(tidyverse)
library(data.table)
library(ape)
library(treeio)
library(ggtree)

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

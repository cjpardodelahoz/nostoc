#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(ggimage)
library(ggpubr)
library(ggnewscale)
library(tidytree)
library(treeio)
library(labdsv)
library(MSCquartets)
source("scripts/r_functions.R")

##### CONFLICT PIE CHARTS AND LIFESTYLE #####

# Load and root the dated tree (wASTRAL topology)
dated_tree <- read.tree("analyses/phylogenetics/set103/divtime/1_part/mcmc/c1/dated.tree")
dated_tree <- dated_tree %>%
  treeio::root(outgroup = c("Aphanizomenon_flos_aquae_NIES_81.fa",
                            "Anabaena_cylindrica_PCC_7122.fa",
                            "Cylindrospermum_stagnale_PCC_7417.fa"), 
               resolve.root = T)
# Get the number of taxa in the tree
ntaxa <- dated_tree$tip.label %>%
  length()
# Load the wAstral vs gene trees Discovista result
discov_out <- read_csv("analyses/phylogenetics/set103/conflict/single_vs_wastral/discovista_out/single.metatable.results.csv")
# Sumarize discov output per bipartition and modify the node column (bipart+ntaxa)
discov_df_top <- condense_discov_out_weak(discov_out) %>%
  select(bipart, percent_concordant:percent_weak_support) %>%
  pivot_longer(percent_concordant:percent_weak_support, 
               names_to = "support", values_to = "percent") %>%
  mutate(node = 
           as.numeric(bipart) + ntaxa)
# Split the discovista output into a list of tables by node number to generate
# the list of pie charts
discov_df_list <- discov_df_top %>%
  group_split(node)
names(discov_df_list) <- discov_df_top %>% 
  pull(node) %>% 
  unique() %>% 
  sort
# Load lifestyle metadata and join it to the tree
lifestyle_metadata <- read_csv("misc_files/set103_lifestyle_metadata.csv") %>%
  mutate(tip_label =
           paste(tip_label, ".fa", sep = ""))
# Get node numbers for major clades
outgroup_node <- MRCA(dated_tree, c("Cylindrospermum_stagnale_PCC_7417.fa", "Anabaena_cylindrica_PCC_7122.fa"))
subclade1_node <- MRCA(dated_tree, c("Nostoc_sp_JC1668.fa", "Nostoc_linckia_z4.fa"))
subclade2_node <- MRCA(dated_tree, c("P12588_bin_4.fa", "NOS_bin_1.fa"))
subclade3a_node <- MRCA(dated_tree, c("P9820_bin_6.fa", "P539_bin_14.fa"))
subclade3b_node <- MRCA(dated_tree, c("P943_bin_5.fa", "Nmoss2.fa"))
# Named vector of colors for piecharts
pie_colors <- c("percent_concordant" = "#40549f", "percent_discordant" = "#ea4753", 
                "percent_uninformative" = "#9b979c", "percent_weak_reject" = "#e8db7e", 
                "percent_weak_support" = "#4599ad")
# Colors for lifestyle categories
lifestyle_colors <- c("bryophyte_associated" = "#C8F32F", "cycad_symbiont" = "#a484f4", 
                      "Free-living" = "#E36A86", "lichenized" = "#51a9a6")
# Generate list of piecharts
conflict_pies <- discov_df_list %>%
  map(~ggplot(.x, aes(x = "", y = percent, fill = support)) +
        geom_bar(stat = "identity") +
        coord_polar("y", start = 0) +
        scale_fill_manual(values = pie_colors) +
        theme_void() +
        theme(legend.position = "none"))
# Plot the dated tree with the piecharts
tree_plot <- ggtree(dated_tree, right = T) %<+% lifestyle_metadata
tree_pies <- tree_plot +
  geom_inset(conflict_pies, width = 0.02, height = 0.02, x = "node") +
  new_scale_color() +
  geom_fruit(geom = geom_point, 
             mapping = aes(y = label, color = lifestyle), shape = "circle", 
             size = 1.3, offset = 0.05) +
  scale_color_manual(values = lifestyle_colors, na.value = "white")+
  geom_cladelab(node = subclade1_node, label = "subclade 1", 
                angle = 270, offset.text= 0.03, 
                offset= 0.075,
                fontsize = 3.5, barsize = 0.8) +
  geom_cladelab(node = subclade2_node, label = "subclade 2", 
                angle = 270, offset.text= 0.03, 
                offset= 0.075,
                fontsize = 3.5, barsize = 0.8) +
  geom_cladelab(node = subclade3a_node, label = "subclade 3a", 
                angle = 270, offset.text = 0.03,
                offset= 0.075,
                fontsize = 3.5, barsize = 0.8) +
  geom_cladelab(node = subclade3b_node, label = "subclade 3b", 
                angle = 270, offset.text = 0.03, 
                offset= 0.075, 
                fontsize = 3.5, barsize = 0.8) +
  geom_cladelab(node = outgroup_node, label = "outgroup", 
                angle = 0, offset.text = 0.03, 
                offset= 0.075, 
                fontsize = 3.5, barsize = 0.8)
ggsave(plot = tree_pies, "document/plots/dated_tree_pies_lifestyle.pdf",
       units = "cm", width = 27, height = 27,device = "pdf")

##### CONFLICT VS TIME #####

# Join dfs with branch lengths and conflict info from discovista
tree_df <- as_tibble(dated_tree) %>%
  left_join(discov_df_top, by = "node") %>%
  filter(!is.na(bipart)) %>%
  mutate(branch.length =
           branch.length*1000)
# Plot support classes vs internode lengths
discordant_vs_time <- tree_df %>%
  filter(support == "percent_discordant") %>%
  ggplot(aes(x = branch.length, y = percent)) +
  geom_point(color = "#ea4753") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Internode length (Myr)", y = "Percent of discordant trees") +
  theme(panel.background = NULL, 
        panel.border = element_rect(fill = "transparent", linewidth = 0.75),
        axis.text = element_text(size = 12, color = "black"))
concordant_vs_time <- tree_df %>%
  filter(support == "percent_concordant") %>%
  ggplot(aes(x = branch.length, y = percent)) +
  geom_point(color = "#40549f") +
  labs(x = "Internode length (Myr)", y = "Percent of concordant trees") +
  theme(panel.background = NULL, 
        panel.border = element_rect(fill = "transparent", linewidth = 0.75),
        axis.text = element_text(size = 12, color = "black"))
weak_support_vs_time <- tree_df %>%
  filter(support == "percent_weak_support") %>%
  ggplot(aes(x = branch.length, y = percent)) +
  geom_point(color = "#4599ad") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Internode length (Myr)", y = "Percent of weakly concordant trees") +
  theme(panel.background = NULL, 
        panel.border = element_rect(fill = "transparent", linewidth = 0.75),
        axis.text = element_text(size = 12, color = "black"))
weak_reject_vs_time <- tree_df %>%
  filter(support == "percent_weak_reject") %>%
  ggplot(aes(x = branch.length, y = percent)) +
  geom_point(color = "#e8db7e") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Internode length (Myr)", y = "Percent of weakly discordant trees") +
  theme(panel.background = NULL, 
        panel.border = element_rect(fill = "transparent", linewidth = 0.75),
        axis.text = element_text(size = 12, color = "black"))
# Arrange plots
main_plots <- ggarrange(ncol = 1, nrow = 2, 
                          concordant_vs_time, weak_reject_vs_time)
suppl_plots <- ggarrange(ncol = 2, nrow = 1, 
                        weak_support_vs_time, discordant_vs_time)
# Save pdfs
ggsave(main_plots, filename = "document/plots/conflict_vs_time_main.pdf", 
       height = 6, width = 3)
ggsave(suppl_plots, filename = "document/plots/conflict_vs_time_suppl.pdf", 
       height = 3, width = 6)

##### LABELED WASTRAL TREE PLOT #####

# Load table with taxon name key
key_df <- read_csv("document/tables/voucher_v1.csv") %>%
  select(genome_id, taxon_name, country) %>%
  mutate(taxon_name = 
           ifelse(taxon_name == genome_id,
                  taxon_name,
                  paste(genome_id, taxon_name, sep = "_"))) %>%
  mutate(genome_id = 
           paste(genome_id, ".fa", sep = "")) %>%
  mutate(taxon_name = 
           paste(taxon_name, country, sep = "_")) %>%
  select(genome_id, taxon_name)
# A toy df for plotting full label
toy_df <- key_df %>% distinct() %>%
  column_to_rownames(var = "taxon_name") %>%
  mutate(genome_id = "not_a_variable")
# Load wASTRAL tree
wastral_tree <- read.tree(file = "analyses/phylogenetics/set103/trees/astral/wastral.tree")
wastral_tree$edge.length[is.na(wastral_tree$edge.length)] <- 1
ggtree(wastral_tree)
# Rename tree tips
tree_key <- wastral_tree$tip.label %>%
  as_tibble() %>%
  left_join(key_df, by = c("value" = "genome_id"))
wastral_tree <- rename_taxa(wastral_tree, data = tree_key)
# Plot the tree with cluster and tip labels
wastra_tree_plot <- ggtree(wastral_tree) +
  geom_tiplab(size = 2) +
  geom_treescale()
tree_plot_tmp <- gheatmap(wastra_tree_plot, toy_df, width = 0.1, offset = 35, colnames = F) +
  theme(legend.position = "none")
# Save the tree to check
ggsave(tree_plot_tmp, filename = "document/plots/wastral_labeled_left.pdf",
       width = 8, height = 12)


##### SIMPLEX PLOTS #####

# Load gene trees
trees <- read.tree(file = "analyses/phylogenetics/set103/trees/astral/ml_gene.trees")
# Load taxon names
taxon_names <- scan(file = "scripts/genome_ids_set103", what = "character")
# Get table of quartet counts from the gene trees
quartet_counts <- quartetTable(trees = trees, taxonnames = taxon_names)

#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse)
library(ggtree)
library(ggimage)
library(ggpubr)
library(tidytree)
library(treeio)
library(labdsv)
library(MSCquartets)
source("scripts/r_functions.R")

##### CONFLICT PIE CHARTS #####

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
# Named vector of colors for piecharts
pie_colors <- c("percent_concordant" = "#40549f", "percent_discordant" = "#ea4753", 
                "percent_uninformative" = "#9b979c", "percent_weak_reject" = "#e8db7e", 
                "percent_weak_support" = "#4599ad")
# Generate list of piecharts
conflict_pies <- discov_df_list %>%
  map(~ggplot(.x, aes(x = "", y = percent, fill = support)) +
        geom_bar(stat = "identity") +
        coord_polar("y", start = 0) +
        scale_fill_manual(values = pie_colors) +
        theme_void() +
        theme(legend.position = "none"))
# Plot the dated tree with the piecharts
tree_plot <- ggtree(dated_tree, right = T)
tree_pies <- tree_plot +
  geom_inset(conflict_pies, width = 0.02, height = 0.02, x = "node")
ggsave(plot = tree_pies, "document/plots/dated_tree_pies.pdf",
       units = "cm", width = 15, height = 27,device = "pdf")

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

##### SIMPLEX PLOTS #####

# Load gene trees
trees <- read.tree(file = "analyses/phylogenetics/set103/trees/astral/ml_gene.trees")
# Load taxon names
taxon_names <- scan(file = "scripts/genome_ids_set103", what = "character")
# Get table of quartet counts from the gene trees
quartet_counts <- quartetTable(trees = trees, taxonnames = taxon_names)

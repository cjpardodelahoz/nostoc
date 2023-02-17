#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse)
library(ape)
source("scripts/r_functions.R")

# Read and tidy DiscoVista output for the L1648+ng+shet gene vs concat conflict analysis

# Load discov output table
discov_out <- read_csv("analyses/phylogenetics/set73/conflict/discovista_out/single.metatable.results.csv")
# Sumarize discov output per bipartition
discov_df_top <- condense_discov_out_weak(discov_out) %>%
  select(bipart, percent_concordant:percent_weak_support) %>%
  pivot_longer(percent_concordant:percent_weak_support, names_to = "support", values_to = "percent") 

# Named vector of colors for piecharts
pie_colors <- c("percent_concordant" = "#40549f", "percent_discordant" = "#ea4753", 
                "percent_uninformative" = "#9b979c", "percent_weak_reject" = "#e8db7e", 
                "percent_weak_support" = "#4599ad")
# Generate pies
conflict_pies <- ggplot(discov_df_top, aes(x = "", y = percent, fill = support)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = pie_colors) +
  theme_void() +
  facet_wrap("bipart")
ggsave(plot = conflict_pies, "analyses/phylogenetics/set73/conflict/conflict_pies.pdf",
       device = "pdf")

# remove long branch from concat tree
concat_tree <- read.tree(file = "analyses/phylogenetics/set73/trees/concat/concat_ng_na.treefile")
concat_tree <- drop.tip(concat_tree, c("P10247_bin.15.fa", "Aphanizomenon_flos_aquae_NIES_81.fa"))
write.tree(concat_tree, file = "analyses/phylogenetics/set73/conflict/concat_tree.tree")

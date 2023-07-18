#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse) 
library(ggtree)
library(tidytree)
library(ggplot2)
library(labdsv)
library(data.table)
library(ggalt)
library(ggpubr)
source("scripts/r_functions.R")

##### SUMMARIZE FASTANI RESULTS #####

# Outgroup taxa
outgroup = c("Aphanizomenon_flos_aquae_NIES_81.fa",
             "Anabaena_cylindrica_PCC_7122.fa",
             "Cylindrospermum_stagnale_PCC_7417.fa")
# Load FastANI output
fastani_df <- read_delim(file = "analyses/species_delimitation/fastani/set12c/fastani_set12c_ql_out", 
                         col_names = FALSE) %>%
  mutate(X1 = 
           str_remove(X1, "analyses/cyano_genomes/set12c/")) %>%
  mutate(X1 = 
           str_remove(X1, "_chromosome")) %>%
  mutate(X2 = 
           str_remove(X2, "analyses/cyano_genomes/set12c/")) %>%
  mutate(X2 = 
           str_remove(X2, "_chromosome")) %>%
  mutate(alignment_fraction = 
           X4/X5) %>%
  rowwise() %>%
  mutate(genome1 = 
           min(c(X1, X2))) %>%
  mutate(genome2 = 
           max(c(X1, X2))) %>%
  ungroup() %>%
  distinct(genome1, genome2, .keep_all = T) %>%
  filter(genome1 %nin% outgroup & genome2 %nin% outgroup) %>%
  select(!c(X1, X2, X4, X5)) %>%
  rename(ani = X3) %>%
  relocate(genome1, genome2)
# Calculate ANI gap
gap_df <- gap(fastani_df, low_lim = 88, hi_lim = 98)

# Plots

# Histogram of pairwise ANIS
ani_hist <- ggplot(data = fastani_df, aes(x = ani)) +
  geom_histogram(fill = "gray40") +
  labs(x = "ANI", y = "Frequency") +
  scale_x_continuous(n.breaks = 13) +
  geom_vline(xintercept = 95, linetype = "dashed", color = "gray40") +
  theme(panel.background = NULL, 
        panel.border = element_rect(fill = "transparent", linewidth = 0.75),
        axis.text = element_text(size = 12, color = "black"))
# ANI vs alignment fraction
ani_vs_aln <- ggplot(data = fastani_df, aes(x = ani, y = alignment_fraction)) +
  geom_point(size = 0.5, alpha = 0.5) +
  labs(x = "ANI", y = "Alignment fraction") +
  scale_x_continuous(n.breaks = 13) +
  geom_vline(xintercept = 95, linetype = "dashed", color = "gray40") +
  theme(panel.background = NULL, 
        panel.border = element_rect(fill = "transparent", linewidth = 0.75),
        axis.text = element_text(size = 12, color = "black"))
# Distribution of ANI gaps
ani_gap_all <- ggplot(gap_df, aes(y = reorder(genomes, gap_low_lim), 
                                  x = gap_low_lim, 
                                  xend = gap_hi_lim)) +
  geom_dumbbell(size_x = 0.2, size_xend = 0.2, size = 0.2) +
  geom_vline(xintercept = 95, linetype = "dashed", color = "gray40") +
  labs(x = "ANI gap span", y = "Genomes") +
  scale_x_continuous(n.breaks = 10, limits = c(80, 98)) +
  theme(panel.background = NULL,
        panel.border = element_rect(fill = "transparent", linewidth = 0.75),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
# Arrange Fig
ani_sum_fig <- ggarrange(ncol = 3, nrow = 1, 
                         ani_hist, ani_vs_aln, ani_gap_all)
# Save plots
ggsave(ani_sum_fig, filename = "document/plots/ani_sum_fig.pdf", 
       width = 11.5, height = 4)


##### CLUSTER GENOMES USING ANI #####

# Load and root the dated tree (wASTRAL topology)
dated_tree <- read.tree("analyses/phylogenetics/set103/divtime/1_part/mcmc/c1/dated.tree")
dated_tree <- dated_tree %>%
  treeio::root(outgroup = c("Aphanizomenon_flos_aquae_NIES_81.fa",
                            "Anabaena_cylindrica_PCC_7122.fa",
                            "Cylindrospermum_stagnale_PCC_7417.fa"), 
               resolve.root = F)
# Converta ANI data to matrix
ani_matrix <- read_delim(file = "analyses/species_delimitation/fastani/set12c/fastani_set12c_ql_out", 
                         col_names = FALSE) %>%
  mutate(X1 = 
           str_remove(X1, "analyses/cyano_genomes/set12c/")) %>%
  mutate(X1 = 
           str_remove(X1, "_chromosome")) %>%
  mutate(X2 = 
           str_remove(X2, "analyses/cyano_genomes/set12c/")) %>%
  mutate(X2 = 
           str_remove(X2, "_chromosome")) %>%
  select(1:3) %>% 
  as.matrix() %>%
  matrify() %>%
  mutate_if(is.character, as.numeric)
ani_matrix[ani_matrix == 0] <- 60
# Make sure that the matrix is symmetric by equating upper and lower triangle
#ani_matrix[lower.tri(ani_matrix)] <- ani_matrix[upper.tri(ani_matrix)]
# cluster the genomes using 95% ANI
ani_95_clusters <- abs(100-ani_matrix) %>%
  as.dist() %>%
  hclust() %>%
  cutree(h = 5) %>%
  enframe(name = "genome", value = "ani_95_cluster") %>%
  arrange(ani_95_cluster) %>%
  mutate(ani_95_cluster = 
           paste("c", ani_95_cluster, sep = ""))
# Save ANI cluster results
write_csv(ani_95_clusters, 
          file = "analyses/species_delimitation/fastani/set12c/ani_95_clusters.csv")



plot_ani_clusters <- function(tree, ani_cluster_df) {
  # Getting tree df
  tree <- dated_tree
  ani_cluster_df <- ani_95_clusters
  tree_df <- dated_tree %>%
    tidytree::as.treedata() %>%
    tidytree::as_tibble()
  # Tree plot 
  tree_plot <- ggtree(tree)
  # Clusters
  clusters <- ani_cluster_df %>%
    dplyr::pull(ani_95_cluster) %>%
    unique()
  n_clusters <-  length(clusters)
  #
  for (i in 1:n_clusters) {
    cluster_genomes <- ani_cluster_df %>%
      dplyr::filter(ani_95_cluster == clusters[i]) %>%
      dplyr::pull(genome)
    node_number <- tidytree::MRCA(tree_df, cluster_genomes) %>%
      dplyr::pull(node)
    tree_plot <- tree_plot +
      geom_cladelab(node = node_number, label = clusters[i],
                    barsize = 4, lineheight = 4)
  }
  
}

ggtree(dated_tree) + 
  geom_cladelab(node=153, label="", lineheight = 2, barsize = 4)


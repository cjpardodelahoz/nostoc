#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse)
library(fastbaps)
library(ape)
library(treeio)
library(ggtree)
library(phytools)
source("scripts/r_functions.R")

#### FASTBAPS CLUSTERING SET103 ####


# Load and root the dated tree (wASTRAL topology)
dated_tree <- read.tree("analyses/phylogenetics/set103/divtime/1_part/mcmc/c1/dated.tree")
# Sequence data
concat_ng <- "analyses/phylogenetics/set103/alignments/concat/concat_ng.fna"
# Import as sparse data matrix and set prior to optimised baps
sparse_data <- import_fasta_sparse_nt(concat_ng)
sparse_data <- optimise_prior(sparse_data, type = "optimise.baps")
# Cluster seqs with FastBAPS
baps_hc <- fast_baps(sparse_data, k.init = 120)
# Get best BAPS clustering scheme
best_partition <- best_baps_partition(sparse_data, baps_hc)
best_partition_tree <- best_baps_partition(sparse_data, dated_tree)
#
plot.df <- data.frame(id = colnames(sparse_data$snp.matrix), fastbaps = best_partition, 
                      stringsAsFactors = FALSE)
plot.df2 <- data.frame(id = dated_tree$tip.label, fastbaps = best_partition_tree, 
                      stringsAsFactors = FALSE)


#### PLOT CLUSTER RESULTS ####

# Trees and cluster data

# Define outgroup taxa to remove them
outgroup <- c("Aphanizomenon_flos_aquae_NIES_81.fa",
              "Anabaena_cylindrica_PCC_7122.fa",
              "Cylindrospermum_stagnale_PCC_7417.fa")
# Load the dated tree (wASTRAL topology) and remove the outgroup
dated_tree <- read.tree("analyses/phylogenetics/set103/divtime/1_part/mcmc/c1/dated.tree")
# Load ANI95 clusters
ani_95_clusters <- read_csv(file = "analyses/species_delimitation/fastani/set12c/ani_95_clusters.csv") %>%
  mutate(ani_95_cluster = 
           str_remove(ani_95_cluster, "c")) %>%
  mutate(ani_95_cluster = 
           as.integer(ani_95_cluster))
# Load PopCOGenT clusters
popcogent_clusters <- read_delim(file = "analyses/species_delimitation/popcogent/set12c/set12c_0.000355362.txt.cluster.tab.txt",
                                 delim = "\t") %>%
  select(Strain, Main_cluster) %>%
  mutate(Main_cluster =
           Main_cluster + 1) %>%
  mutate(Strain = 
           str_replace(Strain, "_chromosome", ".fa"))

# Transform cluster assignments into binary variable for plotting

# Function to transform cluster assignments into binary variable for plotting with
# gheatmap(). This gets around the fact that ggtree doesn't plot clade bars for
# a single taxon, and keeps it at two colors only
# tree_plot     Plot of the tree with the taxa in the order you want them displayed. Generate with ggtree(tree)
# cluster_df    Data frame with the cluster assignations for each taxon in the tree
# label_col     Character string. Column label in the df with the tip labels
# cluster_col   Character string. Column label in the df with cluster assignments
# bin_col       Character string. Label for the new column for the bins
bin_tip_clusters <- function(tree_plot, cluster_df, label_col, cluster_col,
                             bin_col) {
  # Get ordered tip labels
  ordered_tips <- get_taxa_name(tree_plot) %>%
    as_tibble() %>%
    rename(tip_label = value)
  # Get cluster ids order by tip labels as shown on tree plot
  ordered_clusters <- left_join(ordered_tips, cluster_df, 
                                       by = c("tip_label" = label_col)) %>%
    select(all_of(cluster_col)) %>%
    distinct()
  # Vector with sequence of cluster bins
  cluster_bins <- rep(c("black", "gray"), length.out = nrow(ordered_clusters))
  # New df with cluster bins ready for gheatmap plot
  cluster_bin_df <- bind_cols(ordered_clusters, cluster_bins) %>%
    rename(!!bin_col := `...2`) %>%
    right_join(cluster_df, by = cluster_col) %>%
    select(all_of(c(label_col, bin_col)))
}
# Transform cluster assignments for plotting
ani_95_bins <- bin_tip_clusters(tree_plot = ggtree(dated_tree), 
                                cluster_df = ani_95_clusters, 
                                label_col = "genome", 
                                cluster_col = "ani_95_cluster", 
                                bin_col = "ani_95_bin")
popcogent_bins <- bin_tip_clusters(tree_plot = ggtree(dated_tree), 
                                   cluster_df = popcogent_clusters, 
                                   label_col = "Strain", 
                                   cluster_col = "Main_cluster", 
                                   bin_col = "popcogent_bin")
# Join all cluster bins
all_cluster_bins <- left_join(ani_95_bins, popcogent_bins, 
                              by = c("genome" = "Strain")) %>%
  column_to_rownames(var = "genome")

# Plot the tree

# Plot the dated tree with the ANI and popcogent clusters mapped
t <- ggtree(dated_tree)
clusters_plot <- gheatmap(t, all_cluster_bins, width = 0.1, offset = -0.035, colnames = F) +
  scale_fill_manual(values = c("gray" = "#AAA7BF", "black" = "#1176A5")) +
  theme(legend.position = "none")
# Save the plot
ggsave(clusters_plot, filename = "document/plots/tree_clusters.pdf", 
       height = 5, width = 2.5)


# Get ordered tip labels
ordered_tips <- get_taxa_name(ggtree(dated_tree)) %>%
  as_tibble() %>%
  rename(tip_label = value)
ani <- ordered_tips %>%
  left_join(ani_95_clusters, by = c("tip_label" = "genome")) %>%
  distinct(ani_95_cluster) %>%
  add_column(x = seq(1:76)) %>%
  left_join(ani_95_clusters, by = "ani_95_cluster") %>%
  select(genome, x)
popcogent <- ordered_tips %>%
  left_join(popcogent_clusters, by = c("tip_label" = "Strain")) %>%
  distinct(Main_cluster) %>%
  add_column(x = seq(1:107)) %>%
  left_join(popcogent_clusters, by = "Main_cluster") %>%
  select(Strain, x)

f4 <- facet_plot(t, panel = "ANI 95%", data = ani, geom = geom_tile, 
                 aes(x = x), color = "brown")
f5 <- facet_plot(f4, panel = "popcogent", data = popcogent, geom = geom_tile, 
                 aes(x = x), color = "brown")

f4 <- facet_plot(t, panel = "ANI 95%", data = ani_95_clusters, geom = geom_tile, 
                 aes(x = ani_95_cluster), color = "brown")
f5 <- facet_plot(f4, panel = "popcogent", data = popcogent_clusters, geom = geom_tile, 
                 aes(x = Main_cluster), color = "brown")


  
  
  
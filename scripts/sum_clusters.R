#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse)
library(fastbaps)
library(ape)
library(phytools)
library(treeio)
library(ggtree)
library(aplot)
library(tidytree)
library(ggplot2)
library(labdsv)
library(data.table)
library(ggalt)
library(ggpubr)
source("scripts/r_functions.R")

#### PARSE AND SUMMARIZE FASTANI RESULTS #####

## Cluster genomes with ANI 95%

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

## Parse fastani results for histogram and alignment fraction calculation

# Load FastANI output as a data frame and remove duplicate pairwise entries
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

## ANI gap calculation

# Outgroup taxa
outgroup = c("Aphanizomenon_flos_aquae_NIES_81.fa",
             "Anabaena_cylindrica_PCC_7122.fa",
             "Cylindrospermum_stagnale_PCC_7417.fa")
# Get DF for plotting ANI gap with all ANI points
fastani_df2 <- ani_matrix %>%
  rownames_to_column(var = "genome_id") %>%
  pivot_longer(cols = 2:152, names_to = "genome_pair", values_to = "ani") %>%
  filter(genome_id %nin% outgroup) %>%
  filter(ani > 79 & ani < 100)


#### PLOT CLUSTER RESULTS ####

## Tree plot with clusters and species complex clades

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
# Load and root the dated tree (wASTRAL topology)
dated_tree <- read.tree("analyses/phylogenetics/set103/divtime/1_part/mcmc/c1/dated.tree") %>%
  drop.tip(c(outgroup))
# Load PopCOGenT clusters
popcogent_clusters <- read_delim(file = "analyses/species_delimitation/popcogent/set12c/set12c_0.000355362.txt.cluster.tab.txt",
                                 delim = "\t") %>%
  select(Strain, Main_cluster) %>%
  mutate(Main_cluster =
           Main_cluster + 1) %>%
  mutate(Strain = 
           str_replace(Strain, "_chromosome", ".fa"))
# ANI95 clusters as integers
ani_95_clusters <- ani_95_clusters %>%
  mutate(ani_95_cluster = 
           str_remove(ani_95_cluster, "c")) %>%
  mutate(ani_95_cluster = 
           as.integer(ani_95_cluster))
# Genome ids binned by ANI95 for plotting
ani_95_bins <- bin_tip_clusters(tree_plot = ggtree(dated_tree), 
                                cluster_df = ani_95_clusters, 
                                label_col = "genome", 
                                cluster_col = "ani_95_cluster", 
                                bin_col = "ani_95_bin")
# Genome ids binned by popcogent clusters for plotting
popcogent_bins <- bin_tip_clusters(tree_plot = ggtree(dated_tree), 
                                   cluster_df = popcogent_clusters, 
                                   label_col = "Strain", 
                                   cluster_col = "Main_cluster", 
                                   bin_col = "popcogent_bin")
# Join both binned dfs
all_cluster_bins <- left_join(ani_95_bins, popcogent_bins,
                              by = c("genome" = "Strain")) %>%
  column_to_rownames(var = "genome")
# Get species complex nodes
node_a <- MRCA(dated_tree, c("P2081_bin_11.fa", "P8768_bin_2.fa"))
node_b <- MRCA(dated_tree, c("P6447_bin_3.fa",  "Nostoc_sp_Lobaria_pulmonaria_5183.fa"))
node_d <- MRCA(dated_tree, c("P6636_bin_11.fa", "P2213_bin_24.fa"))
node_e <- MRCA(dated_tree, c("S67_bin_3.fa", "P6602_bin_13.fa"))
node_f <- MRCA(dated_tree, c("P2162_bin_10.fa", "P8202_bin_9.fa"))
node_g <- MRCA(dated_tree, c("P6620_bin_4.fa", "P3068_bin_7.fa"))
node_h <- MRCA(dated_tree, c("Nostoc_KVS11.fa", "Nmoss2.fa"))
node_j <- MRCA(dated_tree, c("P2039_bin_23.fa", "S32_bin_15.fa"))
node_k <- MRCA(dated_tree, c("NMS1_bin_5.fa", "P6963_bin_9.fa"))
node_l <- MRCA(dated_tree, c("P9820_bin_6.fa", "Nostoc_commune_NIES_4072.fa"))
node_m <- MRCA(dated_tree, c("P12545_bin_5.fa", "P12564_bin_10.fa"))
node_n <- MRCA(dated_tree, c("P12537_bin_5.fa", "P12591_bin_6.fa"))
node_p <- MRCA(dated_tree, c("P10247_bin_16.fa", "NOS_bin_1.fa"))
# Plot the dated tree with the ANI and popcogent clusters mapped, and species
# complexes highlighted
t <- ggtree(dated_tree) +
  geom_highlight(node = node_a, fill = "black", alpha = 0.3) +
  geom_highlight(node = node_b, fill = "#939495", alpha = 0.25) +
  geom_highlight(node = node_d, fill = "#939495", alpha = 0.25) +
  geom_highlight(node = node_e, fill = "black", alpha = 0.3) +
  geom_highlight(node = node_f, fill = "#939495", alpha = 0.25) +
  geom_highlight(node = node_g, fill = "black", alpha = 0.3) +
  geom_highlight(node = node_h, fill = "#939495", alpha = 0.25) +
  geom_highlight(node = node_j, fill = "#939495", alpha = 0.25) +
  geom_highlight(node = node_k, fill = "black", alpha = 0.3) +
  geom_highlight(node = node_l, fill = "#939495", alpha = 0.25) +
  geom_highlight(node = node_m, fill = "black", alpha = 0.3) +
  geom_highlight(node = node_n, fill = "#939495", alpha = 0.25) +
  geom_highlight(node = node_p, fill = "#939495", alpha = 0.25)
tree_plot <- gheatmap(t, all_cluster_bins, width = 0.1, offset = -0.01, colnames = F) +
  scale_fill_manual(values = c("gray" = "#AAA7BF", "black" = "#1176A5"), na.value = "white") +
  theme(legend.position = "none")

## ANI plots

# Distribution of ANI gaps (all ANIs as dots)
ani_gap_dots <- ggplot(fastani_df2,  aes(y = genome_id, 
                                         x = ani)) +
  geom_point(size = 0.15)  +
  scale_x_continuous(n.breaks = 11, limits = c(80, 100)) +
  geom_vline(xintercept = 95, linetype = "dashed", color = "gray40") +
  labs(x = "ANI") +
  theme(panel.background = NULL,
        panel.border = element_rect(fill = "transparent", linewidth = 0.75),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())
ani_gap_dots_tree <- ani_gap_dots %>% insert_left(tree_plot)
# Histogram of pairwise ANIS
ani_hist <- ggplot(data = fastani_df, aes(x = ani)) +
  geom_histogram(fill = "gray40") +
  labs(x = "ANI", y = "Frequency") +
  scale_x_continuous(n.breaks = 13) +
  geom_vline(xintercept = 95, linetype = "dashed", color = "gray40") +
  theme(panel.background = NULL, 
        panel.border = element_rect(fill = "transparent", linewidth = 0.75),
        axis.text = element_text(size = 12, color = "black")) +
  annotate("rect", xmin = 83, xmax = 96, ymin = 0, ymax = 3100, fill = "grey70", 
           alpha = 0.3)
# ANI vs alignment fraction
ani_vs_aln <- ggplot(data = fastani_df, aes(x = ani, y = alignment_fraction)) +
  geom_point(size = 0.5, alpha = 0.5) +
  labs(x = "ANI", y = "Alignment fraction") +
  scale_x_continuous(n.breaks = 13) +
  geom_vline(xintercept = 95, linetype = "dashed", color = "gray40") +
  theme(panel.background = NULL, 
        panel.border = element_rect(fill = "transparent", linewidth = 0.75),
        axis.text = element_text(size = 12, color = "black")) +
  annotate("rect", xmin = 83, xmax = 96, ymin = 0, ymax = 1, fill = "grey70", 
           alpha = 0.3)
# Arrange ani hist and aln fraction
ani_traditional <- ggarrange(nrow = 1, ncol = 2, ani_hist, ani_vs_aln)
# Save plots
ggsave(ani_gap_dots_tree, filename = "document/plots/ani_gap_dots_tree.pdf", units = "cm", 
       width = 17.4, height = 13)
 ggsave(ani_traditional, filename = "document/plots/ani_traditional.pdf", 
       width = 7.7, height = 4)

## Old stuff for faceting

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






#### PLOT FOR CHECKING ####

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
# Rename cluster bins df
all_cluster_bins_tmp <- all_cluster_bins %>%
  rownames_to_column() %>%
  left_join(key_df, by = c("rowname" = "genome_id")) %>%
  select(taxon_name, ani_95_bin, popcogent_bin, bin_16s) %>%
  column_to_rownames(var = "taxon_name")
# Rename tree tips
tree_key <- dated_tree$tip.label %>%
  as_tibble() %>%
  left_join(key_df, by = c("value" = "genome_id"))
dated_tree_tmp <- rename_taxa(dated_tree, data = tree_key)
# Plot the tree with cluster and tip labels
t_tmp <- ggtree(dated_tree_tmp) +
  geom_tiplab(size = 2)
tree_plot_tmp <- gheatmap(t_tmp, all_cluster_bins_tmp, width = 0.1, offset = 1.1, colnames = F) +
  scale_fill_manual(values = c("gray" = "#AAA7BF", "black" = "#1176A5"), na.value = "white") +
  theme(legend.position = "none")
# Save the tree to check
ggsave(tree_plot_tmp, filename = "document/plots/cluster_tree_to_check.pdf",
       width = 8, height = 12)

# Load 16S alignment
aln_16s <- read.dna("analyses/phylogenetics/set103/alignments/single/16s_aln.fas", 
                    format = "fasta")
aln <- read.dna("seq.fasta", format = "fasta")
find_synapomorphies(aln, c("a"))
find_synapomorphies(aln_16s, c("Nostoc_sp_Peltigera_membranacea_cyanobiont_210A.fa", "P2083_bin_10.fa", "3_bin_5.fa"))

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


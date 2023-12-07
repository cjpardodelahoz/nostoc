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
# shapes for substrate categories
  substrate_shapes <- c("terricolous" = 1, "epiphytic" = 16,
                      "aquatic" = 4)
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
  scale_color_manual(values = lifestyle_colors, na.value = "white") +
  new_scale(new_aes = "shape") +
  geom_fruit(geom = geom_point, 
             mapping = aes(y = label, shape = lichen_substrate),
             offset = 0.09) +
  scale_shape_manual(values = substrate_shapes) +
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

##### MODEL AND PLOT CONFLICT VS TIME #####

# Join dfs with branch lengths and conflict info from discovista
tree_df <- as_tibble(dated_tree) %>%
  left_join(discov_df_top, by = "node") %>%
  filter(!is.na(bipart)) %>%
  mutate(branch.length =
           branch.length*1000)

# Strongly discordant vs time

# Strongly discordant data subset excluding branches = 0 because it is an artifact of the rooting
discordant_data <- filter(tree_df, support == "percent_discordant" & branch.length > 0)
# Linear regression
discordant_model <- lm(percent ~ log(branch.length),
                       data = discordant_data)
# Predict values for regression line
discordant_x_fit <- seq(min(discordant_data$branch.length), 
                        max(discordant_data$branch.length), 
                        length.out = length(discordant_data$branch.length))
discordant_y_fit <- predict(discordant_model, 
                            newdata = data.frame(branch.length = discordant_x_fit))
# Plot strong discordance
discordant_vs_time <- discordant_data %>%
  ggplot(aes(x = branch.length, y = percent)) +
  geom_point(color = "#ea4753") +
  geom_line(aes(x = discordant_x_fit, y = discordant_y_fit), color = "gray20", linetype = 2) +
  annotate("text", label = "P < 0.001 *", x = 400, y = 80, size = 4) +
  annotate("text", label = "R^2 == 0.40", x = 400, y = 70, size = 4, parse = T) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Internode length (Myr)", y = "Percent of discordant trees") +
  theme(panel.background = NULL, 
        panel.border = element_rect(fill = "transparent", linewidth = 0.75),
        axis.text = element_text(size = 12, color = "black"))

# Strongly concordant vs time

# Strongly discordant data subset excluding branches = 0 because it is an artifact of the rooting
concordant_data <- filter(tree_df, support == "percent_concordant" & branch.length > 0)
# Linear regression
concordant_model <- lm(percent ~ log(branch.length),
                       data = concordant_data)
# Predict values for regression line
concordant_x_fit <- seq(min(concordant_data$branch.length), 
                        max(concordant_data$branch.length), 
                        length.out = length(concordant_data$branch.length))
concordant_y_fit <- predict(concordant_model, 
                            newdata = data.frame(branch.length = discordant_x_fit))
# Plot strong concordance vs time with predicted regression line
concordant_vs_time <- concordant_data %>%
  ggplot(aes(x = branch.length, y = percent)) +
  geom_point(color = "#40549f") +
  geom_line(aes(x = concordant_x_fit, y = concordant_y_fit), color = "gray50", linetype = 2) +
  scale_y_continuous(limits = c(0, 100)) +
  annotate("text", label = "P < 0.001 *", x = 500, y = 50, size = 4) +
  annotate("text", label = "R^2 == 0.58", x = 500, y = 40, size = 4, parse = T) +
  labs(x = "Internode length (Myr)", y = "Percent of concordant trees") +
  theme(panel.background = NULL, 
        panel.border = element_rect(fill = "transparent", linewidth = 0.75),
        axis.text = element_text(size = 12, color = "black"))

# Weakly concordant vs time

# Weakly concordant data subset excluding branches = 0 because it is an artifact of the rooting
weak_support_data <- filter(tree_df, support == "percent_weak_support" & branch.length > 0)
# Linear regression to log of branch lengths -> NOT SIGNIFICANT< SKIP PLOTTING
weak_support_model <- lm(percent ~ branch.length,
                        data = weak_support_data)
# Plot weak support vs time
weak_support_vs_time <- weak_support_data %>%
  ggplot(aes(x = branch.length, y = percent)) +
  geom_point(color = "#4599ad") +
  annotate("text", label = "P = 0.205", x = 400, y = 80, size = 4) +
  annotate("text", label = "R^2 == 0.004", x = 400, y = 70, size = 4, parse = T) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Internode length (Myr)", y = "Percent of weakly concordant trees") +
  theme(panel.background = NULL, 
        panel.border = element_rect(fill = "transparent", linewidth = 0.75),
        axis.text = element_text(size = 12, color = "black"))

# Weakly discordant vs time

# Weakly discordant data subset excluding branche = 0 because it is an artifact of the rooting
weak_reject_data <- filter(tree_df, support == "percent_weak_reject" & branch.length > 0)
# Linear regression to log of branch lengths
weak_reject_model <- lm(percent ~ log(branch.length),
                       data = weak_reject_data)
# Predict values for regression line
weak_reject_x_fit <- seq(min(weak_reject_data$branch.length), 
                        max(weak_reject_data$branch.length), 
                        length.out = length(weak_reject_data$branch.length))
weak_reject_y_fit <- predict(weak_reject_model, 
                            newdata = data.frame(branch.length = weak_reject_x_fit))
# Plot weak conflict vs time with predicted regression line
weak_reject_vs_time <- weak_reject_data %>%
  ggplot(aes(x = branch.length, y = percent)) +
  geom_point(color = "#e8db7e") +
  geom_line(aes(x = weak_reject_x_fit, y = weak_reject_y_fit), color = "gray50", linetype = 2) +
  scale_y_continuous(limits = c(0, 100)) +
  annotate("text", label = "P < 0.001 *", x = 400, y = 80, size = 4) +
  annotate("text", label = "R^2 == 0.52", x = 400, y = 70, size = 4, parse = T) +
  labs(x = "Internode length (Myr)", y = "Percent of weakly discordant trees") +
  theme(panel.background = NULL, 
        panel.border = element_rect(fill = "transparent", linewidth = 0.75),
        axis.text = element_text(size = 12, color = "black"))

# Prepare plots for Fig 1B,C and Fig S2A,B

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



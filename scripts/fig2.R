#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse) 
library(ggtree)
library(tidytree)
library(ggplot2)
library(labdsv)
source("scripts/r_functions.R")

##### ANI MATRIX ####

# Load and root the dated tree (wASTRAL topology)
dated_tree <- read.tree("analyses/phylogenetics/set103/divtime/1_part/mcmc/c1/dated.tree")
dated_tree <- dated_tree %>%
  treeio::root(outgroup = c("Aphanizomenon_flos_aquae_NIES_81.fa",
                            "Anabaena_cylindrica_PCC_7122.fa",
                            "Cylindrospermum_stagnale_PCC_7417.fa"), 
               resolve.root = T)
# Load FastANI output
fastani_df <- read_delim(file = "analyses/species_delimitation/fastani/set12c/fastani_set12c_ql_out", 
                         col_names = FALSE) %>%
  mutate(X1 = 
           str_remove(X1, "analyses/cyano_genomes/set12c/")) %>%
  mutate(X2 = 
           str_remove(X2, "analyses/cyano_genomes/set12c/"))
# Converta ANI data to matrix
ani_matrix <- fastani_df %>%
  select(1:3) %>% 
  as.matrix() %>%
  matrify() %>%
  mutate_if(is.character, as.numeric)
ani_matrix[ani_matrix < 80] <- 80
ani_matrix <- ani_matrix %>%
  select(all_of(ordered_tips)) %>%
  as.data.frame()
rownames(ani_matrix) <- ordered_tips
# Plot tree with matrix
astral_tree_plot <- ggtree(astral_tree)
gheatmap(astral_tree_plot, ani_matrix, colnames=F) +
  scale_fill_continuous(type = "viridis")
ggsave(h, filename = "test.pdf", device = "pdf")
#

tree <- read.tree(text = "(((A,B),(C,D)),E);")
tree2 <- ladderize(tree, right = FALSE)
tree$tip.label
#> [1] "A" "B" "C" "D" "E"
tree2$tip.label
#> [1] "A" "B" "C" "D" "E"
plot(tree2)
nodelabels()
tiplabels()

astral_tree1 <- ladderize(astral_tree, right = F)
is_tip <- astral_tree1$edge[,2] <= length(astral_tree1$tip.label)
ordered_tips <- astral_tree1$edge[is_tip, 2]
ordered_tips <- astral_tree1$tip.label[ordered_tips] %>% rev()


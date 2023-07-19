#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse)
library(ggtree)
library(tidytree)
library(treeio)

# Load and root the dated tree (wASTRAL topology)
dated_tree <- read.tree("analyses/phylogenetics/set103/divtime/1_part/mcmc/c1/dated.tree")
dated_tree <- dated_tree %>%
  treeio::root(outgroup = c("Aphanizomenon_flos_aquae_NIES_81.fa",
                            "Anabaena_cylindrica_PCC_7122.fa",
                            "Cylindrospermum_stagnale_PCC_7417.fa"), 
               resolve.root = T)
# Get test node
node <- MRCA(dated_tree, c("Nostoc_sp_Peltigera_membranacea_cyanobiont_N6.fa", "P1574_bin_1.fa"))
test_tree <- tree_subset(dated_tree, node = node, levels_back = 0)
test_df <- data.frame(taxon = test_tree$tip.label) %>%
  mutate(file = 
           paste("analyses/cyano_genomes/set12c/", taxon, sep = "")) %>%
  mutate(file =
           str_replace(file, ".fa", "_chromosome.fa"))
dir.create("analyses/species_delimitation/gubbins/test", recursive = T)
write_delim(test_df, file = "analyses/species_delimitation/gubbins/test/test_df.tsv", 
            delim = "\t", col_names = F)

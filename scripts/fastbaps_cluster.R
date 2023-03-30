#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse)
library(fastbaps)
library(ape)
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

#
gg <- ggtree(dated_tree)

f2 <- facet_plot(gg, panel = "fastbaps", data = plot.df, geom = geom_tile, aes(x = fastbaps), 
                 color = "blue")
f3 <- facet_plot(f2, panel = "fastbaps level 2", data = plot.df2, geom = geom_tile, 
                 aes(x = fastbaps), color = "green")


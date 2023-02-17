#!/usr/bin/env Rscript

# USAGE:
#    get_delta_contigs_taxonomy_layer.R contig_collection_file taxonomy_file out_path
#
# contig_collection_file    The file with a table of split names and status of the contig after graphbin
# taxonomy_file             The tsv output of mmseqs taxonomy for the contigs in the metagenome
# out_path                  Base path for the output file of this script which will be used as visualization layer in anvio

# Load needed packages and functions
library(tidyverse)

# Get arguments from command line
commands <- commandArgs(trailingOnly = T)
# Define variables
contig_collection_file <- commands[1] 
taxonomy_file <- commands[2]
out_path <- commands[3]
bin_ids_graphbin_delta_file <- paste(out_path, "bin_ids_graphbin_delta",
                                     sep = "/")
anvio_layer_file <- paste(out_path, "delta_contigs_taxonomy_layer.tsv", 
                          sep = "/")

# load delta contig collection with split names table and get clean contig names
contig_collection <- read_delim(
  file = contig_collection_file, 
  col_names = T,
  delim = "\t") %>%
  mutate(contig =
           str_remove(split, "_split_00001"))
# Load taxonomy table
contig_taxonomy <- read_delim(file = taxonomy_file,
                              col_names = F,
                              delim = "\t")
# Add taxonomy to contig collection and select the split, bin and taxonomy cols
delta_contigs_taxonomy_layer <- contig_collection %>%
  left_join(contig_taxonomy, 
            by = c("contig" = "X1")) %>% 
  mutate(bin = 
           str_replace(bin, "green.", "added")) %>%
  mutate(bin = 
           str_replace(bin, "red.", "removed")) %>%
  mutate(bin =
           str_replace(bin, "black.", "unchanged")) %>%
  mutate(taxonomy = 
           gsub("'","", X4)) %>%
  mutate(taxonomy = 
           str_replace(taxonomy, "\\.","")) %>%
  mutate(taxonomy = 
           str_replace(taxonomy, "-","_")) %>%
  mutate(taxonomy = 
           gsub(" ","_", taxonomy)) %>% # used gsub here because there were unidentified whitespace characters
  select(split, bin, taxonomy)
# Write anvio layer file
write_delim(delta_contigs_taxonomy_layer, 
            file = anvio_layer_file,
            delim = "\t")
# Write file with list of bin ids in the collectiion
bin_ids <- contig_collection %>%
  pull(bin) %>%
  unique()
write(bin_ids, file = bin_ids_graphbin_delta_file)
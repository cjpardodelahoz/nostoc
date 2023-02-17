#!/usr/bin/env Rscript

# This is a script to get the range of bp from a blast output table. The table
# should be generated with the blast command -outfmt '6 sseqid sseq' 
# Usage:
# seq_range_from_blast.R blast_out_file range_file header_file
#
# blast_out_file  path to blast output table
# range_file      path to desired output file with seq range

library(tidyverse)

# Define variables from command line arguments
args <- commandArgs(trailingOnly = T)
blast_out_file <- args[1]
range_file <- args[2]
header_file <- args[3]
# Load blast output
blast_out <- read_delim(blast_out_file, col_names = F)
# Count number of target seqs and stop if there are more than 1
n_seqs <- blast_out %>%
  pull(X1) %>%
  unique() %>%
  length()
if (n_seqs > 1) {
  print("Error: more than one target sequence")
  stop()
}
# Get all sequence coordinates
seq_positions <- c(blast_out$X2, blast_out$X3)
# Get sequence range and header
seq_range <- paste(min(seq_positions), max(seq_positions), 
                   sep = ":")
seq_header <- blast_out %>%
  pull(X1) %>%
  unique()
# Write files with range and seq header
write(seq_range, file = range_file)
write(seq_header, file = header_file)


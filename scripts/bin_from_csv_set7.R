#!/usr/bin/env Rscript

# Load custom functions
source("scripts/r_functions.R")
library(tidyverse)
library(ape)


# Load sample ids
sample_ids <- scan("scripts/sample_ids_set7", what = "character")
# Load graphbin output
paths <- sapply(sample_ids, function(x) {
  paste(... = "analyses/bins/graphbin/", ...= x, "/graphbin2_output.csv", sep = "")
   }
  )
# Load bin csv and add sample id column and col labels
bin_df <- readr::read_csv(paths, col_names = FALSE, id = "path") %>%
  mutate(sample =
           stringr::str_remove(path, "analyses/bins/graphbin/")) %>%
  mutate(sample =
           stringr::str_remove(sample, "/graphbin2_output.csv"))
colnames(bin_df) <- c("path", "name", "bin", "sample")
# Get paths to contigs
contigs_paths <- sapply(sample_ids, function(x) {
   paste(... = "analyses/assemblies/", ...= x, "/contigs.fasta", sep = "")
     }
    )  
# Bin contigs for all samples
v_bin_from_csv(sample_id = sample_ids, bin_df = bin_df, 
               contigs_file = contigs_paths, out_dir = "analyses/bins/graphbin")

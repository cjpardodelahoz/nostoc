#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse)
source("scripts/r_functions.R")

# Load the nucleotide alignment summary trimmed to no gaps
aln_summary <- read_delim(file = "analyses/phylogenetics/set103/alignments/sumaries/summary_single_na_ng.txt") %>%
  mutate(busco_id = 
           str_remove(Alignment_name, "_.*"))
# Filter to busco loci that are present in more than 135 taxa (90%) and with more than 200 variable sites
kept_busco_ids <- aln_summary %>%
  filter(No_variable_sites >= 200 & No_of_taxa > 135) %>%
  pull(busco_id)
# Write file with kept busco_ids
write(kept_busco_ids, file = "scripts/busco_ids_filtered_set103")

         
#!/usr/bin/env Rscript

# Load custom functions
source("scripts/r_functions.R")
library(tidyverse)
# Load list of taxa
sample_ids_nic <- scan(file = "scripts/sample_ids_nic", 
                                  what = "character")
# Combine checkm output for all bins from nic assemblies
checkm_nic <- sum_checkm(sample_ids_nic, 
                                    checkm_out_dir = "analyses/genome_qc/metabat/nic",
                                    out_file_name = "checkm_out",
                                    bin_dir = "analyses/bins/metabat",
                                    bin_ext = ".fa")
# Filter checkm df to cyano bins
checkm_cyano_nic <- checkm_nic %>%
  dplyr::filter(
    stringr::str_detect(string = marker_lineage, pattern = "Cyanobacteria*.")
  )
# Test which samples have multiple cyano bins
multi_cyano_nic <- checkm_cyano_nic %>%
  pull(sample_id) %>% duplicated()
# Add column indicating samples that have multiple cyanos
checkm_cyano_nic <- checkm_cyano_nic %>%
  tibble::add_column(multi = multi_cyano_nic)
# Write records with multiple cyano bins
checkm_cyano_nic %>%
  filter(multi == TRUE) %>%
  write.table(file = "scripts/multi_cyano_nic", quote = FALSE)
# Get file with cyano bin paths
checkm_cyano_nic %>%
  dplyr::pull(bin_path) %>%
  write(file = "scripts/bin_paths_nic")

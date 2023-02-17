# Load required libraries and custom functions
library(tidyverse)
library(phylotools)
source("scripts/r_functions.R")

# Load busco output for set 1
sample_ids <- scan("scripts/sample_ids_set1", what = "character")
busco_test <- sum_busco(sample_ids = sample_ids, 
                        busco_out_dir = "analyses/genome_qc/set1/busco/by_taxon", 
                        busco_db = "run_cyanobacteria_odb10", 
                        out_file_name = "full_table.tsv")

busco_test_1 <- add_busco_paths(busco_test, 
                                busco_out_dir = "analyses/genome_qc/set1/busco/by_taxon",
                                busco_db = "run_cyanobacteria_odb10") %>%
  add_busco_labels()

cat_busco_seqs(busco = "146at1117",
               busco_df_paths_labels = busco_test_1,
               data_type = "DNA",
               out_dir = "test/test")

sort_busco_seqs(sample_ids = sample_ids,
                busco_out_dir = "analyses/genome_qc/set1/busco/by_taxon", 
                busco_db = "run_cyanobacteria_odb10", 
                out_file_name = "full_table.tsv",
                data_type = "DNA",
                out_dir = "analyses/phylogenetics/set1/seqs")

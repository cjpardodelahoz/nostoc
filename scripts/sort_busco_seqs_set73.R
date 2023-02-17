#!/usr/bin/env Rscript

# Load required libraries and custom functions
library(tidyverse)
library(phylotools)
source("scripts/r_functions.R")

# Load sample ids from set 1
sample_ids <- scan("scripts/genome_ids_set73", what = "character")
# Get busco na seq files
sort_busco_seqs(sample_ids = sample_ids,
                busco_out_dir = "analyses/genome_qc/set73/busco/by_taxon", 
                busco_db = "run_cyanobacteria_odb10", 
                out_file_name = "full_table.tsv",
                data_type = "DNA",
                out_dir = "analyses/phylogenetics/set73/seqs", 
                busco_ids_file = "scripts/busco_ids_set73")
# Get busco aa seq files
sort_busco_seqs(sample_ids = sample_ids,
                busco_out_dir = "analyses/genome_qc/set73/busco/by_taxon", 
                busco_db = "run_cyanobacteria_odb10", 
                out_file_name = "full_table.tsv",
                data_type = "AA",
                out_dir = "analyses/phylogenetics/set73/seqs", 
                busco_ids_file = "scripts/busco_ids_set73")

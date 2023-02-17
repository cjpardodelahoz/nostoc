#!/usr/bin/env Rscript

source("scripts/r_functions.R")
library(tidyverse)

#
genome_ids <- scan("scripts/genome_ids_set8", what = "character") 
# Load and merge busco output set7 and set8
busco_out_set7 <- merge_busco(sample_ids = genome_ids, 
                              busco_out_dir = "analyses/genome_qc/set7/busco/by_taxon/", 
                              busco_db = "run_nostocales_odb10", 
                              out_file_name = "full_table.tsv") 
busco_out_set8 <- merge_busco(sample_ids = genome_ids, 
                         busco_out_dir = "analyses/genome_qc/set8/busco/by_taxon/", 
                         busco_db = "run_nostocales_odb10", 
                         out_file_name = "full_table.tsv") 
# Summarize busco output
busco_sum_set7 <- busco_out_set7 %>% 
  sum_busco()
busco_sum_set8 <- busco_out_set8 %>% 
  sum_busco()
# Rename columns by set
colnames(busco_sum_set7) <- str_c(colnames(busco_sum_set7), "_set7")
colnames(busco_sum_set8) <- str_c(colnames(busco_sum_set8), "_set8")
# join dfs
sum_busco_set7_set8 <- busco_sum_set7 %>%
  full_join(busco_sum_set8, by = c("sample_id_set7" = "sample_id_set8")) %>%
  mutate(complete_diff = 
           Complete_set7 - Complete_set8)
# Write table for genomes that had less complete buscos after graphbin
sum_busco_set7_set8 %>%
  filter(complete_diff > 0) %>%
  write.csv("document/report_2/busco_comparison.csv")

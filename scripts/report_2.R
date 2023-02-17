#!/usr/bin/env Rscript

source("scripts/r_functions.R")
library(tidyverse)
library(ape)

# Load genome ids for set 8 and make a df with sample (specimen) ids
genome_ids <- scan("scripts/genome_ids_set8", what = "character") 
genome_ids_df <- genome_ids %>%
  as_tibble() %>%
  rename(genome_id = value) %>%
  mutate(sample_id = 
           str_remove(genome_id, pattern = "_bin.*"))
# Load voucher primer table and merge with genome ids for set 8
voucher_primer <- read_delim("scripts/voucher_primer_2", delim = "\t") %>%
  right_join(genome_ids_df, by = "sample_id")
# Load checkm output table
checkm_out <- read_delim("analyses/genome_qc/set8/checkm/checkm_out") %>%
  mutate(`Bin Id` =
           paste(`Bin Id`, ".fa", sep = ""))
checkm_colnames <- stringr::str_replace(colnames(checkm_out), " ", "_") %>%
  stringr::str_to_lower() 
colnames(checkm_out) <- checkm_colnames
# Load gunc output table
gunc_out <- read_table("analyses/genome_qc/set8/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv") %>%
  mutate(genome =
           paste(genome, ".fa", sep = ""))
# Load and merge quast output
quast_out <- merge_quast(sample_ids = genome_ids, 
                         quast_out_dir = "analyses/genome_qc/set8/quast") %>%
  mutate(Assembly = paste(Assembly, ".fa", sep = ""))
# Load and merge busco output
busco_out <- merge_busco(sample_ids = genome_ids, 
                         busco_out_dir = "analyses/genome_qc/set8/busco/by_taxon/", 
                         busco_db = "run_nostocales_odb10", 
                         out_file_name = "full_table.tsv") 
# Summarize busco output
busco_sum <- busco_out %>% 
  sum_busco()
# Join qc tables
fullqc_set8 <- voucher_primer %>%
  left_join(checkm_out, by = c("genome_id" = "bin_id")) %>%
  left_join(gunc_out, by = c("genome_id" = "genome")) %>%
  left_join(quast_out, by = c("genome_id" = "Assembly")) %>%
  left_join(busco_sum, by = c("genome_id" = "sample_id"))
# Make report directory and write qc table
system("mkdir -p document/report_2")
write_delim(fullqc_set8, file = "document/report_2/fullqc_set8.csv", delim = ",")

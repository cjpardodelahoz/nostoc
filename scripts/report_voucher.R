#!/usr/bin/env Rscript

source("scripts/r_functions.R")
library(tidyverse)
library(ape)

# Load genome ids for set 12c and 12p and make a df with sample (specimen) ids
genome_ids <- scan("misc_files/genome_ids_set12c", what = "character")
plasmid_ids <- scan("misc_files/genome_ids_set12p", what = "character") %>%
  as_tibble() %>%
  mutate(genome_id = 
           str_remove(value, pattern = "_plasmid.fa")) %>%
  rename(plasmid_filename = value)
genome_ids_df <- genome_ids %>%
  as_tibble() %>%
  rename(chromosome_filename = value) %>%
  mutate(genome_id = 
           str_remove(chromosome_filename, pattern = "_chromosome.fa")) %>%
  left_join(plasmid_ids, by = "genome_id") %>%
  mutate(sample_id = 
           str_remove(genome_id, "_bin.*")) 
# Load voucher primer table and merge with genome ids for set 103a
voucher_df <- read_delim("misc_files/voucher_primer_3.txt", delim = "\t") %>%
  left_join(genome_ids_df, by = "sample_id") %>%
  relocate(sample_id, genome_id)
write_csv(voucher_df, file = "document/tables/voucher_v1.csv")


tree <- read.tree("analyses/phylogenetics/set103/divtime/1_part/mcmc/c1/dated.tree")
tree_key <-  voucher_df %>%
  filter(!is.na(genome_id)) %>%
  mutate(tree_tip = 
           paste(genome_id, ".fa", sep = "")) %>%
  mutate(taxon_name = 
           ifelse(genome_id == sample_id, 
                  yes = taxon_name, 
                  no = paste(genome_id, taxon_name, sep = "_"))) %>%
  select(tree_tip, taxon_name)
tree_new <- rename_taxa(tree, tree_key)
write.tree(tree_new, file = "document/sharing/dated_tree.tree")  

perrin_voucher <-  voucher_df %>%
  filter(!is.na(genome_id)) %>%
  mutate(tree_tip = 
           paste(genome_id, ".fa", sep = "")) %>%
  mutate(taxon_name = 
           ifelse(genome_id == sample_id, 
                  yes = taxon_name, 
                  no = paste(genome_id, taxon_name, sep = "_"))) %>%
  select(taxon_name, host_family, host_genus, substrate_type, lifestyle,
         country,  locality, voucher, substrate_notes, reference)
write_csv(perrin_voucher, 
          file = "document/tables/nostoc_metadata_for_perrin_v1.csv")





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

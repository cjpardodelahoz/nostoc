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
# Load and merge quast output
quast_out <- merge_quast(sample_ids = genome_ids, 
                         quast_out_dir = "analyses/genome_qc/set12c/quast") %>%
  mutate(Assembly = paste(Assembly, ".fa", sep = ""))
# Load and merge busco output
busco_out <- merge_busco(sample_ids = genome_ids, 
                         busco_out_dir = "analyses/genome_qc/set12c/busco/by_taxon/", 
                         busco_db = "run_nostocales_odb10", 
                         out_file_name = "full_table.tsv") 
# Summarize busco output
busco_sum <- busco_out %>% 
  sum_busco()
# Load lifestyle metadata
lifestyle_metadata <- read_csv("misc_files/set103_lifestyle_metadata.csv")
# Load table with plasmid classification info and get plasmid lengths and chromosome depths
depths_and_lenghts <- read.csv(file = "document/tables/contig_class_consensus.csv") %>%
  mutate(genome_id = str_remove(genome_id, ".fa")) %>%
  group_by(genome_id, contig_class) %>%
  summarise(median_depth = median(totalAvgDepth), length = sum(contig_length)) %>%
  ungroup(genome_id, contig_class) %>%
  pivot_wider(names_from = contig_class, values_from = c(length, median_depth)) %>%
  select(genome_id, median_depth_chromosome, length_plasmid) %>%
  rename(chromosome_median_depth = median_depth_chromosome,
         combined_plasmids_length = length_plasmid) %>%
  mutate(chromosome_median_depth = case_when(
    !is.na(chromosome_median_depth) ~ paste(round(chromosome_median_depth), "x", sep = ""),
    !is.na(chromosome_median_depth) ~ NA))
# Generate voucher table joining voucher primer 5, BUSCO, QUAST and Plasmid results
voucher_df <- read_delim("misc_files/voucher_primer_5.txt", delim = "\t") %>%
  right_join(genome_ids_df, by = "sample_id") %>%
  relocate(sample_id, genome_id) %>%
  left_join(quast_out, by = c("chromosome_filename" = "Assembly")) %>%
  left_join(busco_sum, by = c("chromosome_filename" = "sample_id")) %>%
  left_join(lifestyle_metadata, by = c("genome_id" = "tip_label")) %>%
  left_join(depths_and_lenghts, by = "genome_id") %>%
  mutate(taxon_name = case_when(
    str_detect(genome_id, taxon_name) | str_detect(taxon_name, genome_id) ~ taxon_name,
    .default = paste(genome_id, taxon_name, sep = "_")),
    percent_complete = round(percent_complete, digits = 1)) %>%
  select(taxon_name, voucher, site_id, lifestyle, lichen_substrate,
         N50, "Total length", "GC (%)", chromosome_median_depth, percent_complete, combined_plasmids_length, accession, 
         reference) %>%
  rename("Taxon name" = taxon_name,
         "Region/Voucher" = voucher,
         "Site ID" = site_id,
         "Lifestyle"= lifestyle,
         "Lichen substrate" = lichen_substrate,
         "Chromosome N50" = N50,
         "Chromosome total length (bp)" = "Total length",
         "Chromosome median depth" = chromosome_median_depth,
         "BUSCO completenes (%)" = percent_complete,
         "Plasmids total length (bp)" = combined_plasmids_length,
         "NCBI accession" = accession,
         "Reference" = reference)
# The genome labeled P6636 was mislabeled from the sequencing source and the actual code is P6236
# This is to correct in the final table but all previous files will retain P6636
voucher_df[107, 1] <- "P6236_bin_11_Peltigera_fuscopraetextata"
# This is to update the DNA ID of the ponojensis type per Jola's request
voucher_df[151, 1] <- "PL483_bin_4_Peltigera_ponojensis"
write_csv(voucher_df, file = "document/tables/voucher_v3.csv")

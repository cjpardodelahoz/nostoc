#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse)

# Load voucher tables and extract DNA ID, DNA source, rbcLX accession, and
# phylogroup
# Note that some of the special characters in the voucher column may not be 
# parsed correctly
obrien_2013 <- read_csv(file = "misc_files/voucher_obrien2013.csv", 
                        locale=locale(encoding="latin1")) %>%
  rename(dna_id = collection_number) %>%
  mutate(dna_source =paste(genus, species, sep = " ")) %>%
  rename(nostoc_phylogroup = rbclx_cluster) %>%
  add_column(reference = "obrien_et_al_2013") %>%
  select(dna_id, dna_source, rbclx_accession, nostoc_phylogroup, reference)
magain_2017 <- read_csv(file = "misc_files/voucher_magain2017.csv", 
                        locale=locale(encoding="latin1")) %>%
  rename(dna_source = Taxon) %>%
  rename(dna_id = "DNA Number") %>%
  rename(rbclx_accession = "rbcLX") %>%
  rename(nostoc_phylogroup = "rbcLX phylogroup or haplotype") %>%
  add_column(reference = "magain_et_al_2017") %>%
  select(dna_id, dna_source, rbclx_accession, nostoc_phylogroup, reference) %>%
  filter(rbclx_accession != "--") %>%
  mutate(dna_id = na_if(dna_id, "GB"))
miadlikowska_2018 <- read_csv(file = "misc_files/voucher_miadlikowska2018.csv",
                              locale=locale(encoding="latin1")) %>%
  rename(dna_id = "DNA extraction no.*") %>%
  rename(dna_source = "Species name and clade no.**") %>%
  rename(rbclx_accession = "rbcLX") %>%
  rename(nostoc_phylogroup = "Nostoc haplotype****") %>%
  add_column(reference = "miadlikowska_et_al_2018") %>%
  select(dna_id, dna_source, rbclx_accession, nostoc_phylogroup, reference) %>%
  filter(rbclx_accession != "missing" & !is.na(rbclx_accession))
chagnon_2018 <- read_csv(file = "misc_files/voucher_chagnon_2018.csv", 
                         locale=locale(encoding="latin1")) %>%
  rename(dna_id = "Collection Number") %>%
  rename(dna_source = "Species") %>%
  rename(rbclx_accession = "rbcLX accession number") %>%
  rename(nostoc_phylogroup = "rbcLX phylogroup") %>%
  add_column(reference = "chagnon_et_al_2018") %>%
  select(dna_id, dna_source, rbclx_accession, nostoc_phylogroup, reference) %>%
  filter(rbclx_accession != "XXXXXXXX" & rbclx_accession != "YYYYYYYY")
pardodelahoz_2018 <- read_csv(file = "misc_files/voucher_pardodelahoz2018.csv",
                              locale=locale(encoding="latin1")) %>%
  rename(dna_id = "DNA id.") %>%
  rename(dna_source = "Species") %>%
  rename(rbclx_accession = "rbcLX (Nostoc)") %>%
  rename(nostoc_phylogroup = "Nostoc phylogroup or haplotype") %>%
  add_column(reference = "pardodelahoz_et_al_2018") %>%
  select(dna_id, dna_source, rbclx_accession, nostoc_phylogroup, reference) %>%
  filter(rbclx_accession != "-")
magain_2018 <- read_csv(file = "misc_files/voucher_magain2018.csv",
                        locale=locale(encoding="latin1")) %>%
  rename(dna_source = "taxon") %>%
  rename(rbclx_accession = "rbclx") %>%
  rename(nostoc_phylogroup = "rbclx_phylogroup_haplotype") %>%
  add_column(reference = "magain_et_al_2018") %>%
  select(dna_id, dna_source, rbclx_accession, nostoc_phylogroup, reference) %>%
  filter(rbclx_accession != "---")
# Join all rbclx metadata tables
public_rbclx_metadata <- bind_rows(obrien_2013, magain_2017, miadlikowska_2018, 
                                   chagnon_2018, pardodelahoz_2018, magain_2018) %>%
  distinct(rbclx_accession, .keep_all = T)
# Save the list of unique accessions and the metadata
dir.create(path = "analyses/species_delimitation/rbclx/public", recursive = T)
public_rbclx_metadata %>%
  pull(rbclx_accession) %>%
  write(file = "analyses/species_delimitation/rbclx/public/public_rbclx_accessions.txt")
write_csv(public_rbclx_metadata, file = "analyses/species_delimitation/rbclx/public/public_rbclx_metadata.csv")

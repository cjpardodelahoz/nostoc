#!/usr/bin/env Rscript

##### Load required packages and functions ####
library(tidyverse)
# Function to extract part before the first semicolon
extract_before_semicolon <- function(string) {
  parts <- unlist(strsplit(string, ";")) # Split the string
  trimmed_part <- trimws(parts[1])       # Get the first part and trim whitespace
  return(trimmed_part)
}
extract_before_semicolon_v <- Vectorize(extract_before_semicolon, USE.NAMES = F)
# Function to extract part before the first two commas
extract_before_commas <- function(string) {
  parts <- unlist(strsplit(string, ",")) # Split the string
  trimmed_part_1 <- trimws(parts[1])       # Get the first part and trim whitespace
  trimmed_part_2 <- trimws(parts[2])
  result <- paste(trimmed_part_1, trimmed_part_2, sep = ", ")
  return(result)
}
extract_before_commas_v <- Vectorize(extract_before_commas, USE.NAMES = F)
# Function to extract part before the first comma
extract_before_comma <- function(string) {
  parts <- unlist(strsplit(string, ",")) # Split the string
  trimmed_part <- trimws(parts[1])       # Get the first part and trim whitespace
  return(trimmed_part)
}
extract_before_comma_v <- Vectorize(extract_before_comma, USE.NAMES = F)


#### Curate and join metadata tables ####

# Load voucher tables and extract DNA ID, DNA source, rbcLX accession, and
# phylogroup
# Note that some of the special characters in the voucher column may not be 
# parsed correctly
obrien_2013 <- read_csv(file = "misc_files/voucher_obrien2013.csv", 
                        locale=locale(encoding="latin1")) %>%
  rename(dna_id = collection_number) %>%
  mutate(dna_source = paste(genus, species, sep = " ")) %>%
  mutate(dna_source = str_replace(dna_source, "Peltigera", "P.")) %>%
  mutate(dna_source = str_remove(dna_source, " Clade.*")) %>%
  rename(nostoc_phylogroup = rbclx_cluster) %>%
  add_column(reference = "obrien_et_al_2013") %>%
  add_column(region = "Canada, British Columbia") %>%
  select(dna_id, dna_source, rbclx_accession, nostoc_phylogroup, region, reference)
magain_2017 <- read_csv(file = "misc_files/voucher_magain2017.csv", 
                        locale=locale(encoding="latin1")) %>%
  rename(dna_source = Taxon) %>%
  rename(dna_id = "DNA Number") %>%
  rename(rbclx_accession = "rbcLX") %>%
  rename(nostoc_phylogroup = "rbcLX phylogroup or haplotype") %>%
  add_column(reference = "magain_et_al_2017") %>%
  mutate(region = extract_before_semicolon_v(`Voucher/Published source`)) %>%
  mutate(region = str_replace(region, "Qubec", "Quebec")) %>%
  mutate(region = str_replace(region, "Iceland.*", "Iceland")) %>%
  mutate(region = str_remove(region, ", Rudi et al. 1998")) %>%
  select(dna_id, dna_source, rbclx_accession, nostoc_phylogroup, region, reference) %>%
  filter(rbclx_accession != "--") %>%
  mutate(dna_id = na_if(dna_id, "GB"))
miadlikowska_2018 <- read_csv(file = "misc_files/voucher_miadlikowska2018.csv",
                              locale=locale(encoding="latin1")) %>%
  rename(dna_id = "DNA extraction no.*") %>%
  rename(dna_source = "Species name and clade no.**") %>%
  rename(rbclx_accession = "rbcLX") %>%
  rename(nostoc_phylogroup = "Nostoc haplotype****") %>%
  add_column(reference = "miadlikowska_et_al_2018") %>%
  mutate(region = extract_before_commas_v(Voucher)) %>%
  mutate(region = str_replace(region, "Qubec", "Quebec")) %>%
  mutate(region = str_remove(region, " J. Miadliko.*")) %>%
  mutate(region = str_remove(region, ", N. Magain s.n.")) %>%
  mutate(region = str_remove(region, ", H. Kristinsson L-27478")) %>%
  mutate(region = str_remove(region, ", A. Agren 609")) %>%
  mutate(region = str_replace(region, "U.S.", "USA")) %>%
  select(dna_id, dna_source, rbclx_accession, nostoc_phylogroup, region, reference) %>%
  filter(rbclx_accession != "missing" & !is.na(rbclx_accession))
chagnon_2018 <- read_csv(file = "misc_files/voucher_chagnon_2018.csv", 
                         locale=locale(encoding="latin1")) %>%
  rename(dna_id = "Collection Number") %>%
  rename(dna_source = "Species") %>%
  rename(rbclx_accession = "rbcLX accession number") %>%
  rename(nostoc_phylogroup = "rbcLX phylogroup") %>%
  add_column(reference = "chagnon_et_al_2018") %>%
  mutate(region = extract_before_semicolon_v(`Voucher Information / Collection Site`)) %>%
  select(dna_id, dna_source, rbclx_accession, nostoc_phylogroup, region, reference) %>%
  filter(rbclx_accession != "XXXXXXXX" & rbclx_accession != "YYYYYYYY")
pardodelahoz_2018 <- read_csv(file = "misc_files/voucher_pardodelahoz2018.csv",
                              locale=locale(encoding="latin1")) %>%
  rename(dna_id = "DNA id.") %>%
  rename(dna_source = "Species") %>%
  rename(rbclx_accession = "rbcLX (Nostoc)") %>%
  rename(nostoc_phylogroup = "Nostoc phylogroup or haplotype") %>%
  add_column(reference = "pardodelahoz_et_al_2018") %>%
  mutate(region = extract_before_commas_v(Voucher)) %>%
  mutate(region = str_replace(region, "U.S.", "USA")) %>%
  mutate(region = str_replace(region, "columbia", "Columbia")) %>%
  mutate(region = str_replace(region, "BC", "British Columbia")) %>%
  mutate(region = str_replace(region, "Canada, Qu.*", "Canada, Quebec")) %>%
  mutate(region = str_remove(region, ", B. Krzewicka 2512")) %>%
  mutate(region = str_remove(region, ", I. Karlsson s.n.")) %>%
  mutate(region = str_remove(region, ", N. Magain s.n.")) %>%
  mutate(region = str_replace(region, "kray", "Kray")) %>%
  mutate(region = str_replace(region, "Iceland, Langisj.*",  "Iceland, Langisjor")) %>%
  mutate(region = str_replace(region, "ColombiaC", "Columbia")) %>%
  select(dna_id, dna_source, rbclx_accession, nostoc_phylogroup, region, reference) %>%
  filter(rbclx_accession != "-")
magain_2018 <- read_csv(file = "misc_files/voucher_magain2018.csv",
                        locale=locale(encoding="latin1")) %>%
  rename(dna_source = "taxon") %>%
  rename(rbclx_accession = "rbclx") %>%
  rename(nostoc_phylogroup = "rbclx_phylogroup_haplotype") %>%
  add_column(reference = "magain_et_al_2018") %>%
  mutate(region = extract_before_comma_v(`voucher_published source`)) %>%
  mutate(region = str_replace(region, ":", ",")) %>%
  mutate(region = str_remove(region, "; Truong 3991")) %>%
  select(dna_id, dna_source, rbclx_accession, nostoc_phylogroup, region, reference) %>%
  filter(rbclx_accession != "---")
# Join all rbclx metadata tables
public_rbclx_metadata <- bind_rows(obrien_2013, magain_2017, miadlikowska_2018, 
                                   chagnon_2018, pardodelahoz_2018, magain_2018) %>%
  distinct(rbclx_accession, .keep_all = T) %>%
  mutate(dna_source = str_remove(dna_source, " Gyeln.*")) %>%
  mutate(dna_source = str_remove(dna_source, " Zahlbr.*")) %>%
  mutate(dna_source = str_remove(dna_source, " S\u008erus.*")) %>%
  mutate(dna_source = str_remove(dna_source, " C. W. Dodge")) %>%
  mutate(dna_source = str_remove(dna_source, " \\(With.\\).*")) %>%
  mutate(dna_source = str_remove(dna_source, " R\u008as.*")) %>%
  mutate(dna_source = str_remove(dna_source, " Goward.*")) %>%
  mutate(dna_source = str_remove(dna_source, " Vitik.*")) %>%
  mutate(dna_source = str_remove(dna_source, " L. F.*")) %>%
  mutate(dna_source = str_remove(dna_source, " \\(.*")) %>%
  mutate(dna_source = str_remove(dna_source, " M\u009f.*")) %>%
  mutate(dna_source = str_remove(dna_source, " Vain.*")) %>%
  mutate(dna_source = str_remove(dna_source, " R. Sa.*"))
# Save the list of unique accessions and the metadata
dir.create(path = "analyses/species_delimitation/rbclx/public", recursive = T)
public_rbclx_metadata %>%
  pull(rbclx_accession) %>%
  write(file = "analyses/species_delimitation/rbclx/public/public_rbclx_accessions.txt")
write_csv(public_rbclx_metadata, file = "analyses/species_delimitation/rbclx/public/public_rbclx_metadata.csv")

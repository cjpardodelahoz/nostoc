#!/usr/bin/env Rscript

source("scripts/r_functions.R")
library(tidyverse)
library(ape)

##### GENOME GENBANK SUMBMISSION #####

# Metadata for GoLife and ABMI sites
golife_site_data <- read_csv("document/tables/golife_site_metadata.csv") %>%
  rename("Elevation" = "Elevation (mASL)")
abmi_site_data <- read_csv("document/tables/abmi_site_data.csv") %>%
  rename("Site ID" = "Site",
         "Country" = "Region",
         "Date" = "Year",
         "Latitude" = "Public latitud",
         "Longitude" = "Public longitude",
         "Elevation" = "Altitude (mASL)")
# Table with information to submit lichen biosamples
lichen_biosample_info <- read_delim("misc_files/voucher_primer_5.txt", delim = "\t") %>%
  right_join(genome_ids_df, by = "sample_id") %>%
  relocate(sample_id, genome_id) %>%
  filter(str_detect(reference, "study") | is.na(reference)) %>%
  left_join(golife_site_data, by = c("site_id" = "Site ID")) %>%
  left_join(abmi_site_data, by = c("site_id" = "Site ID")) %>%
  mutate(Date.y = as.character(Date.y),
         Longitude.x = as.character(Longitude.x),
         geo_loc_name = case_when(is.na(Country.x) & !is.na(Country.y) ~ Country.y,
                                  !is.na(Country.x) & is.na(Country.y) ~ Country.x),
         collection_date = case_when(is.na(Date.x) & !is.na(Date.y) ~ Date.y,
                                     !is.na(Date.x) & is.na(Date.y) ~ Date.x),
         env_broad_scale = case_when(is.na(Habitat.x) & !is.na(Habitat.y) ~ Habitat.y,
                                     !is.na(Habitat.x) & is.na(Habitat.y) ~ Habitat.x),
         env_medium = "Lichen thallus",
         lat =  case_when(is.na(Latitude.x) & !is.na(Latitude.y) ~ Latitude.y,
                          !is.na(Latitude.x) & is.na(Latitude.y) ~ Latitude.x),
         lon = case_when(is.na(Longitude.x) & !is.na(Longitude.y) ~ Longitude.y,
                         !is.na(Longitude.x) & is.na(Longitude.y) ~ Longitude.x),
         lat_lon = paste(lat, lon, sep = "_")) %>%
  select(sample_id, taxon_name, collection_date, env_broad_scale, env_medium, geo_loc_name, lat_lon)
write_csv(biosample_info, "document/ncbi/lichen_biosample_info.csv")

# Table with MAG biosample info
lichen_biosample_accessions <- read_delim("document/ncbi/lichen_biosample_attributes.tsv")
mag_biosample_info <- read_delim("misc_files/voucher_primer_5.txt", delim = "\t") %>%
  right_join(genome_ids_df, by = "sample_id") %>%
  relocate(sample_id, genome_id) 
# Update the sample names where needed
mag_biosample_info[107, 1] <- "P6236"
mag_biosample_info[151, 1] <- "PL483"
mag_biosample_info[79, 1] <- "P14320"
mag_biosample_info[107, 2] <- "P6236_bin_11"
mag_biosample_info[151, 2] <- "PL483_bin_4"
mag_biosample_info[79, 2] <- "P14320_bin_11"
# Join lichen biosample info and accessions
mag_biosample_info <- mag_biosample_info %>%
  right_join(lichen_biosample_accessions, by = c("sample_id" = "sample_name")) %>%
  mutate("derived-from" = paste("This BioSample is a metagenomic assembly obtained from the lichen metagenome BioSample: ", accession.y, sep = ""),
         bioproject_accession = "PRJNA1066398",
         organism = "Nostoc sp.",
         sample_title = paste("Nostoc sp. cyanobiont of", sample_id, isolation_source, sep = " "),
         sample_name = genome_id) %>%
  select(sample_name, sample_title, bioproject_accession, organism, host, 
         isolation_source, collection_date, geo_loc_name, lat_lon, source_material_id, 
         env_broad_scale, "derived-from")
write_csv(mag_biosample_info, "document/ncbi/mag_biosample_info.csv")

# Table with information for SRA submission
sra_attributes <- lichen_biosample_accessions %>%
  mutate(title = paste("Metagenome of lichen ", sample_name, " ", isolation_source),
         filename = paste(sample_name, "_R1_all.fastq"),
         filename2 = paste(sample_name, "_R2_all.fastq")) %>%
  rename(biosample_accession = accession, 
         library_ID = sample_name) %>%
  select(biosample_accession, library_ID, title, filename, filename2)
# Update file names
sra_attributes[40, 4] <- "P14321_R1_all.fastq"
sra_attributes[40, 5] <- "P14320_R2_all.fastq"
sra_attributes[68, 4] <- "P6636_R1_all.fastq"
sra_attributes[68, 5] <- "P6636_R2_all.fastq"
sra_attributes[112, 4] <- "X2_R1_all.fastq"
sra_attributes[112, 5] <- "X2_R2_all.fastq"
write_csv(sra_attributes, "document/ncbi/sra_attributes.csv")
  

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
# Table with MAG biosample metadata for MAG file submission
mag_attributes <- read_delim(file = "document/ncbi/mag_biosample_attributes.tsv") %>%
  left_join(depths_and_lenghts, by = c("sample_name" = "genome_id")) %>%
  mutate(filename = paste(sample_name, "_chromosome.fa", sep = "")) %>%
  add_column(assembly_method = "metaSPAdes",
             assembly_method_version = "3.14.1",
             sequencing_technology = "Illumina NovaSeq 6000") %>%
  rename(biosample_accession = accession,
         genome_coverage = chromosome_median_depth) %>%
  select(biosample_accession, sample_name, assembly_method, assembly_method_version, genome_coverage, sequencing_technology, filename)
# Update the coverage for the genomes for which the name was updated after the analyses
mag_attributes[41, 5] <- "156x"
mag_attributes[69, 5] <- "99x"
mag_attributes[113, 5] <- "363x"
# Update the filenames for those genomes
mag_attributes[41, 7] <- "P14321_bin_11_chromosome.fa"
mag_attributes[69, 7] <- "P6636_bin_11_chromosome.fa"
mag_attributes[113, 7] <- "X2_bin_4_chromosome.fa"
# Write table with MAG metadata for submission
write_csv(mag_attributes, "document/ncbi/mag_attributes.csv")

##### VOUCHER TABKE ####

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
# This is to update the DNA ID of the ponojensis type per Jola's request
voucher_df[79, 1] <- "P14320_bin_11_Leptogium_sp"
write_csv(voucher_df, file = "document/tables/voucher_v3.csv")



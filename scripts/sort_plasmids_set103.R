#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse)
library(data.table)
library(ape)
source("scripts/r_functions.R")

#### EXPLORE PLASX RESULTS ####

# List of paths to plasx scores
scores_paths <- list.files(path = "analyses/plasmid_detection/set103" , 
                           pattern = "scores.txt", full.names = T, recursive = T)
header_key_paths <- list.files(path = "analyses/plasmid_detection/set103", 
                               pattern = "header_key.txt", full.names = T, 
                               recursive = T)
depths_paths <- list.files(path = "analyses/plasmid_detection/set103/depths", 
                               pattern = "depths.txt", full.names = T, 
                               recursive = T)
# Load plasx scores and header keys
# Note that plasx did not score all contigs, so the header_key files have more 
# entries than the scores files
header_key <- readr::read_delim(file = header_key_paths, delim = "\t", col_names = F, id = "header_key_path")
scores <- readr::read_delim(file = scores_paths, delim = "\t", col_names = T, id = "scores_path")
# Wrangle header_key and scores to create a df with plasx results
# Note the length and coverage values will only be obtained for the Nostoc genomes
# from set10 (i.e., the ones sequenced as part of this study)
scores <- scores %>%
  # extract genome id
  dplyr::mutate(genome_id = 
           stringr::str_remove(scores_path, ".*set103/")) %>%
  dplyr::mutate(genome_id =
           stringr::str_remove(genome_id, "/scores.txt"))
plasx_df <- header_key %>%
  dplyr::rename("contig_code" = "X1", "contig_label" = "X2") %>%
  dplyr::mutate(pwd = 
           stringr::str_remove(header_key_path, "/header_key.txt")) %>%
  # extract contig length
  dplyr::mutate(contig_length =
           stringr::str_remove(contig_label, "_cov.*")) %>%
  dplyr::mutate(contig_length = 
           stringr::str_remove(contig_length, ".*length_")) %>% 
  dplyr::mutate(contig_length =
           as.numeric(contig_length)) %>%
  # extract kmer coverage
  dplyr::mutate(contig_kmer_coverage =
           stringr::str_remove(contig_label, ".*cov_")) %>% 
  dplyr::mutate(contig_kmer_coverage =
           as.numeric(contig_kmer_coverage)) %>%
  dplyr::mutate(genome_id = 
           stringr::str_remove(header_key_path, ".*set103/")) %>%
  dplyr::mutate(genome_id =
           stringr::str_remove(genome_id, "/header_key.txt")) %>%
  # join plasx scores
  dplyr::left_join(scores, by = c("contig_code" = "contig", "genome_id")) %>% 
  # Classify contigs into plasmids and chromosomes with plasx score > 0.5
  dplyr::mutate(plasx_class =
           dplyr::if_else(score > 0.1, 
                          true = "plasmid", false = "chromosome"))
# Distribution of plasx scores
hist_plasx_scores <- plasx_df %>%
  dplyr::pull(score) %>%
  hist(main = "Distribution of plasx scores",
       xlab = "plasx score", 
       xlim = c(0, 1), breaks = 200)
# Check out which contig lengths don't get scored by plasx
hist_unscored_length <- plasx_df %>% 
  dplyr::filter(is.na(score)) %>%
  dplyr::pull(contig_length) %>%
  hist(breaks = 100, main = "Size distribution of unscored contigs", 
       xlab = "contig length (bp)")

#### EXPLORE DEEPLASMID RESULTS ####

# List of paths to deeplasmid predictions and scores
predictions_paths <- list.files(path = "analyses/plasmid_detection/set103" , 
                           pattern = "predictions.txt", full.names = T, recursive = T)
# Load deeplasmid predictions and scores
deeplasmid_df <- readr::read_delim(file = predictions_paths, delim = ",", 
                            col_names = T, id = "predictions_path") %>%
  # Remove the +/- from the deeplasmid scores
  dplyr::mutate(deeplasmid_score = 
                  as.numeric(str_remove(conf, " .*"))) %>%
  # Get the genome id of the contigs
  dplyr::mutate(genome_id = 
                  str_remove(predictions_path, "analyses/plasmid_detection/set103/")) %>%
  dplyr::mutate(genome_id =
                  str_remove(genome_id, "/outPR.*"))
# Distribution of deeplasmid scores
hist_deeplasmid_scores <- deeplasmid_df %>%
  dplyr::pull(deeplasmid_score) %>%
  hist(main = "Distribution of deeplasmid scores",
       xlab = "deeplasmid score", 
       xlim = c(0, 1), breaks = 200) 

#### CONTIG CLASSIFICATION WITH DEEPLASMID AND PLASX ####

# Read and bind the depth tables
depth_tables <- list()
for (i in 1:length(depths_paths)) {
  depth_tables[[depths_paths[i]]] <- read_delim(file = depths_paths[i], 
                                                delim = "\t", col_names = T) %>%
    select(1:3)
}
depths <- bind_rows(depth_tables, .id = "depth_path") %>%
  mutate(genome_id =
           stringr::str_remove(depth_path, "analyses/plasmid_detection/set103/depths/")) %>%
  mutate(genome_id =
           stringr::str_remove(genome_id, "_depths.txt")) %>%
  mutate(contigName = 
           str_remove(contigName, "_cov.*")) %>%
  select(2:5)
# Join the plasx_df with the deeplasmid results and the contig depths for
# classification
contig_class <- left_join(plasx_df, deeplasmid_df, 
                          by = c("contig_code" = "name", "genome_id" = "genome_id")) %>%
  mutate(contig_label_trimmed = 
           str_remove(contig_label, "_cov.*")) %>%
  left_join(depths, by = c("genome_id" = "genome_id", "contig_label_trimmed" = "contigName")) %>%
  select(contig_code, contig_label, genome_id, contig_length, contig_kmer_coverage, 
         totalAvgDepth, score, deeplasmid_score, plasx_class, pred)
# Get table with summary of median contig depth per genome
depth_sum <- contig_class %>% 
  group_by(genome_id) %>%
  filter(plasx_class == "chromosome") %>%
  summarise(median_depth = median(totalAvgDepth, na.rm = T))
# Classify contigs using results from both plasx and deeplasmid. For contigs
# Deeplasmid identified as plasmid but plasx didn't, we'll only consider them
# as plasmid if they have aberrant coverage, defined as deviating by >20X from the
# median coverage of chromosome contigs according to plasx
contig_class_consensus <- contig_class %>%
  left_join(depth_sum, by = "genome_id") %>%
  mutate(contig_class = 
           if_else(condition = plasx_class == "chromosome" & pred == "PLASMID" &
                     abs(totalAvgDepth - median_depth) > 20, 
                   true = "plasmid", false = plasx_class))

#### How much sequence length per genome is plasmid? ####

# Summarise contig origin length
contig_origin_summary <- plasx_df %>%
  tidyr::pivot_wider(names_from = contig_origin, 
                     values_from = contig_length) %>%
  dplyr::rename(unassigned = "NA") %>%
  dplyr::group_by(genome_id) %>%
  dplyr::summarise(chromosome_length = sum(chromosome, na.rm = T), 
            plasmid_length = sum(plasmid, na.rm = T), 
            unassigned_length = sum(unassigned, na.rm = T))

plasx_df %>%
  filter(genome_id == "JL33_bin_16.fa" & contig_origin == "chromosome") %>%
  pull(contig_kmer_coverage) %>%
  hist(breaks = 40)

#### Sort contig names by origin for genome refinement ####

# Get df for plasmid and chromosome contig sorting
contig_class_sum <- contig_class_consensus %>%
  select(contig_label, contig_class, genome_id) %>%
  mutate(genome_id = 
           stringr::str_remove(genome_id, ".fa"))
# Vector with contig paths. I used the files with the original contig names to
# avoid losing information present in the headers
contig_paths <- list.files(path = "analyses/cyano_genomes/set103", 
                           pattern = ".fa", recursive = T, full.names = T)
# Vector with genome ids for sorting (without ".fa")
genome_ids <- contig_class_sum %>%
  pull(genome_id) %>%
  unique()
# Sort contigs
bin_from_df_v2(sample_id = genome_ids, 
               contigs_file = contig_paths, 
               bin_df = contig_class_sum, 
               out_dir = "analyses/plasmid_detection/set103",
               out_suffix_1 = ".fa")

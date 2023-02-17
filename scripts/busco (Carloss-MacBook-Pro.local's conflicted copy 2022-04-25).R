# Load required libraries and custom functions
source("scripts/r_functions.R")
library(tidyverse)

# Load busco output for set 1
sample_ids <- c("3_bin.5.fa", "8274_bin.7.fa")
busco_test <- sum_busco(sample_ids = sample_ids, 
                        busco_out_dir = "analyses/genome_qc/set1/busco/by_taxon", 
                        busco_db = "run_cyanobacteria_odb10", 
                        out_file_name = "full_table.tsv")
# Get paths for busco loci
path_test_1 <- busco_test %>%
  filter(status == "Complete") %>%
  mutate(na_seq_path =
           paste("analyses/genome_qc/set1/busco/by_taxon", "/", sample_id, 
                 "/", "run_cyanobacteria_odb10", "/", 
                 "busco_sequences/single_copy_busco_sequences", "/",
                 busco_id, ".fna", sep = "")
         ) %>%
  pull(na_seq_path)
paths <- path_test_1[1:3]
seqs_test <- read.multiFASTA(paths)
write.FASTA(seqs_test, file = "test.fa")


analyses/genome_qc/set1/busco/by_taxon/3_bin.5.fa/run_cyanobacteria_odb10/busco_sequences/single_copy_sequences/105at1117.fna

analyses/genome_qc/set1/busco/by_taxon/3_bin.5.fa/run_cyanobacteria_odb10/busco_sequences/single_copy_busco_sequences/105at1117.fna

seq_base_path <- c("busco_sequences/single_copy_")

busco_test <- read_busco(sample_ids = "3_bin.5.fa", 
                         busco_out_dir = "analyses/genome_qc/set1/busco/by_taxon", 
                         busco_db = "run_cyanobacteria_odb10", 
                         out_file_name = "full_table.tsv") 
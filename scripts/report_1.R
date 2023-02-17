#!/usr/bin/env Rscript

library(tidyverse)
library(ape)
library(treeio)

# Load voucher primer table
voucher_primer_1 <- read_delim("scripts/voucher_primer_1", delim = "\t")
# Load checkm output table
checkm_set1_out <- read_delim("analyses/genome_qc/set1/checkm/checkm_out") %>%
  mutate(`Bin Id` =
           paste(`Bin Id`, ".fa", sep = "")
         )
checkm_colnames <- stringr::str_replace(colnames(checkm_set1_out), " ", "_") %>%
  stringr::str_to_lower() 
colnames(checkm_set1_out) <- checkm_colnames
# Load gunc output table
gunc_set1_out <- read_table("analyses/genome_qc/set1/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv") %>%
  mutate(genome =
           paste(genome, ".fa", sep = ""))
# Join qc tables and write table
fullqc_set1 <- voucher_primer_1 %>%
  right_join(checkm_set1_out, by = c("sample_id" = "bin_id")) %>%
  left_join(gunc_set1_out, by = c("sample_id" = "genome"))
write_delim(fullqc_set1, file = "document/report_1/fullqc_set1.csv", delim = ",")
# Load na concatenated tree of set1 taxa
concat_set1_ng_na <- read.tree("analyses/phylogenetics/set1/trees/concat/concat_ng_na1.treefile")
# Change tip labels to taxa names
concat_set1_ng_na <- rename_taxa(concat_set1_ng_na, voucher_primer_1, key = 1, value = 3)
# Save updated tree
write.tree(concat_set1_ng_na, file = "document/report_1/concat_set1_ng_na.tree")

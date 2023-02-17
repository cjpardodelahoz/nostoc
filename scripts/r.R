# Load custom functions
source("scripts/r_functions.R")
# Load list of taxa
sample_ids_6240_6256_7535 <- scan(file = "scripts/sample_ids_6240_6526_7535", 
                                  what = "character")
# Combine checkm output for all bins from orders 6240_6256_7535
checkm_6240_6256_7535 <- sum_checkm(sample_ids_6240_6256_7535, 
                                    checkm_out_dir = "analyses/genome_qc/metabat/6240_6526_7535",
                                    out_file_name = "checkm_out",
                                    bin_dir = "analyses/bins/metabat")
# Filter checkm df to cyano bins
checkm_cyano_6240_6256_7535 <- checkm_6240_6256_7535 %>%
  dplyr::filter(
    stringr::str_detect(string = marker_lineage, pattern = "Cyanobacteria*.")
  )
# Test which samples have multiple cyano bins
multi_cyano_6240_6256_7535 <- checkm_cyano_6240_6256_7535 %>%
  pull(sample_id) %>% duplicated()
# Add column indicating samples that have multiple cyanos
checkm_cyano_6240_6256_7535 <- checkm_cyano_6240_6256_7535 %>%
  dplyr::add_column(multi = multi_cyano_6240_6256_7535)
# Show records with multiple cyano bins
checkm_cyano_6240_6256_7535 %>%
  filter(multi == TRUE)

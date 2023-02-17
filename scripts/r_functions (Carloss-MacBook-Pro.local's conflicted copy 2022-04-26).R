# Function to read the checkm output from a single genome in a run
read_checkm <- function(sample_id, checkm_out_dir, out_file_name) {
  checkm_out_file <- paste(checkm_out_dir, sample_id, out_file_name, sep = "/")
  checkm_out <- readr::read_delim(checkm_out_file, delim = "\t") %>%
    dplyr::mutate(sample_id = sample_id)
  return(checkm_out)
}

# Function to combine the checkm output from multiple genomes of the same run into
# a single table
sum_checkm <- function(sample_ids, checkm_out_dir, out_file_name, bin_dir, bin_ext) {
  checkm_out_list <- lapply(sample_ids, read_checkm, 
                            checkm_out_dir = checkm_out_dir, 
                            out_file_name = out_file_name)
  checkm_sum <- dplyr::bind_rows(checkm_out_list)
  new_colnames <- stringr::str_replace(colnames(checkm_sum), " ", "_") %>%
    stringr::str_to_lower()
  colnames(checkm_sum) <- new_colnames
  checkm_sum <- checkm_sum %>%
    dplyr::mutate(bin_path = 
                    paste(bin_dir, "/", sample_id, "/", bin_id, 
                          bin_ext, sep = "")
    )
  return(checkm_sum)
}

# Function to read the busco output from a single genome in a run
read_busco <- function(sample_id, busco_out_dir, busco_db, out_file_name) {
  busco_out_file <- paste(busco_out_dir, sample_id, busco_db, out_file_name, sep = "/")
  busco_out <- readr::read_delim(busco_out_file, skip = 2, delim = "\t") %>%
    dplyr::mutate(sample_id = sample_id)
  return(busco_out)
}

# Function to combine the busco output from multiple genomes of the same run into
# a single table
sum_busco <- function(sample_ids, busco_out_dir, busco_db, out_file_name) {
  busco_out_list <- lapply(sample_ids, read_busco, 
                            busco_out_dir = busco_out_dir, 
                            out_file_name = out_file_name,
                            busco_db = busco_db)
  busco_sum <- dplyr::bind_rows(busco_out_list)
  new_colnames <- stringr::str_remove(colnames(busco_sum), "# ") %>%
    stringr::str_replace(" ", "_") %>%
    stringr::str_to_lower()
  colnames(busco_sum) <- new_colnames
  return(busco_sum)
}

# Function to add paths to complete single copy busco to busco df
add_busco_paths <- function(busco_df, busco_out_dir, busco_db) {
  seq_base_path <- c("busco_sequences/single_copy_busco_sequences")
  busco_new_df <- busco_df %>%
  filter(status == "Complete") %>%
    mutate(na_seq_path =
             paste(busco_out_dir, "/", sample_id, 
                   "/", busco_db, "/", 
                   seq_base_path, "/",
                   busco_id, ".fna", sep = "")
           ) %>%
    mutate(aa_seq_path =
             paste(busco_out_dir, "/", sample_id, 
                   "/", busco_db, "/", 
                   seq_base_path, "/",
                   busco_id, ".faa", sep = "")
           )
  return(busco_new_df)
}


path_test_1 <- busco_test %>%
  filter(status == "Complete") %>%
  mutate(na_seq_path =
           paste("analyses/genome_qc/set1/busco/by_taxon", "/", sample_id, 
                 "/", "run_cyanobacteria_odb10", "/", 
                 "busco_sequences/single_copy_busco_sequences", "/",
                 busco_id, ".fna", sep = "")
  ) %>%
  pull(na_seq_path)
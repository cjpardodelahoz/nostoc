#!/usr/bin/env Rscript

# Operators
`%nin%` = Negate(`%in%`) # not in

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
# sample_id      Name of the genome that matches the output directory of the BUSCO run.
# busco_out_dir  Base path to directory with genome-specific BUSCO output
# busco_db       Name of the BUSCO database used in the run. This is used to build the path to the results E.g., "run_cyanobacteria_odb1o"
# out_file_name  Name of the BUSCO output talbe file. Typically "full_table.tsv".
read_busco <- function(sample_id, busco_out_dir, busco_db, out_file_name) {
  busco_out_file <- paste(busco_out_dir, sample_id, busco_db, out_file_name, sep = "/")
  busco_out <- readr::read_delim(busco_out_file, skip = 2, delim = "\t") %>%
    dplyr::mutate(sample_id = sample_id)
  return(busco_out)
}

# Function to combine the busco output from multiple genomes of the same run into
# a single table
# sample_ids      Character vector with names of the genomes that match the output directory of the BUSCO run.
# busco_out_dir   Base path to directory with genome-specific BUSCO output
# busco_db        Name of the BUSCO database used in the runs. This is used to build the path to the results E.g., "run_cyanobacteria_odb1o"
# out_file_name   Name of the BUSCO output table file. Typically "full_table.tsv".
merge_busco <- function(sample_ids, busco_out_dir, busco_db, out_file_name) {
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
# busco_df        Data frame with busco results as obtained with merge_busco()
# busco_out_dir   Base path to directory with genome-specific BUSCO output
# busco_db        Name of the BUSCO database used in the runs. This is used to build the path to the results E.g., "run_cyanobacteria_odb1o"
add_busco_paths <- function(busco_df, busco_out_dir, busco_db) {
  seq_base_path <- c("busco_sequences/single_copy_busco_sequences")
  busco_df_paths <- busco_df %>%
  filter(status == "Complete") %>%
    dplyr::mutate(na_seq_path =
             paste(busco_out_dir, "/", sample_id, 
                   "/", busco_db, "/", 
                   seq_base_path, "/",
                   busco_id, ".fna", sep = "")
           ) %>%
    dplyr::mutate(aa_seq_path =
             paste(busco_out_dir, "/", sample_id, 
                   "/", busco_db, "/", 
                   seq_base_path, "/",
                   busco_id, ".faa", sep = "")
           )
  return(busco_df_paths)
}

# Function to add seq labels to the filtered busco df
# busco_df_paths  Data frame with busco results after paths have been added with add_busco_paths()
add_busco_labels <- function(busco_df_paths) {
  paths <- busco_df_paths %>%
    pull(na_seq_path)
  seq_label <- lapply(paths, ape::read.FASTA) %>%
    lapply(., names) %>%
    unlist()
  busco_df_paths_labels <- busco_df_paths %>%
    tibble::add_column(seq_label)
  return(busco_df_paths_labels)
}

# Function to generate multifasta file for one busco locus
# busco                   ID of a single busco locus
# busco_df_paths_labels   Data frame with busco results after paths and labels have been added with add_busco_paths() and add_busco_labels()
# data_type               "DNA" or "AA"
# out_dir                 Path to directory to write the multifasta file
cat_busco_seqs <- function(busco, busco_df_paths_labels, data_type, out_dir) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = T)
  }
  filtered_df <- busco_df_paths_labels %>%
    dplyr::filter(busco_id == busco)
  ref_table <- filtered_df %>%
    dplyr::select(seq_label, sample_id)
  if (data_type == "DNA") {
    path_file <- paste(out_dir, "/", busco, "_na_paths", sep = "")
    tmp_file <- paste(out_dir, "/", busco, ".fna.tmp", sep = "")
    seq_file <- paste(out_dir, "/", busco, ".fna", sep = "")
    cat_commnad <- paste("cat $(cat ", path_file, ") > ", 
                         tmp_file, sep = "")
    rm_commnad1 <- paste("rm ", out_dir, "/", "*_na_paths", sep = "")
    rm_commnad2 <- paste("rm ", out_dir, "/", ".fna", sep = "")
    filtered_df %>%
      dplyr::pull(na_seq_path) %>%
      write(path_file)
    system(cat_commnad)
    phylotools::rename.fasta(infile = tmp_file, 
                             ref_table = ref_table, 
                             outfile = seq_file)
  } else {
    path_file <- paste(out_dir, "/", busco, "_aa_paths", sep = "")
    tmp_file <- paste(out_dir, "/", busco, ".faa.tmp", sep = "")
    seq_file <- paste(out_dir, "/", busco, ".faa", sep = "")
    cat_commnad <- paste("cat $(cat ", path_file, ") > ", 
                         tmp_file, sep = "")
    filtered_df %>%
      dplyr::pull(aa_seq_path) %>%
      write(path_file)
    system(cat_commnad)
    phylotools::rename.fasta(infile = tmp_file, 
                             ref_table = ref_table, 
                             outfile = seq_file)
  }
}

# Function to generate multifasta files for all busco loci across all taxa
# This is a wrapper for the full pipeline
# sample_ids      Text file with a list of the names of the genome samples. This should match the names of the sample-specific directories that contain the busco output
# busco_out_dir   Base path to directory with genome-specific BUSCO output
# busco_db        Name of the BUSCO database used in the runs. This is used to build the path to the results E.g., "run_cyanobacteria_odb1o"
# out_file_name   Name of the file with the raw busco output. Typically "full_table.tsv"
# data_type       "DNA" or "AA". Make sure that busco outputs the data type that you want to consolidate.
# out_dir         Path to directory to write the multifasta files
# busco_ids_file  Path to the file where you want the output file with the list of busco ids that were processed.
sort_busco_seqs <- function(sample_ids, busco_out_dir, busco_db, out_file_name,
                            data_type, out_dir, busco_ids_file) {
  rm_commnad1 <- paste("rm ", out_dir, "/", "*_paths", sep = "")
  rm_commnad2 <- paste("rm ", out_dir, "/", "*.tmp", sep = "")
  busco_df <- merge_busco(sample_ids = sample_ids, 
                          busco_out_dir = busco_out_dir, 
                          busco_db = busco_db, 
                          out_file_name = out_file_name) %>% 
    add_busco_paths(busco_out_dir = busco_out_dir, busco_db = busco_db) %>%
    add_busco_labels()
  busco_ids <- busco_df %>%
    dplyr::pull(busco_id) %>%
    unique()
  output <- lapply(busco_ids, cat_busco_seqs, 
                   busco_df_paths_labels = busco_df,
                   data_type = data_type,
                   out_dir = out_dir)
  if (!file.exists(busco_ids_file)) {
    write(busco_ids, file = busco_ids_file)
  }
  system(rm_commnad1)
  system(rm_commnad2)
}

# Function to summarize busco output. Note: the busco output table must be read
# with the funciton read_busc()
sum_busco <- function(busco_df) {
  busco_df %>% group_by(sample_id) %>%
    dplyr::distinct(busco_id, .keep_all = T) %>%
    dplyr::count(status) %>% 
    tidyr::pivot_wider(names_from = status, values_from = n) %>%
    dplyr::mutate(n_busco = 
             sum(Complete, Duplicated, Fragmented, Missing, na.rm = T)) %>%
    dplyr::mutate(percent_complete =
             (Complete/n_busco)*100) %>%
    dplyr::mutate(percent_duplicated =
             (Duplicated/n_busco)*100) %>%
    dplyr::mutate(percent_fragmented =
             (Fragmented/n_busco)*100) %>%
    dplyr::mutate(percent_missing =
             (Missing/n_busco)*100)
}

# Function read quast output from a single genome
# sample_id      The name of the assembly-specific directory with the quast output
# quast_out_dir  Base path to the directory with the quast output for all assemblies
read_quast <- function(sample_id, quast_out_dir) {
  quast_file <- paste(quast_out_dir, sample_id, "transposed_report.tsv", sep = "/")
  quast_out <- readr::read_delim(quast_file, delim = "\t")
  return(quast_out)
}

# Merge quast output tables
# sample_ids     Vector with the names of the assembly-specific directories with the quast output
# quast_out_dir  Base path to the directory with the quast output for all assemblies
merge_quast <- function(sample_ids, quast_out_dir) {
  quast_out_list <- lapply(sample_ids, read_quast, quast_out_dir = quast_out_dir)
  merged_quast_out <- dplyr::bind_rows(quast_out_list)
}

# Function to bin contigs from a dataframe
# sample_id        Id of the metagenomic sample from which the bins were generated
# bin_df           Data frame with contig labels (name), bin labels (bin), and sample ids (sample)
# contigs_file     Path to contigs to bin in FASTA format
# out_dir          Base path for the output directory
bin_from_df <- function(sample_id, bin_df, contigs_file, out_dir) {
  contigs <- ape::read.FASTA(contigs_file)
  bin_df <- bin_df %>%
    dplyr::filter(sample == sample_id)
  bin_n <- bin_df %>%
   dplyr::pull(bin) %>%
   unique()
 for (n in bin_n) {
   bin_file <- paste(out_dir, "/", sample_id, "/", sample_id, "_bin.", n, ".fa", sep = "")
   contig_labels <- bin_df %>%
     dplyr::filter(bin == n) %>%
     dplyr::pull(name)
   binned_contigs <- contigs[contig_labels]
   ape::write.FASTA(binned_contigs, file = bin_file)
 }
}

# Function to bin contigs from a dataframe for a single sample
# sample_id        Character vector with ids of the metagenomic samples from which the bins were generated
# bin_df           Data frame with contig labels (name), bin labels (bin), and sample ids (sample)
# contigs_file     Character vector with paths to contigs to bin in FASTA format
# out_dir          Base path for the output directory
v_bin_from_csv <- Vectorize(bin_from_df, 
                            vectorize.args = c("sample_id", "contigs_file"))

# Second version of function to bin contigs from a dataframe for a single sample.
# This one has more flexibility with output names and colum names of the input db
# sample_id        Id of the metagenomic sample from which the bins were generated
# bin_df           Data frame with contig labels (name), bin labels (bin), and sample ids (sample)
# contigs_file     Path to contigs to bin in FASTA format
# out_dir          Base path for the output directory
# out_suffix_1     Suffix to add to the output file names. By default, they will be "out_dir/sample_id/sample_id_binid.fa". If a suffix is provided, the output will be "out_dir/sample_idSUFFIX1/sample_id_binid.fa". The parent directory must exist.
# out_suffix_2     Suffix to add to the output file names. By default, they will be "out_dir/sample_id/sample_id_binid.fa". If a suffix is provided, the output will be "out_dir/sample_id/sample_id_SUFFIX2binid.fa". The parent directory must exist.
bin_from_df_single <- function(sample_id, bin_df, contigs_file, out_dir, 
                               out_suffix_1 = "", out_suffix_2 = "") {
  colnames(bin_df) <- c("name", "bin", "sample")
  bin_df <- bin_df %>%
    tidyr::replace_na(list(bin = "unassigned"))
  contigs <- ape::read.FASTA(contigs_file)
  bin_df <- bin_df %>%
    dplyr::filter(sample == sample_id)
  bin_n <- bin_df %>%
    dplyr::pull(bin) %>%
    unique()
  for (n in bin_n) {
    bin_file <- paste(out_dir, "/", sample_id, out_suffix_1, "/", 
                      sample_id, "_", out_suffix_2, n, ".fa", sep = "")
    contig_labels <- bin_df %>%
      dplyr::filter(bin == n) %>%
      dplyr::pull(name)
    binned_contigs <- contigs[contig_labels]
    ape::write.FASTA(binned_contigs, file = bin_file)
  }
}

# Vectorized function to bin contigs from a dataframe
# sample_id        Character vector with ids of the metagenomic samples from which the bins were generated
# bin_df           Data frame with contig labels (name), bin labels (bin), and sample ids (sample)
# contigs_file     Character vector with paths to contigs to bin in FASTA format
# out_dir          Base path for the output directory
# out_suffix       Suffix to add to the output file names. By default, they will be "out_dir/sample_id/sample_id_binid.fa". If a suffix is provided, the output will be "out_dir/sample_id/sample_id_SUFFIXbinid.fa"
bin_from_df_v2 <- Vectorize(bin_from_df_single, 
                            vectorize.args = c("sample_id", "contigs_file"))

# Function to replace and condense discovista conflict scoring keeping weak support category
condense_discov_out_weak <- function(discov_out) {
  # Keep only bipartition columns
  discov_new <- select(discov_out, 4:ncol(discov_out))
  # Replace conflict scoring
  discov_new <- apply(discov_new, MARGIN = 2, str_replace, pattern = "Missing", replacement = "uninformative")
  discov_new <- apply(discov_new, MARGIN = 2, str_replace, pattern = "Compatible \\(Weak Rejection\\)", replacement = "weak_reject")
  discov_new <- apply(discov_new, MARGIN = 2, str_replace, pattern = "Weak Support", replacement = "weak_support")
  discov_new <- apply(discov_new, MARGIN = 2, str_replace, pattern = "Strong Rejection", replacement = "discordant")
  discov_new <- apply(discov_new, MARGIN = 2, str_replace, pattern = "Strong Support", replacement = "concordant")
  # Make the df
  discov_new <- as.data.frame(discov_new) %>%
    pivot_longer(everything(), names_to = "bipart", values_to = "support") %>%
    count(bipart, support)  %>%
    pivot_wider(names_from = support, values_from = n) %>%
    mutate_if(is.integer, ~replace(., is.na(.), 0)) %>%
    mutate(no_loci =
             rowSums(across(concordant:weak_support))) %>%
    mutate(percent_concordant =
             (concordant/no_loci)*100) %>%
    mutate(percent_discordant = 
             (discordant/no_loci)*100) %>%
    mutate(percent_uninformative = 
             (uninformative/no_loci)*100) %>%
    mutate(percent_weak_reject =
             (weak_reject/no_loci)*100) %>%
    mutate(percent_weak_support =
             (weak_support/no_loci)*100)
  return(discov_new)
}

# Function to get a starter file with all bipartitions from a phylogenetic tree tree.
# The output can be used to generate the clade definition file required to run discovista
# tree    A rooted phylogenetic tree in newick format
get_biparts <- function(tree) {
  # Get the number of tips in the tree
  n_tips <- tree$tip.label %>% length()
  # Calculate the maximum number of internal nodes
  max_nodes <- n_tips - 1
  # Create the vector to store the bipartitions
  biparts <- character()
  # Iterate on all internal nodes of the tree
  for (i in 1:max_nodes) {
    node <- n_tips + i
    # Subset the tree to each internal node to get the taxa on one side of each bipartition
    tree_subset <- treeio::tree_subset(tree = tree, node = node, levels_back = 0)
    # Callapse taxa labels into a single character strin separated by commas
    biparts[i] <- tree_subset$tip.label %>%
      paste(collapse = ",")
  }
  return(biparts)
}


# Function to replace empty UFboot values for 0 in iqtree output treefiles
# treefiles   A character vector witht the path(s) to the trees. This function expect the treefile names to end in ".treefile" as this is the regular iqtree output
replace_empty_ufboot <- function(treefiles) {
  # Iterate over each treefile in the input vector
  for (treefile in treefiles) {
    # Load the tree and convert to tibble_tree
    tree <- ape::read.tree(file = treefile)
    # Replace empty UFBoot values with 0
    tree$node.label[tree$node.label == ""] <- "0"
    # Get the path for the new treefile
    new_treefile <- stringr::str_replace(treefile, ".treefile", "_noempty.treefile")
    # Write the edited tree
    write.tree(tree, file = new_treefile)
  }
}


# Function to calculate ANI GAP (Greatest Average nucleotide Identity gaP)
# fastani_df  Data frame where the first two columns have the genome names (must be named geonme1 and genome2), and the third column (named ani) has the ANI between genome1 and genome 2.
# low_lim     ANI lower limit to consider to calculate GAP
# hi_lim      ANI upper limit to consider to calculate GAP
gap <- function(fastani_df, low_lim, hi_lim) {
  genomes <- fastani_df %>% dplyr::pull(1) %>%
    base::unique()
  gap <- numeric()
  gap_low_lim <- numeric()
  gap_hi_lim <- numeric()
  for (i in 1:length(genomes)) {
    genome <- genomes[i]
    anis <- fastani_df %>% dplyr::filter(genome1 == genome | 
                                           genome2 == genome) %>% 
      dplyr::pull(ani) %>% 
      base::subset(. > low_lim & . <= hi_lim) %>%
      sort()
    diffs <- diff(anis)
    gap[i] <- max(diffs)
    if (is.finite(gap[i])) {
      gap_low_lim_index <- which(diffs == gap[i]) %>%
        base::unique()
      gap_low_lim[i] <- anis[gap_low_lim_index]
      gap_hi_lim[i] <- anis[gap_low_lim_index+1]
    } else {
      gap_low_lim[i] <- NA
      gap_hi_lim[i] <- NA
    }
  }
  out <- dplyr::bind_cols(genomes, gap, gap_low_lim, gap_hi_lim) %>%
    dplyr::mutate(`...2` =
                    dplyr::if_else(gap == -Inf, NA, gap))
  colnames(out) <- c("genomes", "gap",  "gap_low_lim", "gap_hi_lim")
  out
}


# Function to find synapomorphies for sets of taxa in a fasta sequence alignment
# alignment   DNAbin object with sequence alignment. Can obtain with ape::read.fasta()
# taxa set    Character vector with the names of the taxa in the alingment for which to find the synappomorphies
find_synapomorphies <- function(alignment, taxa_set) {
  # Convert alignment to data frame
  aln_df <- as.character(alignment) %>%
    as.data.frame()
  synapomorphies <- character()
  if (all(taxa_set %in% rownames(aln_df))) {
    other_taxa <- setdiff(rownames(aln_df), taxa_set)
  } else {
    stop("some of the taxa in the set are not part of the alignment")
  }
  for (site in names(aln_df)) {
    # Get the character states for the target taxa set
    unique_target_states <- aln_df[taxa_set, site] %>%
      unique()
    # Get the caracter states for the off-target set
    unique_other_states <- aln_df[other_taxa, site] %>%
      unique()
    # Check if there is only one unique value in the target set
    if (length(unique_target_states) == 1) {
      # Check to see if it is a synapomorphy
      if (!(unique_target_states %in% unique_other_states)) {
        synapomorphies[site] <- unique_target_states
      }
    }
  }
  # Check to see if any synapomorphies were found
  if (length(synapomorphies) > 0) {
    # Rename the apomorphic sites
    names(synapomorphies) <- stringr::str_replace(names(synapomorphies), 
                                                  "V", "site_")
  } else {
    synapomorphies <- c("none")
  }
  return(synapomorphies)
}

# Function to obtain a(x) from Degnan and Rosenberg 2016 Plos Biol.
deg_ros_a_of_x <- function(x) {
  num <- 3*(exp(1)^(2*x)) - 2
  den <- 18*((exp(1)^(3*x)) - (exp(1)^(2*x)))
  a_of_x <- log(2/3 + (num/den))
  return(a_of_x)
}


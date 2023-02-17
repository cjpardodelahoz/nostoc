#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse)
library(ape)
source("scripts/r_functions.R")

# Get the list of treefiles
treefiles <- list.files(path = "analyses/phylogenetics/set103/trees/single", 
                        pattern = ".*.treefile", 
                        full.names = T)
# Replace empty ufboot values for 0s
replace_empty_ufboot(treefiles)
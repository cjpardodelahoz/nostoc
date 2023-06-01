#!/usr/bin/env Rscript

# Load required packages and functions
library(tidyverse)

# Load voucher tables
# Note that some of the special characters in the voucher column may not be 
# parsed correctly
obrien_2013 <- read_delim(file = "misc_files/voucher_obrien2013.txt", 
                          skip = 1)
magain_2017 <- read_csv(file = "misc_files/voucher_magain2017.csv", 
                        locale=locale(encoding="latin1")) 
miadlikowska_2018 <- read_csv(file = "misc_files/voucher_miadlikowska2018.xlsx",
                              locale=locale(encoding="latin1")) # HAVE TO REMOVE SUPERSCRIPTS FROM SOURCE BEFORE BRINGING IN
chagnon_2018 <- read_csv(file = "misc_files/voucher_chagnon_2018.csv", 
                         locale=locale(encoding="latin1"))
pardodelahoz_2018 <- read_csv(file = "misc_files/voucher_pardodelahoz2018.csv",
                              locale=locale(encoding="latin1")) # RBCLX ACCESSIONS HAVE ASTERISKS
magain_2018 <- 
                         
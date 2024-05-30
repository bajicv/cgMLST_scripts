################################################################################
# Create a hashed profile for provided hash-cgMLST output dir
#
# Author: Vladimir Bajić
# Date: 2024-05-30
#
# Description:
#   This script takes as input directory with json files produced by hash-cgMLST
#   It combines all of them into single tsv profile
#   where empty hashes are replaced with 0
#
# Usage:
#
# To see help message
#   Rscript --vanilla jsons2tsv_profile4hash-cgMLST.R --help
#
# To make hashed profile from json files
#   Rscript --vanilla jsons2tsv_profile4hash-cgMLST.R -i /path_to/hashcgmlst_out_dir/ -o hashed_profile.tsv
#
################################################################################

# Load libraries ---------------------------------------------------------------
library(optparse)
suppressMessages(library(tidyverse))
suppressMessages(library(jsonlite))


# Making option list -----------------------------------------------------------
option_list <- list(
    make_option(c("-i", "--input_dir"), type = "character", help = "Path to the directory with json files outputed by hash-cgMLST", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", help = "Output name", metavar = "character")
)

# Parsing options --------------------------------------------------------------
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check the provided option and execute the corresponding code -----------------
if (is.null(opt$i)) {
    print_help(opt_parser)
    stop("Provide path to the directory with json files outputed by hash-cgMLST.")
}

if (is.null(opt$o)) {
    print_help(opt_parser)
    stop("Provide output name.")
}

# Find all json files in directory ---------------------------------------------
files <- dir(opt$i, pattern = "/*.json", full.names = TRUE)

# Read all json files in one table ---------------------------------------------
hashed_profile <-
    files %>%
    map_df(~ fromJSON(.)) %>%
    select(-file, -scheme) %>%
    as_tibble() %>%
    select(name, alleles) %>%
    mutate(
        hash = unlist(alleles),
        locus = names(alleles)
    ) %>%
    select(-alleles) %>%
    rename(`#Sample` = "name") %>%
    pivot_wider(names_from = "locus", values_from = "hash")

# Replace empty hashes with 0s -------------------------------------------------
hashed_profile[hashed_profile == ""] <- "0"

# Save hash-cgMLST_profile tsv table -------------------------------------------
write_tsv(hashed_profile, opt$o)

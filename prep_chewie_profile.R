################################################################################
# Prepare ChewBBACA's profile outputted by AlleleCall
#
# Author: Vladimir Bajić
# Date: 2024-04-04
#
# Description:
# This script
#   - Removes "INF-" from allele names
#   - Replaces all other ChewBBACA's non-numeric allele call categories with 0s
#
# Usage:
#
# To see help message
#   Rscript --vanilla prep_chewie_profile.R --help
#
# To prepare chewBBACA's profile and save prepared profile in same dir where is input
#   Rscript --vanilla prep_chewie_profile.R -i results_alleles.tsv
#
# To prepare chewBBACA's profile and save prepared profile with specific name/path
#   Rscript --vanilla prep_chewie_profile.R -i results_alleles.tsv -o my_dir/results_alleles_prepared.tsv
#
################################################################################

# Libraries --------------------------------------------------------------------
suppressMessages(library(tidyverse))
library(optparse)

# Making option list -----------------------------------------------------------
option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "Path to the chewBBACA cgMLST profile before AlleleCallEvaluator", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", help = "Output name prefix", metavar = "character")
)

# Parsing options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check the provided option and execute the corresponding code -----------------
if (is.null(opt$i)) {
    print_help(opt_parser)
    stop("Provide path to the chewBBACA profile file.")
}

# Check if output is specified and if not use input to define it
if (is.null(opt$o)) {
    opt$o <- tools::file_path_sans_ext(opt$i)
    cat("Output not specified.\nOutput will be saved as: ", opt$o, "\n", sep = "")
}

# Load data --------------------------------------------------------------------
cat("Loading profile ...\n")
profile <- read_tsv(opt$i, show_col_types = FALSE, col_types = cols(.default = "c"))

# Set problematic allele calls to 0s -------------------------------------------
profile_prepared <-
    profile %>%
    mutate(across(.cols = -1, ~ str_replace(., "^[[:space:]]+", ""))) %>%
    mutate(across(.cols = -1, ~ str_replace(., "^INF-", ""))) %>%
    mutate(across(.cols = -1, ~ str_replace(., "[[:alpha:]]+", "0")))

# Save output ------------------------------------------------------------------
cat("Saving profile as: ", paste0(opt$o, "_prepared.tsv", "\n"))
write_tsv(profile_prepared, paste0(opt$o, "_prepared.tsv"))

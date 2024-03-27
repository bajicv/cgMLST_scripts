################################################################################
# Compare two cgMLST profiles for same set of samples and loci
#
# Author: Vladimir Bajić
# Date: March 2024
#
# Description:
# This script
#   - Compares two cgMLST profiles for the same set of samples and loci
#   - Outputs comparison summary per locus (suffix: _LCOMP.tsv)
#   - Outputs comparison summary per sample (suffix: _SCOMP.tsv)
#
# Categories used in summary:
#   S   - same allele call in both profiles
#   D   - different allele calls (excluding differences caused by 0s)
#   M   - missing allele calls in both profiles
#   M1  - missing allele call in the first profile but not in the second profile
#   M2  - missing allele call in the second profile but not in the first profile
#
# Usage:
#
# To see help message
#   Rscript --vanilla compare_two_cgMLST_profiles.R --help
#
# To compare two profiles
#   Rscript --vanilla compare_two_cgMLST_profiles.R --f first.profile --s second.profile -o output_path
#
################################################################################

# Libraries  -------------------------------------------------------------------
suppressMessages(library(tidyverse))
library(optparse)

# To suppress summarise information
options(dplyr.summarise.inform = FALSE)

# Making option list -----------------------------------------------------------
option_list <- list(
    make_option(c("-f", "--first_profile"), type = "character", help = "Path to the first cgMLST profile", metavar = "character"),
    make_option(c("-s", "--second_profile"), type = "character", help = "Path to the second cgMLST profile", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", help = "Output name prefix", metavar = "character")
)

# Parsing options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check the provided option and execute the corresponding code -----------------
if (is.null(opt$f)) {
    print_help(opt_parser)
    stop("Provide path to the first profile file.")
}

if (is.null(opt$s)) {
    print_help(opt_parser)
    stop("Provide path to the second profile file.")
}

if (is.null(opt$o)) {
    print_help(opt_parser)
    stop("Output path must be provided.")
}

# ------------------------------------------------------------------------------


# Load data
cat("Loading profile 1 ...\n")
profile_1 <- read_tsv(opt$f, show_col_types = FALSE) %>% rename_with(.cols = 1, ~"Sample")

cat("Loading profile 2 ...\n")
profile_2 <- read_tsv(opt$s, show_col_types = FALSE) %>% rename_with(.cols = 1, ~"Sample")

# Test if samples are the same
if (sum(!c(profile_1$Sample %in% profile_2$Sample)) + sum(!c(profile_2$Sample %in% profile_1$Sample)) == 0) {
    cat("Samples in both profiles are identical.\n")
} else {
    cat("WARNING: Samples in both profiles are NOT identical.")
    cat("Number of samples in profile 1: ", nrow(profile_1), "\n")
    cat("Number of samples in profile 2: ", nrow(profile_2), "\n")
}

# Test if loci are the same
if (sum(!c(names(profile_1) %in% names(profile_2))) + sum(!c(names(profile_2) %in% names(profile_1))) == 0) {
    cat("Loci in both profiles are identical.\n")
} else {
    cat("WARNING: Loci in both profiles are NOT identical.\n")
    cat("Number of loci in profile 1: ", ncol(profile_1), "\n")
    cat("Number of loci in profile 2: ", ncol(profile_2), "\n")
}

# Convert to long formats ------------------------------------------------------
lp1 <- pivot_longer(profile_1, cols = -1, names_to = "Locus", values_to = "Allele_1")
lp2 <- pivot_longer(profile_2, cols = -1, names_to = "Locus", values_to = "Allele_2")


# Join two profiles ------------------------------------------------------------
comp <- full_join(lp1, lp2, by = c("Sample", "Locus"))

# Make classifications ---------------------------------------------------------
classified_comp <- comp %>%
    mutate(
        Classification = case_when(
            Allele_1 == Allele_2 & Allele_1 != 0 & Allele_2 != 0 ~ "S",
            Allele_1 != Allele_2 & Allele_1 != 0 & Allele_2 != 0 ~ "D",
            Allele_1 == Allele_2 & Allele_1 == 0 & Allele_2 == 0 ~ "M",
            Allele_1 != Allele_2 & Allele_1 == 0 & Allele_2 != 0 ~ "M1",
            Allele_1 != Allele_2 & Allele_1 != 0 & Allele_2 == 0 ~ "M2",
            TRUE ~ "Unknown"
        )
    )

# List of expected Classification column names ---------------------------------
classifications <- c("S", "D", "M", "M1", "M2") %>% map_dfr(~ tibble(!!.x := logical()))

# Comparison per Sample --------------------------------------------------------
s_comp <-
    classified_comp %>%
    group_by(Sample, Classification) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = "Classification", values_from = "n", values_fill = 0) %>%
    bind_rows(., classifications) %>%
    replace(is.na(.), 0) %>%
    select("Sample", "S", "D", "M", "M1", "M2")

# Comparison per Locus ---------------------------------------------------------
l_comp <-
    classified_comp %>%
    group_by(Locus, Classification) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = "Classification", values_from = "n", values_fill = 0) %>%
    bind_rows(., classifications) %>%
    replace(is.na(.), 0) %>%
    select("Locus", "S", "D", "M", "M1", "M2")

# Save output ------------------------------------------------------------------
cat("Saving the tables as: ", paste0(opt$o, "_profile_{LCOMP,SCOMP}.tsv", "\n"))
write_tsv(l_comp, paste0(opt$o, "_profile_LCOMP.tsv"))
write_tsv(s_comp, paste0(opt$o, "_profile_SCOMP.tsv"))

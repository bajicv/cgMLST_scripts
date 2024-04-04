################################################################################
# Compare two cgMLST profiles for same set of samples and loci
#
# Author: Vladimir Bajić
# Date: 2024-04-04
#
# Description:
# This script
#   - Compares two cgMLST profiles for the same set of samples and loci
#   - Outputs general comparison summary (suffix: _LCOMP.tsv)
#   - Outputs comparison summary per locus (suffix: _LCOMP.tsv)
#   - Outputs comparison summary per sample (suffix: _SCOMP.tsv)
#
# Categories used in summary:
#   S   - same allele call in both profiles
#   D   - different allele calls (excluding differences caused by 0s)
#   M   - missing allele calls in both profiles
#   M1  - missing allele call in the first profile but not in the second profile
#   M2  - missing allele call in the second profile but not in the first profile
#   X1  - missing loci in the first profile
#   X2  - missing loci in the second profile
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

# To suppress summarize information
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

# DEV TEST ONLY
opt$f <- "/scratch/Projekte/MF1_GE/Research_Projects/cgMLST_benchmarking/Salmonella/1_analysis/5_seqsphere/Salmonella_Testdaten_FG11.profile_NoLastTab_0s.tsv"
opt$s <- "/scratch/projekte/MF1_GE/Research_Projects/cgMLST_benchmarking/Salmonella/1_analysis/2_hash_cgMLST/Salmonella_Testdaten_FG11_hash-cgMLST.profile"

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


# Load data --------------------------------------------------------------------
cat("Loading profile 1 ...\n")
profile_1 <- read_tsv(opt$f, show_col_types = FALSE) %>% rename_with(.cols = 1, ~"Sample")

cat("Loading profile 2 ...\n")
profile_2 <- read_tsv(opt$s, show_col_types = FALSE) %>% rename_with(.cols = 1, ~"Sample")

# Calculate basic stats --------------------------------------------------------

## Number of samples in each profile
nsp1 <- nrow(profile_1)
nsp2 <- nrow(profile_2)

## Number of loci in each profile
nlp1 <- ncol(profile_1)
nlp2 <- ncol(profile_2)

# Test if samples are the same -------------------------------------------------
if (sum(!c(profile_1$Sample %in% profile_2$Sample)) + sum(!c(profile_2$Sample %in% profile_1$Sample)) == 0) {
    cat("Same sample IDs present in both profiles.\n")
    if (nsp1 != nsp2) {
        cat("WARNING: Different number of samples in profiles! \n")
    }
    if (nsp1 == nsp2) {
        cat("Same number of samples in both profiles.\n")
    }
} else {
    cat("WARNING: Samples in both profiles are NOT identical!")
}

# Print number of samples in each profile --------------------------------------
cat("Number of samples in profile 1: ", nsp1, "\n")
cat("Number of samples in profile 2: ", nsp2, "\n")

# Test if samples are duplicated and list those that are -----------------------
if (sum(duplicated(profile_1$Sample)) > 0) {
    cat("Duplicated samples detected in profile 1: ", profile_1$Sample[duplicated(profile_1$Sample)], "\n")
}
if (sum(duplicated(profile_2$Sample)) > 0) {
    cat("Duplicated samples detected in profile 2: ", profile_2$Sample[duplicated(profile_2$Sample)], "\n")
}

# Test if loci are the same ----------------------------------------------------

lp1_notin_lp2 <- sum(!c(names(profile_1) %in% names(profile_2)))
lp2_notin_lp1 <- sum(!c(names(profile_2) %in% names(profile_1)))

if (lp1_notin_lp2 + lp2_notin_lp1 == 0) {
    cat("Loci in both profiles are identical.\n")
} else {
    cat("WARNING: Loci in both profiles are NOT identical!\n")
    if (lp1_notin_lp2 > 0) {
        cat("Number of loci from profile 1 not present in profile 2:", lp1_notin_lp2, "\n")
    }
    if (lp2_notin_lp1 > 0) {
        cat("Number of loci from profile 2 not present in profile 1:", lp2_notin_lp1, "\n")
    }
}

# Print number of loci in each profile -----------------------------------------
cat("Number of loci in profile 1: ", nlp1, "\n")
cat("Number of loci in profile 2: ", nlp2, "\n")

# Test if loci are duplicated and list those that are --------------------------
if (sum(duplicated(names(profile_1))) > 0) {
    cat("Duplicated loci detected in profile 1: ", names(profile_1)[duplicated(names(profile_1))], "\n")
}
if (sum(duplicated(names(profile_2))) > 0) {
    cat("Duplicated loci detected in profile 2: ", names(profile_2)[duplicated(names(profile_2))], "\n")
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
            is.na(Allele_1) & !is.na(Allele_2) ~ "X1",
            !is.na(Allele_1) & is.na(Allele_2) ~ "X2",
            TRUE ~ "Unknown"
        )
    )

# List of expected Classification column names ---------------------------------
classifications <- c("S", "D", "M", "M1", "M2", "X1", "X2") %>% map_dfr(~ tibble(!!.x := logical()))


# General comparison -----------------------------------------------------------
# if SUM is different from Total_COMP something is wrong!
g_comp <-
    classified_comp %>%
    group_by(Classification) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = "Classification", values_from = "n", values_fill = 0) %>%
    bind_rows(., classifications) %>%
    replace(is.na(.), 0) %>%
    mutate(
        Total_COMP = nrow(classified_comp),
        SUM = S + D + M + M1 + M2 + X1 + X2,
        pct_S = S / Total_COMP * 100,
        pct_D = D / Total_COMP * 100,
        pct_M = M / Total_COMP * 100,
        pct_M1 = M1 / Total_COMP * 100,
        pct_M2 = M2 / Total_COMP * 100,
        pct_X1 = X1 / Total_COMP * 100,
        pct_X2 = X2 / Total_COMP * 100,
        Profile_1 = opt$f,
        Profile_2 = opt$s,
        NLP1 = nlp1,
        NLP2 = nlp2,
        NSP1 = nsp1,
        NSP2 = nsp2,
        LP1_notin_LP2 = lp1_notin_lp2,
        LP2_notin_LP1 = lp2_notin_lp1
    ) %>%
    relocate(Profile_1, Profile_2, NLP1, NLP2, NSP1, NSP2, LP1_notin_LP2, LP2_notin_LP1, Total_COMP, SUM, S, D, M, M1, M2, X1, X2, pct_S, pct_D, pct_M, pct_M1, pct_M2, pct_X1, pct_X2)

# Comparison per Sample --------------------------------------------------------
s_comp <-
    classified_comp %>%
    group_by(Sample, Classification) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = "Classification", values_from = "n", values_fill = 0) %>%
    bind_rows(., classifications) %>%
    replace(is.na(.), 0) %>%
    select("Sample", "S", "D", "M", "M1", "M2", "X1", "X2") %>%
    arrange(S, desc(D), desc(M), desc(M1), desc(M2), desc(X1), desc(X2))

# Comparison per Locus ---------------------------------------------------------
l_comp <-
    classified_comp %>%
    group_by(Locus, Classification) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = "Classification", values_from = "n", values_fill = 0) %>%
    bind_rows(., classifications) %>%
    replace(is.na(.), 0) %>%
    select("Locus", "S", "D", "M", "M1", "M2", "X1", "X2") %>%
    arrange(S, desc(D), desc(M), desc(M1), desc(M2), desc(X1), desc(X2))

# Save output ------------------------------------------------------------------
cat("Saving the tables as: ", paste0(opt$o, "_profile_{GCOMP,LCOMP,SCOMP}.tsv", "\n"))
write_tsv(g_comp, paste0(opt$o, "_profile_GCOMP.tsv"))
write_tsv(l_comp, paste0(opt$o, "_profile_LCOMP.tsv"))
write_tsv(s_comp, paste0(opt$o, "_profile_SCOMP.tsv"))

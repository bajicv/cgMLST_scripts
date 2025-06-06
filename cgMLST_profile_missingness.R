################################################################################
# Calculate and plot cgMLST profile missingness
#
# Author: Vladimir Bajić
# Date: June 2025
#
# Description:
# This script
#   - Calculates missingness per sample and per locus based on cgMLST profile
#   - Plots histograms for missingness per locus (LMISS) and missingness per sample (SMISS)
#     (file suffix: _profile_miss.png)
#   - Outputs two new profiles QC-ed for missing loci.
#     1. Profile without missing alleles
#        (file suffix: _QC_LFMISS_0.profile)
#     2. Profile with loci missing alleles in less than 10% of samples
#        (file suffix: _QC_LFMISS_01.profile)
#   - Outputs two tables with LMISS (.lmiss) and SMISS (.smiss)
#
# Usage:
#
# To see help message
#   Rscript --vanilla cgMLST_profile_missingness.R --help
#
# To plot and calculate missingness and save files using same prefix as input
#   Rscript --vanilla cgMLST_profile_missingness.R -i my_profile.profile
#
# To plot and calculate missingness and save files using specific output
#   Rscript --vanilla cgMLST_profile_missingness.R -i my_profile.profile -o my_out
#
################################################################################

# Libraries  -------------------------------------------------------------------
suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
library(optparse)


# Making option list -----------------------------------------------------------
option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "Path to the input profile file [tsv]", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", help = "Output name prefix", metavar = "character")
)
# Parsing options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# Check the provided option and execute the corresponding code -----------------

if (is.null(opt$i)) {
    print_help(opt_parser)
    stop("Input file must be provided.\n")
}

# Check if output is specified and if not use input to define it
if (is.null(opt$o)) {
    opt$o <- tools::file_path_sans_ext(opt$i)
    cat("Output not specified.\nOutput will be saved as: ", opt$o, sep = "")
}

# Load data --------------------------------------------------------------------
profile <- read_tsv(opt$i, col_types = cols(.default = col_character()), show_col_types = FALSE) %>% rename_with(.cols = 1, ~"SAMPLE")

# Make long format of profile --------------------------------------------------
profile_long <-
    profile %>%
    pivot_longer(-1, names_to = "LOCUS", values_to = "ALLELE")

# Calculate missingness --------------------------------------------------------

# Extract number of samples and loci -------------------------------------------
n_samples <- nrow(profile)
n_loci <- ncol(profile)

## LMISS: missingness by locus. Number of samples missing per locus.
#
#  Column explanations:
#  LOCUS       Locus ID
#  N_MISS      Number of missing samples per locus
#  F_MISS      Frequency of missing samples per locus

LMISS <-
    profile_long %>%
    filter(ALLELE == 0) %>%
    group_by(LOCUS) %>%
    summarise(N_MISS = n(), F_MISS = N_MISS / n_samples) %>%
    mutate(QC = case_when(F_MISS >= 0.1 ~ "FAIL", F_MISS < 0.1 ~ "PASS"))

## SMISS: missingness by sample. Number of loci missing per sample.
#
#  Column explanations:
#  SAMPLE      Sample ID
#  N_MISS      Number of missing loci per sample
#  F_MISS      Frequency of missing loci per sample

SMISS <-
    profile_long %>%
    filter(ALLELE == 0) %>%
    group_by(SAMPLE) %>%
    summarise(N_MISS = n(), F_MISS = N_MISS / n_loci) %>%
    mutate(QC = case_when(F_MISS >= 0.1 ~ "FAIL", F_MISS < 0.1 ~ "PASS"))


## Find loci missing in more than 10% of samples
LFMISS_01 <-
    LMISS %>%
    filter(F_MISS >= 0.1) %>%
    arrange(desc(F_MISS)) %>%
    pull(LOCUS)


## Find samples missing more than 10% of loci
SFMISS_01 <-
    SMISS %>%
    filter(F_MISS >= 0.1) %>%
    arrange(desc(F_MISS)) %>%
    pull(SAMPLE)

# Save missingness tables ------------------------------------------------------
write_tsv(SMISS, paste0(opt$o, ".smiss"))
write_tsv(LMISS, paste0(opt$o, ".lmiss"))

# Plotting histograms ----------------------------------------------------------

## LMISS plot
p_lmiss <-
    LMISS %>%
    ggplot(aes(x = N_MISS, fill = QC)) +
    geom_histogram(binwidth = 1) +
    xlab("Number of missing samples") +
    ylab("Locus count") +
    ggtitle(
        label = "LMISS histogram",
        subtitle = paste0("No. Samples = ", n_samples, "\n10% of Samples  = ", n_samples / 10, "\nNo. loci missing in ≥10% of samples = ", length(LFMISS_01))
    ) +
    geom_vline(xintercept = n_samples / 10, color = "#f46d43", linetype = "dashed") +
    scale_fill_manual(name = "QC", values = c("PASS" = "#66bd63", "FAIL" = "#f46d43")) +
    theme_bw() +
    theme(legend.position = "none")


## SMISS plot
p_smiss <-
    SMISS %>%
    ggplot(aes(x = N_MISS, fill = QC)) +
    geom_histogram() +
    xlab("Number of missing loci") +
    ylab("Sample count") +
    ggtitle(
        label = "SMISS histogram",
        subtitle = paste0("No. Loci = ", n_loci, "\n10% of Loci  = ", n_loci / 10, "\nNo. samples missing ≥10% of loci = ", length(SFMISS_01))
    ) +
    geom_vline(xintercept = n_loci / 10, color = "#f46d43", linetype = "dashed") +
    scale_fill_manual(name = "QC", values = c("PASS" = "#66bd63", "FAIL" = "#f46d43")) +
    theme_bw() +
    theme(legend.position = "none")

## Arrange plots in one
p_miss <- ggarrange(p_smiss, p_lmiss)

## Annotate plot with main title indicating file on which it was performed
p_miss_ann <- annotate_figure(p_miss,
    top = text_grob(paste0("Missingness: ", tools::file_path_sans_ext(basename(opt$i)), "\n"),
        color = "Black",
        face = "bold",
        size = 16
    )
)

## Save plot
ggsave(paste0(opt$o, "_profile_miss.png"), p_miss_ann, width = 10, height = 5, bg = "white")



# 100% and 90% QC for loci -----------------------------------------------------

## Find loci missing in at least one sample
LFMISS_0 <-
    LMISS %>%
    filter(F_MISS > 0) %>%
    arrange(desc(F_MISS)) %>%
    pull(LOCUS)

## Find loci missing in more than 10% of samples
LFMISS_01 <-
    LMISS %>%
    filter(F_MISS >= 0.1) %>%
    arrange(desc(F_MISS)) %>%
    pull(LOCUS)

## Make a new conservative profiles where each locus is present in all samples
profile_QC_LFMISS_0 <-
    profile %>%
    select(names(profile)[!names(profile) %in% LFMISS_0])

## Make new profile QC-ed for loci missing in more than 10% of samples
profile_QC_LFMISS_01 <-
    profile %>%
    select(names(profile)[!names(profile) %in% LFMISS_01])


## Write QC-ed profiles --------------------------------------------------------
write_tsv(profile_QC_LFMISS_0, paste0(opt$o, "_QC_LFMISS_0.profile"))
write_tsv(profile_QC_LFMISS_01, paste0(opt$o, "_QC_LFMISS_01.profile"))

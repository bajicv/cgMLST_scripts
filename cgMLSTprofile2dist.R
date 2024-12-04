################################################################################
# Calculate and plot distance matrix from cgMLST profile
#
# Author: Vladimir Bajić
# Date: 2024-12-04
#
# Description:
# This script
#  - calculate the hamming distance matrix (with zeroes excluded) from
#    tsv formatted cgMLST allele call table (e.g. ChewBBACA or hash-cgmlst).
#  - accepts profiles with allele names that are either numeric, characters, or hashes.
#  - can produce heat plot and histogram depicting pairwise allele distances between samples [--plot].
#
# Note:
# If providing chewBBACA profile, you should converts/masks all non-integer
# classifications in the profile to 0 and remove all the INF- prefixes
# BEFORE using this script.
#
# Usage:
#
# To see help message
#   Rscript --vanilla cgMLSTprofile2dist.R --help
#
# To calculate distance matrix
#   Rscript --vanilla cgMLSTprofile2dist.R -i input.profile -o output_path
#
# To calculate distance matrix and plot histogram and heat map with pairwise allele distances
#   Rscript --vanilla cgMLSTprofile2dist.R -i input.profile -o output_path --plot
#
################################################################################

# Libraries  -------------------------------------------------------------------
suppressMessages(library(tidyverse))
library(optparse)
library(tidyheatmaps)

# To suppress summarize information
options(dplyr.summarise.inform = FALSE)

# Making option list -----------------------------------------------------------
option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "Path to the input cgMLST profile [tsv]", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", help = "Output name prefix", metavar = "character"),
    make_option(c("-p", "--plot"), action = "store_true", help = "Output plotts. [default: %default]", default = FALSE)
)

# Parsing options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the input is provided -----------------------------------------------
if (is.null(opt$i)) {
    print_help(opt_parser)
    stop("Path to the input cgMLST profile.")
}

# Check if output is specified and if not use input to define it
if (is.null(opt$o)) {
    opt$o <- tools::file_path_sans_ext(opt$i)
    cat("Output not specified.\nOutput will be saved as: ", opt$o, sep = "")
}

# ------------------------------------------------------------------------------


# Load data --------------------------------------------------------------------
cat("Loading profile ...\n")
profile <- read_tsv(opt$i, show_col_types = FALSE, col_types = cols(.default = "c")) %>% rename_with(.cols = 1, ~"Sample")

# Calculate basic stats --------------------------------------------------------

## Number of samples in profile
nsp <- nrow(profile)

## Number of loci in profile
nlp <- ncol(profile)

# Print number of samples in each profile --------------------------------------
cat("Number of samples in profile: ", nsp, "\n")

# Test if samples are duplicated and list those that are -----------------------
if (sum(duplicated(profile$Sample)) > 0) {
    cat("Duplicated samples detected: ", profile$Sample[duplicated(profile$Sample)], "\n")
}

# Print number of loci in each profile -----------------------------------------
cat("Number of loci in profile: ", nlp, "\n")

# Test if loci are duplicated and list those that are --------------------------
if (sum(duplicated(names(profile))) > 0) {
    cat("Duplicated loci detected: ", names(profile)[duplicated(names(profile))], "\n")
}

# Convert to long formats ------------------------------------------------------
lp <- pivot_longer(profile, cols = -1, names_to = "Locus", values_to = "Allele_1")

# Functions --------------------------------------------------------------------

## Function to calculate pairwise distances between samples in a given cgMLST profile
profile_pwdist_df <- function(profile) {
    pivot_longer(profile, cols = -1, names_to = "Locus", values_to = "Allele") %>%
        left_join(., ., by = "Locus", relationship = "many-to-many") %>%
        select(Locus, Sample.x, Sample.y, Allele.x, Allele.y) %>%
        filter(Allele.x != 0 & Allele.y != 0) %>%
        mutate(
            Classification = case_when(
                Allele.x == Allele.y ~ "0",
                Allele.x != Allele.y ~ "1"
            )
        ) %>%
        group_by(Sample.x, Sample.y) %>%
        summarise(dist = sum(Classification == 1)) %>%
        pivot_wider(., names_from = "Sample.y", values_from = "dist") %>%
        rename(Sample = "Sample.x")
}

## Function to calculate pairwise distances between samples in a given cgMLST profile
profile_pwdist_long_tibble <- function(profile) {
    pivot_longer(profile, cols = -1, names_to = "Locus", values_to = "Allele") %>%
        left_join(., ., by = "Locus", relationship = "many-to-many") %>%
        select(Locus, Sample.x, Sample.y, Allele.x, Allele.y) %>%
        filter(Allele.x != 0 & Allele.y != 0) %>%
        mutate(
            Classification = case_when(
                Allele.x == Allele.y ~ "0",
                Allele.x != Allele.y ~ "1"
            )
        ) %>%
        group_by(Sample.x, Sample.y) %>%
        summarise(dist = sum(Classification == 1))
}

## Calculate allelic pairwise distances between samples ------------------------
cat("Calculating pairwise allele distances between samples in the profile. \n")

## In data frame format (to safe as tables)
pwdist <- profile_pwdist_df(profile)

## In long format (to use for plotting)
pwdist_long <- profile_pwdist_long_tibble(profile)

## Save tables with allelic pairwise distances between samples -----------------
cat("Saving matrix with pairwise allele distances between samples in profile. \n")
write_tsv(pwdist, paste0(opt$o, "_PWAD.tsv"))

################################################################################
# Check if flag for plots is specified and if so produce plots

if (opt$p) {
    # Heat plots of pairwise distances -----------------------------------------
    cat("Plotting heat plot with pairwise allele distances.\n")

    ## Define colors for heat plots
    my_palette_dist <- colorRampPalette(c("white", "#40004b"))(6000)

    ## Heat plot of pairwise differences between samples in P1
    rg <- max(abs(pwdist_long$dist))

    p_pwdist <-
        tidy_heatmap(
            pwdist_long,
            rows = Sample.x,
            columns = Sample.y,
            values = dist,
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            display_numbers = FALSE,
            clustering_distance_rows = "euclidean",
            clustering_distance_cols = "euclidean",
            clustering_method = "complete",
            main = "Pairwise allele distances between samples",
            colors = my_palette_dist,
            color_legend_n = max(pwdist_long$dist),
            silent = TRUE
        )

    # Save heat plot -----------------------------------------------------------
    cat("Saving heat plot.\n")
    ggsave(paste0(opt$o, "_heatplot_PWAD.png"), p_pwdist, width = 10, height = 10)

    # Histograms of pairwise allele distances ----------------------------------
    cat("Plotting histogram of pairwise allele distances.\n")

    ## Histogram of pairwise distances between samples in Profile
    phist <-
        pwdist_long %>%
        ggplot(aes(x = dist)) +
        geom_histogram(binwidth = 1, fill = "black") +
        theme_bw() +
        ggtitle("Histogram of pairwise allele distances between samples")

    # Save histogram ---------------------------------------------------------------
    cat("Saving histogram.\n")
    ggsave(paste0(opt$o, "_histograms.png"), phist, width = 10, height = 5)
}
cat("Done.\n")
################################################################################

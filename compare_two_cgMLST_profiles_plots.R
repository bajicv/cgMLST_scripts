################################################################################
# Plot results of compare_two_cgMLST_profiles.R
#
# Author: Vladimir Bajić
# Date: 2024=04-03
#
# Description:
# This script
#   - Outputs a plot summarizing comparison per locus (suffix: _LCOMP.png)
#   - Outputs a plot summarizing comparison per sample (suffix: _SCOMP.png)
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
#   Rscript --vanilla compare_two_cgMLST_profiles_plot.R --help
#
# To run the script it is enough to provide path to either LCOMP or SCOMP (assuming that LCOMP and SCOMP are in the same directory)
#   Rscript --vanilla compare_two_cgMLST_profiles_plot.R -l LCOMP.tsv
#   Rscript --vanilla compare_two_cgMLST_profiles_plot.R -s SCOMP.tsv
#
# To run the script with specific output path, threshold and input paths for LCOMP and SCOMP
#   Rscript --vanilla compare_two_cgMLST_profiles_plot.R -l LCOMP.tsv -s SCOMP.tsv - o out_path -t 95
#
################################################################################

# Libraries  -------------------------------------------------------------------
suppressMessages(library(tidyverse))
suppressMessages(library(ggpubr))
library(optparse)

# Making option list -----------------------------------------------------------
option_list <- list(
    make_option(c("-l", "--lcomp"), type = "character", help = "Path to LCOMP file [tsv].", metavar = "character"),
    make_option(c("-s", "--scomp"), type = "character", help = "Path to SCOMP file [tsv].", metavar = "character"),
    make_option(c("-t", "--threshold"), type = "numeric", help = "Threshould value [1-100]. Plot only loci with the same allele calls in less than provided threshould value. Default value = 90", metavar = "numeric", default = 90),
    make_option(c("-o", "--output"), type = "character", help = "Output name prefix.", metavar = "character")
)

# Parsing options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# Check the provided option and execute the corresponding code -----------------
if (is.null(opt$l) & is.null(opt$s)) {
    print_help(opt_parser)
    stop("Provide path to LCOMP or SCOMP file.")
}

# Check if path to SCOMP is specified and if not use LCOMP to define it
if (is.null(opt$s)) {
    cat("SCOMP not specified.\nSCOMP will be searched in the same folder where is LCOMP.\n")
    opt$s <- str_replace(opt$l, "_LCOMP.tsv$", "_SCOMP.tsv")
}

# Check if path to LCOMP is specified and if not use SCOMP to define it
if (is.null(opt$l)) {
    cat("LCOMP not specified.\nSCOMP will be searched in the same folder where is SCOMP.\n")
    opt$l <- str_replace(opt$s, "_SCOMP.tsv$", "_LCOMP.tsv")
}

# Check if output is specified and if not use input to define it
if (is.null(opt$o)) {
    opt$o <- dirname(opt$l)
    cat("Output not specified.\nOutput will be saved in: ", opt$o, "\n", sep = "")
}

# ------------------------------------------------------------------------------

# Load data --------------------------------------------------------------------
cat("Reading in the data.\n")
LCOMP <- read_tsv(opt$l, show_col_types = FALSE)
SCOMP <- read_tsv(opt$s, show_col_types = FALSE)


# Plot SAMPLE ------------------------------------------------------------------

my_colors <-
    c(
        "M" = "#1f78b4",
        "M1" = "#ffff99",
        "M2" = "#fdbf6f",
        "S" = "#33a02c",
        "D" = "#e31a1c",
        "X1" = "#969696",
        "X2" = "#252525"
    )

cat("Plotting SCOMP.\n")
p_s <- SCOMP %>%
    arrange(desc(S), D, M, M1, M2, X1, X2) %>%
    pivot_longer(cols = -1, names_to = "Type", values_to = "Count") %>%
    mutate(Sample = forcats::fct_inorder(Sample)) %>%
    ggplot(aes(x = Count, y = Sample, fill = Type)) +
    geom_bar(stat = "identity") +
    labs(
        title = "Comparisons per SAMPLE",
        x = "Loci",
        y = "Sample"
    ) +
    scale_fill_manual(values = my_colors) +
    theme_bw()

# Plot LOCUS -------------------------------------------------------------------
cat("Plotting LCOMP.\n")
p_l <- LCOMP %>%
    arrange(desc(S), D, M, M1, M2, X1, X2) %>%
    filter(S < nrow(SCOMP) * (opt$t) / 100) %>%
    pivot_longer(cols = -1, names_to = "Type", values_to = "Count") %>%
    mutate(Locus = forcats::fct_inorder(Locus)) %>%
    ggplot(aes(x = Count, y = Locus, fill = Type)) +
    geom_bar(stat = "identity") +
    labs(
        title = paste0("Comparisons per LOCUS"),
        subtitle = paste0("Shown are only loci with the same allele calls in less than ", opt$t, "% of samples"),
        x = "Samples",
        y = "Locus",
    ) +
    scale_fill_manual(values = my_colors) +
    theme_bw()

# Save plots -------------------------------------------------------------------
cat("Saving plots.\n")
ggsave(paste0(opt$o, "/", tools::file_path_sans_ext(basename(opt$s)), ".png"), p_s, width = 10, height = 20)
ggsave(paste0(opt$o, "/", tools::file_path_sans_ext(basename(opt$l)), ".png"), p_l, width = 10, height = 20)

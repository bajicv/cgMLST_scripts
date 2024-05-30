################################################################################
# Compare two cgMLST profiles for same set of samples and loci
#
# Author: Vladimir Bajić
# Date: 2024-04-09
#
# Description:
# This script
#   - Compares two cgMLST profiles for the same set of samples and loci
#   - Outputs a table with general comparison summary (suffix: _GCOMP.tsv)
#   - Outputs a table with comparison summary per locus (suffix: _LCOMP.tsv)
#   - Outputs a table with comparison summary per sample (suffix: _SCOMP.tsv)
#   - Outputs a long format table with comparison summary per sample including clarifications (suffix: _CCOMP.tsv)
#   - Outputs a table with allelic pairwise distances between samples in profile 1 (suffix: _PWD_P1.tsv)
#   - Outputs a table with allelic pairwise distances between samples in profile 2 (suffix: _PWD_P2.tsv)
#   - Outputs a table with differences between allelic pairwise distances between samples in two profiles (suffix: _DIFF_PWD_P1_P2.tsv)
#   - Outputs a plot summarizing general comparison between two profiles (suffix: _GCOMP.png)
#   - Outputs a plot summarizing comparison per locus (suffix: _LCOMP.png)
#   - Outputs a plot summarizing comparison per sample (suffix: _SCOMP.png)
#   - Outputs histograms with pairwise distances between samples in both profiles, and differences between PWDs between the two profiles (suffix: _histograms.png)
#   - Outputs heatplots with pairwise distances between samples in both profiles (suffix: _PWD_P1.tsv; _PWD_P2.tsv), and differences between PWDs between the two profiles (suffix: _DIFF_PWD_P1_P2.tsv)
#
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
# To compare two profiles with specific output and threshold value for LCOMP plot
#   Rscript --vanilla compare_two_cgMLST_profiles.R --f first.profile --s second.profile -o output_path -t 95
#
################################################################################

# Libraries  -------------------------------------------------------------------
suppressMessages(library(tidyverse))
library(optparse)
library(tidyheatmaps)
library(ggpubr)

# To suppress summarize information
options(dplyr.summarise.inform = FALSE)

# Making option list -----------------------------------------------------------
option_list <- list(
    make_option(c("-f", "--first_profile"), type = "character", help = "Path to the first cgMLST profile", metavar = "character"),
    make_option(c("-s", "--second_profile"), type = "character", help = "Path to the second cgMLST profile", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", help = "Output name prefix", metavar = "character"),
    make_option(c("-t", "--threshold"), type = "numeric", help = "Threshould value [1-100]. Plot only loci with the same allele calls in less than provided threshould value. Default value = 90", metavar = "numeric", default = 90)
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


# Load data --------------------------------------------------------------------
cat("Loading profile 1 ...\n")
profile_1 <- read_tsv(opt$f, show_col_types = FALSE, col_types = cols(.default = "c")) %>% rename_with(.cols = 1, ~"Sample")

cat("Loading profile 2 ...\n")
profile_2 <- read_tsv(opt$s, show_col_types = FALSE, col_types = cols(.default = "c")) %>% rename_with(.cols = 1, ~"Sample")

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
# if SUM is different from TotaLCOMP something is wrong!
GCOMP <-
    classified_comp %>%
    group_by(Classification) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = "Classification", values_from = "n", values_fill = 0) %>%
    bind_rows(., classifications) %>%
    replace(is.na(.), 0) %>%
    mutate(
        TotaLCOMP = nrow(classified_comp),
        SUM = S + D + M + M1 + M2 + X1 + X2,
        pct_S = S / TotaLCOMP * 100,
        pct_D = D / TotaLCOMP * 100,
        pct_M = M / TotaLCOMP * 100,
        pct_M1 = M1 / TotaLCOMP * 100,
        pct_M2 = M2 / TotaLCOMP * 100,
        pct_X1 = X1 / TotaLCOMP * 100,
        pct_X2 = X2 / TotaLCOMP * 100,
        Profile_1 = opt$f,
        Profile_2 = opt$s,
        NLP1 = nlp1,
        NLP2 = nlp2,
        NSP1 = nsp1,
        NSP2 = nsp2,
        LP1_notin_LP2 = lp1_notin_lp2,
        LP2_notin_LP1 = lp2_notin_lp1
    ) %>%
    relocate(Profile_1, Profile_2, NLP1, NLP2, NSP1, NSP2, LP1_notin_LP2, LP2_notin_LP1, TotaLCOMP, SUM, S, D, M, M1, M2, X1, X2, pct_S, pct_D, pct_M, pct_M1, pct_M2, pct_X1, pct_X2)

# Comparison per Sample --------------------------------------------------------
SCOMP <-
    classified_comp %>%
    group_by(Sample, Classification) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = "Classification", values_from = "n", values_fill = 0) %>%
    bind_rows(., classifications) %>%
    replace(is.na(.), 0) %>%
    select("Sample", "S", "D", "M", "M1", "M2", "X1", "X2") %>%
    arrange(S, desc(D), desc(M), desc(M1), desc(M2), desc(X1), desc(X2)) %>%
    ungroup()

# Comparison per Locus ---------------------------------------------------------
LCOMP <-
    classified_comp %>%
    group_by(Locus, Classification) %>%
    summarise(n = n()) %>%
    pivot_wider(names_from = "Classification", values_from = "n", values_fill = 0) %>%
    bind_rows(., classifications) %>%
    replace(is.na(.), 0) %>%
    select("Locus", "S", "D", "M", "M1", "M2", "X1", "X2") %>%
    arrange(S, desc(D), desc(M), desc(M1), desc(M2), desc(X1), desc(X2)) %>%
    ungroup()

# Save output ------------------------------------------------------------------
cat("Saving the tables as: ", paste0(opt$o, "_profile_{GCOMP,LCOMP,SCOMP,raw_CCOMP}.tsv", "\n"))
write_tsv(classified_comp, paste0(opt$o, "_raw_CCOMP.tsv"))
write_tsv(GCOMP, paste0(opt$o, "_profile_GCOMP.tsv"))
write_tsv(LCOMP, paste0(opt$o, "_profile_LCOMP.tsv"))
write_tsv(SCOMP, paste0(opt$o, "_profile_SCOMP.tsv"))



# Define colors for the plot ---------------------------------------------------
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

# Plot GENERAL -----------------------------------------------------------------
cat("Plotting GCOMP\n")
p_g <- GCOMP %>%
    select(S, D, M, M1, M2, X1, X2) %>%
    pivot_longer(cols = everything(), names_to = "Type", values_to = "Count") %>%
    mutate(COMP = "P1_vs_P2") %>%
    ggplot(aes(x = Count, y = COMP, fill = Type)) +
    geom_bar(stat = "identity") +
    labs(
        title = "General Comparisons",
        subtitle = paste0(
            "\nProfile_1: ", GCOMP$Profile_1, "\n",
            "Profile_2: ", GCOMP$Profile_2, "\n\n",
            "NLP1: ", GCOMP$NLP1, "\n",
            "NLP2: ", GCOMP$NLP2, "\n\n",
            "NSP1: ", GCOMP$NSP1, "\n",
            "NSP2: ", GCOMP$NSP2, "\n\n",
            "LP1_notin_LP2: ", GCOMP$LP1_notin_LP2, "\n",
            "LP2_notin_LP1: ", GCOMP$LP2_notin_LP1, "\n\n",
            "TotaLCOMP: ", GCOMP$TotaLCOMP, "\n\n"
        ),
        x = "Total number of comparisons",
        y = ""
    ) +
    scale_fill_manual(values = my_colors) +
    theme_bw()

# Plot SAMPLE ------------------------------------------------------------------
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
    mutate(Locus = forcats::fct_inorder(Locus)) %>%
    filter(S < nrow(SCOMP) * (opt$t) / 100) %>%
    pivot_longer(cols = -1, names_to = "Type", values_to = "Count") %>%
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
cat("Saving GCOMP, SCOMP, and LCOMP plots.\n")
ggsave(paste0(opt$o, "_profile_GCOMP.png"), p_g, width = 20, height = 10)
ggsave(paste0(opt$o, "_profile_SCOMP.png"), p_s, width = 10, height = 20)
ggsave(paste0(opt$o, "_profile_LCOMP.png"), p_l, width = 10, height = 20)


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
cat("Calculating allelic pairwise distances between samples in profiles. \n")

## In data frame format (to safe as tables)
pwdist_p1 <- profile_pwdist_df(profile_1)
pwdist_p2 <- profile_pwdist_df(profile_2)

## In long format (to use for plotting)
pwdist_p1_long <- profile_pwdist_long_tibble(profile_1)
pwdist_p2_long <- profile_pwdist_long_tibble(profile_2)

## Calculate difference between P1 and P2
cat("Calculating differences between two profiles' allelic pairwise distances between samples. \n")
diff_p1_p2 <-
    full_join(pwdist_p1_long, pwdist_p2_long, by = join_by(Sample.x, Sample.y)) %>%
    mutate(diff = dist.x - dist.y) %>%
    select(-dist.x, -dist.y)

## Save tables with allelic pairwise distances between samples -----------------
cat("Saving tables with allelic pairwise distances between samples in profiles. \n")
write_tsv(pwdist_p1, paste0(opt$o, "_PWD_P1.tsv"))
write_tsv(pwdist_p2, paste0(opt$o, "_PWD_P2.tsv"))

cat("Saving table with differences between allelic pairwise distances between samples in two profiles. \n")
diff_p1_p2 %>%
    pivot_wider(names_from = Sample.y, values_from = diff) %>%
    rename(Sample = "Sample.x") %>%
    write_tsv(., paste0(opt$o, "_DIFF_PWD_P1_P2.tsv"))

# Heat plots of pairwise distances ---------------------------------------------
cat("Plotting heat plots with pairwise distances.\n")

## Define colors for heat plots
my_palette_dist <- colorRampPalette(c("white", "#40004b"))(6000)
my_palette_diff <- colorRampPalette(c("#1b7837", "white", "#40004b"))(6000)

## Heat plot of pairwise differences between samples in P1
rg <- max(abs(pwdist_p1_long$dist))
p_pwdist_p1 <- tidy_heatmap(pwdist_p1_long,
    rows = Sample.x,
    columns = Sample.y,
    values = dist,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    main = "Pairwise distances between samples in P1",
    colors = my_palette_dist,
    color_legend_n = max(pwdist_p1_long$dist),
    silent = TRUE
)

## Heat plot of pairwise differences between samples in P2
p_pwdist_p2 <- tidy_heatmap(pwdist_p2_long,
    rows = Sample.x,
    columns = Sample.y,
    values = dist,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    main = "Pairwise distances between samples in P2",
    colors = my_palette_dist,
    color_legend_n = max(pwdist_p2_long$dist),
    silent = TRUE
)

## Heat plot of differences between P1 and P2
rg <- max(abs(diff_p1_p2$diff))
p_diff_p1_p2 <-
    diff_p1_p2 %>%
    tidy_heatmap(
        .,
        rows = Sample.x,
        columns = Sample.y,
        values = diff,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        display_numbers = FALSE,
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        clustering_method = "complete",
        main = "P1-P2",
        colors = my_palette_diff,
        color_legend_min = -rg,
        color_legend_max = rg,
        color_legend_n = rg * 2,
        silent = TRUE
    )

# Save heat plots --------------------------------------------------------------
cat("Saving heat plots.\n")
ggsave(paste0(opt$o, "_heatplot_PWD_P1.png"), p_pwdist_p1, width = 10, height = 10)
ggsave(paste0(opt$o, "_heatplot_PWD_P2.png"), p_pwdist_p2, width = 10, height = 10)
ggsave(paste0(opt$o, "_heatplot_DIFF_PWD_P1_P2.png"), p_diff_p1_p2, width = 10, height = 10)

# Histograms of pairwise distances ---------------------------------------------
cat("Plotting histograms of pairwise distances.\n")

## Histogram of pairwise distances between samples in Profile 1
phist_p1 <-
    pwdist_p1_long %>%
    ggplot(aes(x = dist)) +
    geom_histogram(binwidth = 1, fill = "black") +
    theme_bw() +
    ggtitle("PWD between samples in Profile 1")

## Histogram of pairwise distances between samples in Profile 2
phist_p2 <-
    pwdist_p2_long %>%
    ggplot(aes(x = dist)) +
    geom_histogram(binwidth = 1, fill = "black") +
    theme_bw() +
    ggtitle("PWD between samples in Profile 2")

## Histogram of differences between pairwise distances between samples in Profile 1 and Profile 2
phist_p1_p2_diff <-
    diff_p1_p2 %>%
    ggplot(aes(x = diff)) +
    geom_histogram(binwidth = 1, fill = "darkred") +
    theme_bw() +
    ggtitle("Difference in PWDs between P1 and P2") +
    geom_vline(xintercept = 0, color = "gray", linetype = "dashed")

## Arrange all histograms in one plot
p_histograms <-
    ggarrange(
        ggarrange(phist_p1, phist_p2, nrow = 2, align = "h"),
        phist_p1_p2_diff,
        nrow = 1, ncol = 2,
        common.legend = TRUE,
        legend = "bottom"
    )

# Save histograms --------------------------------------------------------------
cat("Saving histograms.\n")
ggsave(paste0(opt$o, "_histograms.png"), p_histograms, width = 10, height = 5)

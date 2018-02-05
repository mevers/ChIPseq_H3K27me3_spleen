#!/usr/bin/env Rscript

# The R command line script estimates copy numbers of
# ribosomal and mitochondrial DNA from a comparison of read coverages.
#
# Author: Maurits Evers (maurits.evers@anu.edu.au)
# Original date: 04/07/2017
# Last change: 28/07/2017

# Todo:
#  1. Add principle component analysis
#  2. Show coverage per chromosome

## ------------------------------------------------------------------------
## Start clock
t0 <- Sys.time();


## ------------------------------------------------------------------------
## Load libraries
suppressMessages(library("optparse"));    # Python-style command line args
suppressMessages(library("ggplot2"));     # Plotting
suppressMessages(library("ggrepel"));     # Repelling labels
suppressMessages(library("hrbrthemes"));  # ggplot2 colour scheme


## ------------------------------------------------------------------------
## Parse command line arguments
option_list <- list(
    make_option(
        c("-i", "--input"),
        type = "character",
        default = NULL,
        help = "Tab-separated input file of read counts (e.g. from deeptools)",
        metavar = "character"),
    make_option(
        c("-o", "--outdir"),
        type = "character",
        default = ".",
        help = "Output directory [default %default]",
        metavar = "character"),
    make_option(
        c("--scale"),
        action = "store_true",
        type = "logical",
        default = TRUE,
        help = "Scale entries to have unit variance [default %default]",
        metavar = "character"),
    make_option(
        c("--center"),
        action = "store_true",
        type = "logical",
        default = TRUE,
        help = "Shift entries to be zero centered [default %default]",
        metavar = "character"),
    make_option(
        c("--maxN"),
        type = "integer",
        default = 0,
        help = "Number of entries used for PCA [default %default]",
        metavar = "integer")
);
opt_parser <- OptionParser(option_list = option_list);
args <- parse_args(opt_parser);
if (is.null(args$input)) {
    print_help(opt_parser);
    stop("Must give input file.\n", call. = FALSE);
}


## ------------------------------------------------------------------------
## Timestamp function
ts <- function() {
    return(format(Sys.time(), "[%a %b %d %Y %H:%M:%S]"));
}


## ------------------------------------------------------------------------
## Custom function
ggpca <- function(df, x, y, var.pca) {
    gg <- ggplot(
        df,
        aes_string(x = x, y = y, label = "sample", colour = "condition"));
    gg <- gg + geom_point(size = 3, alpha = 0.6);
    gg <- gg + stat_ellipse(type = "t", linetype = 2);
    gg <- gg + geom_text_repel(show.legend = FALSE);
    gg <- gg + theme_bw();
    gg <- gg + scale_colour_ipsum(name = "Condition");
    gg <- gg + labs(
        x = sprintf("%s (%3.2f %%)", x, var.pca[x] * 100),
        y = sprintf("%s (%3.2f %%)", y, var.pca[y] * 100),
        title = sprintf(
            "PCA based on %s",
            basename(input)),
        subtitle = ifelse(
            isTRUE(scale),
            "Unit variance scaling per sample",
            "No unit variance scaling"),
        caption = sprintf(
            "[Maurits Evers (maurits.evers@anu.edu.au), %s]",
            format(Sys.Date(), format = "%B %Y")));
    gg <- gg + theme(
      rect = element_rect(fill = "transparent"),
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", colour = NA),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 16),
      strip.text.y = element_text(angle = 45),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      plot.title = element_text(size = 16),
      legend.box = "vertical",
      legend.key = element_blank(),
      legend.position = "bottom");
    lapply(c("pdf", "png"), function(z) ggsave(
        file = sprintf(
            "%s/plot.PCA.%s_vs_%s.%s.%s",
            outdir,
            x,
            y,
            tools::file_path_sans_ext(basename(input)),
            z),
        gg,
        bg = "transparent",
        width = 8,
        height = 6));

}


## ------------------------------------------------------------------------
## Global variables
input <- args$input;
outdir <- gsub("/$", "", args$outdir);
scale <- args$scale;
center <- args$center;
maxN <- args$maxN;
cat(sprintf("%s Parameter summary\n", ts()));
cat(sprintf(" input          = %s\n", input));
cat(sprintf(" outdir         = %s\n", outdir));
cat(sprintf(
    " scale          = %s\n",
    ifelse(isTRUE(scale), "TRUE", "FALSE")));
cat(sprintf(
    " center         = %s\n",
    ifelse(isTRUE(center), "TRUE", "FALSE")));
cat(sprintf(" maxN           = %s\n", maxN));


## ------------------------------------------------------------------------
# Check if input files and output directory exists
if (!file.exists(input)) {
    stop(
        sprintf("Input file %s does not exists.\n", input),
        call. = FALSE);
}
if (!file.exists(outdir)) {
    stop(
        sprintf("Output directory %s does not exists.\n", outdir),
        call. = FALSE);
}


## ------------------------------------------------------------------------
# Read in data
cat(sprintf(
    "%s Reading data.\n",
    ts()));
data <- read.delim(input, header = TRUE);
sel <- seq(4, ncol(data));
data <- data[order(apply(data[, sel], 1, var), decreasing = TRUE), ];
if (maxN == 0 | maxN > nrow(data)) {
    maxN <- nrow(data);
}
data.maxN <- data[1:maxN, ];


## ------------------------------------------------------------------------
# Perform PCA
cat(sprintf(
    "%s Performing PCA.\n",
    ts()));
df <- data.maxN[apply(data.maxN[, sel], 1, var) > 0, sel];
if (isTRUE(scale)) {
    df <- scale(df);
}
pca <- prcomp(
    t(df),
    center = center,
    scale. = FALSE);
var.pca <- pca$sdev^2;
var.pca <- var.pca / sum(var.pca);
names(var.pca) <- colnames(pca$x);


## ------------------------------------------------------------------------
# Prepare for plotting
df <- as.data.frame(pca$x);
df$sample <- gsub("(^X.|.$)", "", rownames(df));
df$condition <- ifelse(
    grepl("(input|control)", df$sample),
    "input",
    "ChIP");


## ------------------------------------------------------------------------
# Plot PC1 vs. PC2
cat(sprintf(
    "%s Producing plot.\n",
    ts()));
ggpca(df, "PC1", "PC2", var.pca);
ggpca(df, "PC2", "PC3", var.pca);


## ------------------------------------------------------------------------
# Done
cat(sprintf(
    "%s All done.\n",
    ts()));

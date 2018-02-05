#!/usr/bin/env Rscript

# The R command line script plots the coverage distribution per
# reference sequence
#
# Author: Maurits Evers (maurits.evers@anu.edu.au)
# Original date: 10/08/2017
# Last change: 10/08/2017

## ------------------------------------------------------------------------
## Start clock
t0 <- Sys.time();


## ------------------------------------------------------------------------
## Load libraries
suppressMessages(library(optparse));    # Python-style command line args
suppressMessages(library(ggplot2));     # Plotting
suppressMessages(library(reshape2));    # Reshaping dataframes
suppressMessages(library(hrbrthemes));  # ggplot2 color scheme


## ------------------------------------------------------------------------
## Parse command line arguments
option_list <- list(
    make_option(
        c("-i", "--input"),
        type = "character",
        default = NULL,
        help = "(Comma-separated list of) input file(s)",
        metavar = "character"),
    make_option(
        c("-o", "--outdir"),
        type = "character",
        default = ".",
        help = "Output directory [default %default]",
        metavar = "character"),
    make_option(
        c("--skipZeros"),
        action = "store_true",
        type = "logical",
        default = FALSE,
        help = "Skip regions with zero coverage [default %default]"),
    make_option(
        c("--normaliseToAutosomeCov"),
        action = "store_true",
        type = "logical",
        default = FALSE,
        help = "Normalise coverage to average autosomal genome coverae [default %default]")

);
opt_parser <- OptionParser(option_list = option_list);
args <- parse_args(opt_parser);
if (is.null(args$input)) {
    print_help(opt_parser);
    stop("At least one coverage file must be supplied.\n", call. = FALSE);
}


## ------------------------------------------------------------------------
## Timestamp function
ts <- function() {
    return(format(Sys.time(), "[%a %b %d %Y %H:%M:%S]"));
}


## ------------------------------------------------------------------------
## Custom function
calculateCoverage <- function(
    fn,
    outdir = ".") {
    df <- read.table(fn, stringsAsFactors = FALSE);
    # Normalise to length of window to coverage per bp
    cov <- df[, 5] / (df[, 3] - df[, 2]);
    names(cov) <- factor(df$V1, levels = unique(df$V1));
    return(cov);
}


## ------------------------------------------------------------------------
## Global variables
fn <- unlist(strsplit(args$input, ","));
outdir <- gsub("/$", "", args$outdir);
skipZeros <- args$skipZeros;
normaliseToAutosomeCov <- args$normaliseToAutosomeCov;
cat(sprintf("%s Parameter summary\n", ts()));
for (i in 1:length(fn)) {
    cat(sprintf(" Input                  = %s\n", fn[i]));
}
cat(sprintf(" outdir                 = %s\n", outdir));
cat(sprintf(
    " skipZeros              = %s\n",
    ifelse(isTRUE(skipZeros), "TRUE", "FALSE")));
cat(sprintf(
    " normaliseToAutosomeCov = %s\n",
    ifelse(isTRUE(normaliseToAutosomeCov), "TRUE", "FALSE")));


## ------------------------------------------------------------------------
# Check if input files and output directory exists
for (i in 1:length(fn)) {
    if (!file.exists(fn[i])) {
        stop(
            sprintf("Input file %s does not exists.\n", fn[i]),
            call. = FALSE);
    }
}
if (!file.exists(outdir)) {
    stop(
        sprintf("Output directory %s does not exists.\n", outdir),
        call. = FALSE);
}


## ------------------------------------------------------------------------
# Determine coverage distribution and estimate CN
cat(sprintf(
    "%s Extracting coverage profile per reference sequence per sample\n",
    ts()));
flagZeros <- ifelse(isTRUE(skipZeros), "noZeros", "withZeros");
cov <- sapply(fn, function(x) calculateCoverage(x, outdir = outdir));
colnames(cov) <- gsub(
    "(cov.genome_w\\d+_s\\d+.|.bed)", "",
    basename(fn));


## ------------------------------------------------------------------------
# Remove entries that are all zero
if (isTRUE(skipZeros)) {
    cov <- cov[which(rowSums(cov) > 0), ];
}


## ------------------------------------------------------------------------
# Remove entries that are all zero
if (isTRUE(normaliseToAutosomeCov)) {
    meanCovAutosome <- colMeans(cov[which(rownames(cov) %in% seq(1, 22)), ]);
    medianCovAutosome <- apply(
        cov[which(rownames(cov) %in% seq(1, 22)), ],
        2,
        median);
    d <- diag(1.0/medianCovAutosome);
    colnames(d) <- names(meanCovAutosome);
    cov <- cov %*% d;
}


## ------------------------------------------------------------------------
# Prepare for plotting
cat(sprintf(
    "%s Producing plot.\n",
    ts()));
df <- melt(cov);
df$condition <- ifelse(
    grepl("(input|control)", df$Var2),
    "input",
    "ChIP");
df$Var2 <- factor(df$Var2, levels = sort(as.character(unique(df$Var2))));


## ------------------------------------------------------------------------
# Plot
gg <- ggplot(df, aes(x = Var1, y = value, fill = condition));
gg <- gg + geom_boxplot(outlier.shape = NA, width = 0.7);
gg <- gg + facet_wrap(~ Var2, nrow = length(unique(df$condition)));
gg <- gg + theme_bw();
gg <- gg + scale_y_log10();
gg <- gg + geom_hline(yintercept = 1, col = "red", linetype = 3);
gg <- gg + labs(
    x = "Reference",
    y = ifelse(
        isTRUE(normaliseToAutosomeCov),
        "log10(Coverage normalised to median autosomal coverage)",
        "log10(Coverage per bp)"),
    title = ifelse(
        isTRUE(normaliseToAutosomeCov),
        "Normalised coverage per reference sequence per sample",
        "Coverage per bp per reference sequence per sample"),
    subtitle = sprintf(
        "Note: Coverage (%s) is shown on a log10 scale",
        ifelse(
            isTRUE(normaliseToAutosomeCov),
            "normalised to median autosomal genome coverage",
            "in bp")),
    caption = sprintf(
        "[Maurits Evers (maurits.evers@anu.edu.au), %s]",
        format(Sys.Date(), format = "%B %Y")));
gg <- gg + scale_fill_ipsum(name = "Condition");
gg <- gg + theme(
  rect = element_rect(fill = "transparent"),
  panel.background = element_rect(fill = "transparent"),
  plot.background = element_rect(fill = "transparent", colour = NA),
  strip.background = element_blank(),
  strip.text.x = element_text(size = 16),
  strip.text.y = element_text(angle = 45),
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  plot.title = element_text(size = 16),
  legend.box = "vertical",
  legend.key = element_blank(),
  legend.position = "bottom");
lapply(c("pdf", "png"), function(x) ggsave(
    file = sprintf(
        "%s/plot.%s.%s",
        outdir,
        "coverage",
        x),
    gg,
    bg = "transparent",
    width = 18,
    height = 12));


## ------------------------------------------------------------------------
# Table
cat(sprintf(
    "%s Producing table.\n",
    ts()));
tab <- do.call(
    data.frame,
    aggregate(value ~ Var1 + Var2, function(x)
        c(median = median(x), mad = mad(x)), data = df));
tab <- dcast(data = tab, Var1 ~ Var2, value.var = "value.median");
rownames(tab) <- tab[, 1];
tab <- tab[, -1];
write.csv(
    tab,
    file = sprintf("%s/coverage_per_bp.csv", outdir),
    quote = FALSE);
#enrichment <- tab[nrow(tab), ] / colMeans(
#    tab[which(rownames(tab) %in% c(seq(1, 22), "X", "Y")), ]);


## ------------------------------------------------------------------------
# Done
cat(sprintf(
    "%s All done.\n",
    ts()));

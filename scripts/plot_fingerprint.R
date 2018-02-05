#!/usr/bin/env Rscript

# The R command line script shows a fingerprint plot
# based on genome read coverage table generated from
# deeptools plotFingerprint
#
# Author: Maurits Evers (maurits.evers@anu.edu.au)
# Original date: 03/08/2017
# Last change: 03/08/2017


## ------------------------------------------------------------------------
## Start clock
t0 <- Sys.time();


## ------------------------------------------------------------------------
## Load libraries
suppressMessages(library("optparse"));    # Python-style command line args
suppressMessages(library("reshape2"));    # Data wrangling
suppressMessages(library("ggplot2"));     # Plotting
suppressMessages(library("hrbrthemes"));  # ggplot2 colour scheme


## ------------------------------------------------------------------------
## Parse command line arguments
option_list <- list(
    make_option(
        c("-i", "--input"),
        type = "character",
        default = NULL,
        help = "Tab-separated input file of read counts
                (from deeptools plotFingerprint)",
        metavar = "character"),
    make_option(
        c("-o", "--outdir"),
        type = "character",
        default = ".",
        help = "Output directory [default %default]",
        metavar = "character")
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
## Global variables
input <- args$input;
outdir <- gsub("/$", "", args$outdir);
cat(sprintf("%s Parameter summary\n", ts()));
cat(sprintf(" input          = %s\n", input));
cat(sprintf(" outdir         = %s\n", outdir));


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
colnames(data) <- gsub("(X.|.$)", "", colnames(data));


## ------------------------------------------------------------------------
# Split data based on samples
samples <- colnames(data);
l <- lapply(samples, function(x) data[, colnames(data) == x]);
names(l) <- samples;


## ------------------------------------------------------------------------
# Sort read coverage per sample by increasing order
cat(sprintf(
    "%s Calculating normalised rank and cumulative sum based on coverage.\n",
    ts()));
l <- lapply(l, function(x) x[order(x, decreasing = FALSE)]);


## ------------------------------------------------------------------------
# Calculate normalised rank and cumulative sum
l <- lapply(l, function(x) cbind.data.frame(
    normRank = rank(x, ties.method = "first") / length(x),
    normCumsum = cumsum(x) / sum(x)
));


## ------------------------------------------------------------------------
# Prepare for plotting
df <- melt(l, id.vars = colnames(l[[1]]));
df$L1 <- factor(df$L1, levels = sort(samples));
df$condition <- ifelse(
    grepl("(input|control)", df$L1),
    "input",
    "ChIP");
#df2 <- lapply(split(df, df$condition), function(x) {
#    ret <- aggregate(
#        normCumsum ~ .,
#        function(y) c(mean = mean(y), sd = sd(y)),
#        data = subset(x, select = -c(L1)));
#    do.call(data.frame, ret);
#});
#df2 <- melt(df2, id.vars = colnames(df2[[1]]));
#df2 <- subset(df2, select = -c(L1));
#df2$l <- ifelse(
#    is.na(df2[, "normCumsum.sd"]),
#    df2[, "normCumsum.mean"],
#    df2[, "normCumsum.mean"] - 1.96 * df2[, "normCumsum.sd"]);
#df2$h <- ifelse(
#    is.na(df2[, "normCumsum.sd"]),
#    df2[, "normCumsum.mean"],
#    df2[, "normCumsum.mean"] + 1.96 * df2[, "normCumsum.sd"]);

## ------------------------------------------------------------------------
# Plot
cat(sprintf(
    "%s Producing plot.\n",
    ts()));
gg <- ggplot(df, aes(x = normRank, y = normCumsum));
gg <- gg + geom_abline(intercept = 0, slope = 1, linetype = 3);
gg <- gg + geom_line(aes(colour = condition));
#gg <- gg + geom_line(aes(colour = condition), size = 1.5);
#gg <- gg + geom_line(
#    data = df2,
#    aes(x = normRank, y = normCumsum.mean, colour = condition),
#    alpha = 0.8);
#gg <- gg + geom_ribbon(
#    data = df2,
#    aes(x = normRank, ymin = l, ymax = h),
#    alpha = 0.6);
gg <- gg + facet_wrap(~ L1, nrow = length(unique(df$condition)));
gg <- gg + theme_bw();
gg <- gg + labs(
    x = "Fraction of genome",
    y = "Fraction of total reads",
    title = "Fingerprints of different samples",
    subtitle = sprintf(
        "Normalised rank of genome coverage vs. normalised cumulative sum of read coverage"),
    caption = sprintf(
        "[Maurits Evers (maurits.evers@anu.edu.au), %s]",
        format(Sys.Date(), format = "%B %Y")));
gg <- gg + scale_colour_ipsum(name = "Condition");
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
lapply(c("pdf", "png"), function(x) ggsave(
    file = sprintf(
        "%s/plot.%s.%s",
        outdir,
        tools::file_path_sans_ext(basename(input)),
        x),
    gg,
    bg = "transparent",
    width = 12,
    height = 8));


## ------------------------------------------------------------------------
# Done
cat(sprintf(
    "%s All done.\n",
    ts()));

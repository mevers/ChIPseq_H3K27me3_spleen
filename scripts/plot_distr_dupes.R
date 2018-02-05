#!/usr/bin/env Rscript

# The R command line script shows a metagene profile
# based on a score matrix from deeptools computeMatrix
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
suppressMessages(library("Rsamtools"));   # Read BAM files
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
        help = "Inut BAM file; can be comma-separated list of
                multiple BAM files",
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
fn <- unlist(strsplit(args$input, ","));
outdir <- gsub("/$", "", args$outdir);
cat(sprintf("%s Parameter summary\n", ts()));
for (i in 1:length(fn)) {
    cat(sprintf(" Input          = %s\n", fn[i]));
}
cat(sprintf(" outdir         = %s\n", outdir));


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
# Read in BAM files
cat(sprintf(
    "%s Reading BAM files.\n",
    ts()));
bam <- lapply(fn, BamFile);
stopifnot(length(unique(lapply(bam, seqlengths))) == 1);
refsize <- seqlengths(bam[[1]]);
reads <- lapply(bam, scanBam);


## ------------------------------------------------------------------------
# Get reference sequence information from reads and tabulate
refID <- lapply(reads, function(x) x[[1]]$rname);
refID <- lapply(refID, function(x) {
    levels(x) <- names(refsize);
    return(x);
})
refTab <- lapply(refID, function(x) as.data.frame(table(x)));


## ------------------------------------------------------------------------
# Normalise to total number reads
# Normalise to size of reference sequence
cat(sprintf(
    "%s Tabulate data.\n",
    ts()));
refTab <- lapply(refTab, function(x) {
    x$fracTotal <- x$Freq / sum(x$Freq);
    x$refsize <- refsize;
    x$DPK <- x$Freq / refsize * 1000;
    x;
})
names(refTab) <- gsub("(.bam|.dupes)", "", basename(fn));


## ------------------------------------------------------------------------
# Reformat data
df <- melt(refTab, id.vars = colnames(refTab[[1]]));
df$condition <- ifelse(
    grepl("(input|control)", df$L1),
    "input",
    "ChIP");


## ------------------------------------------------------------------------
# Plot
cat(sprintf(
    "%s Producing plot..\n",
    ts()));
gg <- ggplot(df, aes(x = x, y = fracTotal));
gg <- gg + geom_bar(aes(fill = condition), stat = "identity");
gg <- gg + facet_wrap(~ L1, nrow = length(unique(df$condition)));
gg <- gg + theme_bw();
gg <- gg + scale_fill_ipsum(name = "Condition");
gg <- gg + labs(
    x = "Reference sequence",
    y = "Fraction of total duplicates",
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
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  plot.title = element_text(size = 16),
  legend.box = "vertical",
  legend.key = element_blank(),
  legend.position = "bottom");
lapply(c("pdf", "png"), function(x) ggsave(
    file = sprintf(
        "%s/plot.distr_dupes_frac.%s",
        outdir,
        x),
    gg,
    bg = "transparent",
    width = 12,
    height = 8));


## ------------------------------------------------------------------------
# Plot
gg <- ggplot(df, aes(x = x, y = log10(DPK)));
gg <- gg + geom_bar(aes(fill = condition), stat = "identity");
gg <- gg + facet_wrap(~ L1, nrow = length(unique(df$condition)));
gg <- gg + theme_bw();
gg <- gg + scale_fill_ipsum(name = "Condition");
gg <- gg + labs(
    x = "Reference sequence",
    y = "log10(DPK)",
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
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  plot.title = element_text(size = 16),
  legend.box = "vertical",
  legend.key = element_blank(),
  legend.position = "bottom");
lapply(c("pdf", "png"), function(x ) ggsave(
    file = sprintf(
        "%s/plot.distr_dupes_DPK.%s",
        outdir,
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

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
suppressMessages(library("jsonlite"));    # Parse JSON
suppressMessages(library("data.table"));  # Read large text files
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
        help = "Gzipped matrix score file
                (from deeptools computeMatrix)",
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
## Custom function
colSe <- function(x) apply(x, 2, function(y) sd(y) / sqrt(length(y)));


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
# Line 1 is JSON header
# Note: read.delim(...) from base R is too slow
cat(sprintf(
    "%s Reading data.\n",
    ts()));
runDef <- readLines(gzfile(input), n = 1);
runDef <- jsonlite::fromJSON(gsub("@", "", runDef));
d <- fread(sprintf("zcat < %s", input), skip = 1L);


## ------------------------------------------------------------------------
# Reformat data
region <- as.data.frame(d[, 1:6]);
score <- as.data.frame(d[, 7:ncol(d)]);
sample.id <- cut(
    seq(1, ncol(score)),
    breaks = runDef$sample_boundaries,
    labels = runDef$sample_labels);
colnames(score) <- as.character(sample.id);
pos <- runDef$`bin size` / 2 + c(
    seq(-runDef$upstream, -1, by = 10),
    seq(0, runDef$body - 1, by = 10),
    seq(runDef$body, runDef$body + runDef$downstream - 1, by = 10));
l <- lapply(levels(sample.id), function(x)
    ret <- score[, which(colnames(score) == x)]);
names(l) <- levels(sample.id);


## ------------------------------------------------------------------------
# Calculate mean score +- s.e.
cat(sprintf(
    "%s Calculating mean scores and uncertainties.\n",
    ts()));
score.meta <- lapply(l, function(x)
    cbind.data.frame(pos = pos, mean = colMeans(x), se = colSe(x)));


## ------------------------------------------------------------------------
# Prepare for profile plotting
df <- melt(score.meta, id.vars = colnames(score.meta[[1]]));
df$L1 <- factor(
    df$L1,
    levels = sort(unique(as.character(df$L1))));
df$condition <- ifelse(
    grepl("(input|control)", df$L1),
    "input",
    "ChIP");
df$l <- df$mean - 1.96 * df$se;
df$h <- df$mean + 1.96 * df$se;


## ------------------------------------------------------------------------
# Plot profile
cat(sprintf(
    "%s Producing plot..\n",
    ts()));
gg <- ggplot(df, aes(x = pos, y = mean));
gg <- gg + geom_ribbon(
    aes(ymin = l, ymax = h, fill = condition),
    alpha = 0.6);
gg <- gg + geom_line(aes(colour = condition));
gg <- gg + facet_wrap(
    ~ L1,
    nrow = length(unique(df$condition)),
    scale = "free_y");
gg <- gg + theme_bw();
gg <- gg + scale_x_continuous(
    breaks = c(-2000, -1000, 0, 1000, 2000, 3000),
    labels = c("-2000", "-1000", "TSS", "TES", "+1000", "+2000"));
gg <- gg + labs(
    x = "Position",
    y = "Score (mean RPKM)",
    title = "Meta-gene score profile",
    subtitle = "Scores are averaged across RefSeq genes; bands denote 95% confidence intervals",
    caption = sprintf(
        "[Maurits Evers (maurits.evers@anu.edu.au), %s]",
        format(Sys.Date(), format = "%B %Y")));
gg <- gg + scale_colour_ipsum(name = "Condition");
gg <- gg + scale_fill_ipsum(name = "Condition");
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
gg <- gg + guides(fill = FALSE);
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


## ------------------------------------------------------------------------
# Clustering

#   df <- l[[1]];
#   colnames(df) <- pos;
#   #quant <- quantile(df, probs = c(0.025, 0.975));
#   #sel <- which(rowMeans(df) >= quant[1] & rowMeans(df) <= quant[2]);
#
#   # Keep rowMeans within the 95% quantile
#   #df <- df[sel, ];
#
#   df[df == 0] <- 1e-3;
#   df <- log10(df);
#
#   getIC <- function(cluster) {
#       m = ncol(cluster$centers);
#       n = length(cluster$cluster);
#       k = nrow(cluster$centers);
#       D = cluster$tot.withinss;
#       return(data.frame(
#           AIC = D + 2 * m * k,
#           BIC = D + log(n) * m * k));
#   }
#   maxk <- 20;
#   km <- lapply(seq(1, maxk), function(x) kmeans(df, centers = x));
#   # Calculate BIC and AIC
#   df.IC <- do.call(rbind, lapply(km, getIC));
#   df.IC$k <- 1:maxk;
#   bestk <- df.IC[which.min(df.IC$BIC), 3];
#   bestk;
#
#
#   # Plot
#   gg <- ggplot(data = df.IC, aes(x = k, y = BIC));
#   gg <- gg + geom_point();
#   gg <- gg + geom_line();
#   gg <- gg + geom_vline(xintercept = bestk, colour = "black");
#   gg <- gg + scale_x_continuous(breaks = pretty(IC$k, n = maxk - 1));
#   gg <- gg + theme_bw();
#   gg <- gg + theme(
#       rect = element_rect(fill = "transparent"),
#       panel.background = element_rect(fill = "transparent"),
#       plot.background = element_rect(fill = "transparent"),
#       strip.background = element_blank(),
#       panel.border = element_rect(colour = "black"),
#       legend.position = "bottom",
#       legend.key = element_blank());
#   gg <- gg + labs(
#       title = "Cluster number k vs. Bayesian Information Criterion (BIC) in k-means clustering");
#   gg;
#
#   df$cluster <- km[[5]]$cluster;
#   df <- df[order(df$cluster, -rowMeans(df[, -ncol(df)])), ];
#
#   df$cluster <- kmeans(df, center = 5)$cluster;
#   df <- df[order(df$cluster), ];
#
#   fcol <- function(x, n = 5){
#       heat.colors(n)[cut(x,n)]
#   }
#   rcol <- t(cbind.data.frame(fcol(df$cluster)));
#   rownames(rcol) <- c("k-means cluster");
#   hm <- heatmap.3(
#       as.matrix(df[, -ncol(df)]),
#       Rowv = NULL,
#       Colv = NULL,
#       dendrogram = "none",
#       trace = "none",
#       scale = "none",
#       RowSideColors = rcol,
#       labRow = TRUE,
#       col = viridis,
#       symbreaks = FALSE,
#       key = TRUE,
#       density.info = "density",
#       margins = c(8, 8),
#       lhei= c(1, 5),
#       lwid = c(1, 5, 0.75));
#   png("plot.png", width = 10, height = 12);
#   eval(hm$call);
#   dev.off();
#
#   #hc <- hclust(dist(tmp, method = "euclidean"));

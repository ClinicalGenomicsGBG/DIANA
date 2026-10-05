#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
    stop("Usage: plot_genomic_regions.R <gviz_data> <sample_id> <bam_file> <output_egfr_coverage> <output_idh1_coverage> <output_idh2_coverage> <output_tertp_coverage> <cytoband_file>")
}

gviz_data_path <- args[1]
sample_id <- args[2]
bam_file <- args[3]
egfr_output <- args[4]
idh1_output <- args[5]
idh2_output <- args[6]
tertp_output <- args[7]
cytoband_file <- args[8]

# Load required libraries
suppressPackageStartupMessages({
    library(Gviz)
    library(GenomicRanges)
    library(Rsamtools)
    library(BSgenome)
    library(BSgenome.Hsapiens.UCSC.hg38)
})

# Verify BAM file exists
if (!file.exists(bam_file)) {
    stop(sprintf("BAM file not found: %s", bam_file))
}

# Detect chromosome naming convention from BAM header once.
# bam_chr() produces the right name for AlignmentsTrack / plotTracks / DataTrack.
# ucsc_chr() always adds "chr" — required by BSgenome.Hsapiens.UCSC.hg38 and
# the cytoband file (both use UCSC naming regardless of BAM convention).
chr_prefix <- tryCatch({
    targets <- names(scanBamHeader(bam_file)[[1]]$targets)
    if (any(grepl("^chr", targets))) "chr" else ""
}, error = function(e) {
    message("Could not read BAM header for chr detection, assuming UCSC (chr): ",
            conditionMessage(e))
    "chr"
})
message("Chromosome prefix detected: '", chr_prefix, "'")

bam_chr  <- function(chrom) paste0(chr_prefix, chrom)
ucsc_chr <- function(chrom) paste0("chr", chrom)

# Load Gviz object
gvizobject <- load(gviz_data_path)

# Function to create custom ideogram with local band data
createCustomIdeogram <- function(chromosome) {
    # Read the local cytoband file from passed argument
    bands_df <- read.table(
        cytoband_file,  # Use passed file path
        sep="\t",
        col.names=c("chrom", "chromStart", "chromEnd", "name", "gieStain"),
        stringsAsFactors=FALSE
    )

    # Filter for the requested chromosome
    bands_df <- bands_df[bands_df$chrom == ucsc_chr(chromosome), ]

    # Create ideogram track with full band information
    itrack <- IdeogramTrack(
        genome = "hg38",
        chromosome = bam_chr(chromosome),
        bands = bands_df,
        showId = TRUE,
        showBandId = TRUE,
        size = 2
    )
    return(itrack)
}

# Plot EGFR coverage
pdf(egfr_output, width=10, height=6)
tryCatch({
    chr7 <- bam_chr("7")
    itrack <- createCustomIdeogram("7")
    gtrack <- GenomeAxisTrack()
    Sample_track <- AlignmentsTrack(bam_file, name = "EGFR")
    ht <- HighlightTrack(
        trackList = Sample_track,
        start = c(55142193, 55154167),
        width = 10,
        chromosome = chr7
    )
    plotTracks(
        list(itrack, gtrack, EGFR_annot, ht),
        from = 55019017,
        to = 55211628,
        chromosome = chr7,
        #main = paste("EGFR Coverage -", sample_id),
        cex = 0.9,
        cex.mismatch = 0.5
    )
}, error = function(e) {
    message("EGFR plot failed: ", conditionMessage(e))
    plot.new()
    text(0.5, 0.5, paste("EGFR plot error:\n", conditionMessage(e)), cex = 0.8)
})
dev.off()

# Plot IDH1 p.R132
pdf(idh1_output, width=10, height=6)
chr2 <- bam_chr("2")
itrack <- createCustomIdeogram("2")
gtrack <- GenomeAxisTrack()
sTrack <- SequenceTrack(Hsapiens, chromosome = ucsc_chr("2"))
Sample_track <- AlignmentsTrack(bam_file, name = "IDH1 p.R132", reverseStacking = TRUE, showMismatch = FALSE)
ht <- HighlightTrack(
    trackList = list(sTrack, Sample_track),
    start = c(208248387),
    width = 2,
    chromosome = chr2
)
plotTracks(
    list(itrack, gtrack, ht),
    chromosome = chr2,
    from = 208248370,
    to = 208248405,
    #main = paste("IDH1 p.R132 -", sample_id),
    cex = 0.9,
    cex.mismatch = 0.5
)
dev.off()

# Plot IDH2 p.R172
pdf(idh2_output, width=10, height=6)
# IDH2 p.R172 (hg38)
chr15 <- bam_chr("15")
itrack <- createCustomIdeogram("15")
gtrack <- GenomeAxisTrack()
sTrack <- SequenceTrack(Hsapiens, chromosome = ucsc_chr("15"))
Sample_track <- AlignmentsTrack(bam_file, name = "IDH2 p.R172", reverseStacking = TRUE, showMismatch = FALSE)

ht <- HighlightTrack(
  trackList = list(sTrack, Sample_track),
  start = c(90088605),   # first base of the R172 codon
  width = 2,             # or use 3 to cover the full codon
  chromosome = chr15
)

plotTracks(
  list(itrack, gtrack, ht),
  chromosome = chr15,
  from = 90088587,       # start - 18 bp
  to   = 90088622,       # start + 17 bp
  cex = 0.9,
  cex.mismatch = 0.5
)
dev.off()


# Plot TERTp
pdf(tertp_output, width=10, height=6)
chr5 <- bam_chr("5")
itrack <- createCustomIdeogram("5")
gtrack <- GenomeAxisTrack()
sTrack <- SequenceTrack(Hsapiens, chromosome = ucsc_chr("5"))
Sample_track <- AlignmentsTrack(bam_file, name = "TERTp", reverseStacking = TRUE)
ht <- HighlightTrack(
    trackList = list(sTrack, Sample_track),
    start = c(1295113, 1295135),
    width = 0,
    chromosome = chr5,
    name = "TERTp"
)
plotTracks(
    list(itrack, gtrack, ht),
    chromosome = chr5,
    from = 1295103,
    to = 1295145,
    #main = paste("TERTp -", sample_id),
    cex = 0.9,
    cex.mismatch = 0.5
)
dev.off()

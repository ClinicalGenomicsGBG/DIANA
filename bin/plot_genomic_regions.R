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
    library(BSgenome)
    library(BSgenome.Hsapiens.UCSC.hg38)
})

normalize_chr <- function(chr) {
    if (grepl("^chr", chr)) chr else paste0("chr", chr)
}

safe_plot_to_pdf <- function(output_file, title_text, plot_expr) {
    pdf(output_file, width=10, height=6)
    ok <- TRUE
    err <- NULL

    tryCatch(
        eval(plot_expr),
        error = function(e) {
            ok <<- FALSE
            err <<- conditionMessage(e)
        }
    )

    if (!ok) {
        plot.new()
        title(main = title_text)
        text(0.5, 0.5, labels = paste("Coverage plot failed:\n", err), cex = 0.8)
    }

    dev.off()
}

# Load Gviz object
gvizobject <- load(gviz_data_path)

# Verify BAM file exists
if (!file.exists(bam_file)) {
    stop(sprintf("BAM file not found: %s", bam_file))
}

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
    chrom <- normalize_chr(chromosome)
    bands_df <- bands_df[bands_df$chrom == chrom, ]
    
    # Create ideogram track with full band information
    itrack <- IdeogramTrack(
        genome = "hg38",
        chromosome = chrom,
        bands = bands_df,
        showId = TRUE,
        showBandId = TRUE,
        size = 2
    )
    return(itrack)
}

# Plot EGFR coverage
safe_plot_to_pdf(egfr_output, paste("EGFR Coverage -", sample_id), quote({
    chr <- "chr7"
    itrack <- createCustomIdeogram(chr)
    gtrack <- GenomeAxisTrack()
    Sample_track <- AlignmentsTrack(bam_file, name = "EGFR")
    ht <- HighlightTrack(
        trackList = Sample_track,
        start = c(55142193, 55154167),
        width = 10,
        chromosome = chr
    )
    plotTracks(
        list(itrack, gtrack, EGFR_annot, ht),
        from = 55019017,
        to = 55211628,
        chromosome = chr,
        cex = 0.9,
        cex.mismatch = 0.5
    )
}))

# Plot IDH1 p.R132
safe_plot_to_pdf(idh1_output, paste("IDH1 p.R132 -", sample_id), quote({
    chr <- "chr2"
    itrack <- createCustomIdeogram(chr)
    gtrack <- GenomeAxisTrack()
    sTrack <- SequenceTrack(Hsapiens, chromosome = chr)
    Sample_track <- AlignmentsTrack(bam_file, name = "IDH1 p.R132", reverseStacking = TRUE)
    ht <- HighlightTrack(
        trackList = list(sTrack, Sample_track),
        start = c(208248387),
        width = 2,
        chromosome = chr
    )
    plotTracks(
        list(itrack, gtrack, ht),
        chromosome = chr,
        from = 208248370,
        to = 208248405,
        cex = 0.9,
        cex.mismatch = 0.5
    )
}))

# Plot IDH2 p.R172
safe_plot_to_pdf(idh2_output, paste("IDH2 p.R172 -", sample_id), quote({
        chr <- "chr15"
        itrack <- createCustomIdeogram(chr)
        gtrack <- GenomeAxisTrack()
        sTrack <- SequenceTrack(Hsapiens, chromosome = chr)
        Sample_track <- AlignmentsTrack(bam_file, name = "IDH2 p.R172", reverseStacking = TRUE)

        ht <- HighlightTrack(
            trackList = list(sTrack, Sample_track),
            start = c(90088605),
            width = 2,
            chromosome = chr
        )

        plotTracks(
            list(itrack, gtrack, ht),
            chromosome = chr,
            from = 90088587,
            to   = 90088622,
            cex = 0.9,
            cex.mismatch = 0.5
        )
}))


# Plot TERTp
safe_plot_to_pdf(tertp_output, paste("TERTp -", sample_id), quote({
    chr <- "chr5"
    itrack <- createCustomIdeogram(chr)
    gtrack <- GenomeAxisTrack()
    sTrack <- SequenceTrack(Hsapiens, chromosome = chr)
    Sample_track <- AlignmentsTrack(bam_file, name = "TERTp", reverseStacking = TRUE)
    ht <- HighlightTrack(
        trackList = list(sTrack, Sample_track),
        start = c(1295113, 1295135),
        width = 0,
        chromosome = chr,
        name = "TERTp"
    )
    plotTracks(
        list(itrack, gtrack, ht),
        chromosome = chr,
        from = 1295103,
        to = 1295145,
        cex = 0.9,
        cex.mismatch = 0.5
    )
}))
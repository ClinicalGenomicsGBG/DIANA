#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
    stop("Usage: plot_genomic_regions_v2.R <gviz_data> <sample_id> <bam_file> <output_egfr_coverage> <output_idh1_coverage> <output_idh2_coverage> <output_tertp_coverage> <cytoband_file>")
}

gviz_data_path <- args[1]
sample_id      <- args[2]
bam_file       <- args[3]
egfr_output    <- args[4]
idh1_output    <- args[5]
idh2_output    <- args[6]
tertp_output   <- args[7]
cytoband_file  <- args[8]

suppressPackageStartupMessages({
    library(Gviz)
    library(GenomicRanges)
    library(Rsamtools)
    library(GenomicAlignments)
})

if (!file.exists(bam_file)) stop(sprintf("BAM file not found: %s", bam_file))

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

# Load BSgenome — optional; fallback omits the SequenceTrack row.
has_bsgenome <- tryCatch({
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
    TRUE
}, error = function(e) {
    message("BSgenome.Hsapiens.UCSC.hg38 not available — sequence row will be omitted: ",
            conditionMessage(e))
    FALSE
})

gvizobject <- load(gviz_data_path)

# Cytoband file uses UCSC chr names; IdeogramTrack chromosome uses BAM convention.
createCustomIdeogram <- function(chrom) {
    bands_df <- read.table(
        cytoband_file,
        sep = "\t",
        col.names = c("chrom", "chromStart", "chromEnd", "name", "gieStain"),
        stringsAsFactors = FALSE
    )
    bands_df <- bands_df[bands_df$chrom == ucsc_chr(chrom), ]
    IdeogramTrack(
        genome     = "hg38",
        chromosome = bam_chr(chrom),
        bands      = bands_df,
        showId     = TRUE,
        showBandId = TRUE,
        size       = 2
    )
}

# ── EGFR ─────────────────────────────────────────────────────────────────────
# Use Rsamtools DataTrack (coverage histogram) instead of AlignmentsTrack so
# sequenceLayer is never triggered for the wide EGFR window.
EGFR_FROM <- 55019017
EGFR_TO   <- 55211628

pdf(egfr_output, width = 10, height = 6)
tryCatch({
    chr7   <- bam_chr("7")
    itrack <- createCustomIdeogram("7")
    gtrack <- GenomeAxisTrack()

    message("EGFR_annot class: ", paste(class(EGFR_annot), collapse = ", "))
    egfr_annot_tracks <- if (inherits(EGFR_annot, "SequenceTrack")) list() else list(EGFR_annot)

    bam_obj  <- Rsamtools::BamFile(bam_file)
    gr       <- GRanges(chr7, IRanges(EGFR_FROM, EGFR_TO))
    egfr_cov <- tryCatch({
        cov <- GenomicAlignments::coverage(bam_obj,
                   param = Rsamtools::ScanBamParam(which = gr))[[chr7]]
        as.numeric(cov[EGFR_FROM:EGFR_TO])
    }, error = function(e) {
        message("EGFR coverage computation failed: ", conditionMessage(e)); NULL
    })

    if (!is.null(egfr_cov)) {
        cov_pos   <- EGFR_FROM:EGFR_TO
        cov_track <- DataTrack(
            data       = egfr_cov,
            start      = cov_pos,
            end        = cov_pos,
            chromosome = chr7,
            genome     = "hg38",
            name       = "EGFR",
            type       = "h",
            col        = "steelblue",
            fill       = "steelblue"
        )
        ht <- HighlightTrack(
            trackList  = cov_track,
            start      = c(55142193, 55154167),
            width      = 10,
            chromosome = chr7
        )
        plotTracks(
            c(list(itrack, gtrack), egfr_annot_tracks, list(ht)),
            from       = EGFR_FROM,
            to         = EGFR_TO,
            chromosome = chr7,
            cex        = 0.9
        )
    } else {
        message("Plotting EGFR annotation only (coverage unavailable)")
        plotTracks(
            c(list(itrack, gtrack), egfr_annot_tracks),
            from       = EGFR_FROM,
            to         = EGFR_TO,
            chromosome = chr7,
            cex        = 0.9
        )
    }
}, error = function(e) {
    message("EGFR plot failed: ", conditionMessage(e))
    plot.new()
    text(0.5, 0.5, paste("EGFR plot error:\n", conditionMessage(e)), cex = 0.8)
})
dev.off()

# ── IDH1 p.R132 ──────────────────────────────────────────────────────────────
# showMismatch=FALSE prevents sequenceLayer from being called — avoids the
# off-limits crash with long ONT reads at this 35 bp window.
# Falls back to pileup without SequenceTrack if BSgenome is unavailable.
pdf(idh1_output, width = 10, height = 6)
tryCatch({
    chr2   <- bam_chr("2")
    itrack <- createCustomIdeogram("2")
    gtrack <- GenomeAxisTrack()

    idh1_plotted <- FALSE
    if (has_bsgenome) {
        tryCatch({
            sTrack       <- SequenceTrack(Hsapiens, chromosome = ucsc_chr("2"))
            Sample_track <- AlignmentsTrack(bam_file, name = "IDH1 p.R132",
                                            reverseStacking = TRUE,
                                            showMismatch    = FALSE)
            ht <- HighlightTrack(
                trackList  = list(sTrack, Sample_track),
                start      = c(208248387),
                width      = 2,
                chromosome = chr2
            )
            plotTracks(
                list(itrack, gtrack, ht),
                chromosome = chr2,
                from       = 208248370,
                to         = 208248405,
                type       = "pileup",
                cex        = 0.9
            )
            idh1_plotted <- TRUE
        }, error = function(e) {
            message("IDH1 SequenceTrack attempt failed, using pileup fallback: ",
                    conditionMessage(e))
        })
    }
    if (!idh1_plotted) {
        Sample_track <- AlignmentsTrack(bam_file, name = "IDH1 p.R132",
                                        reverseStacking = TRUE)
        ht <- HighlightTrack(
            trackList  = Sample_track,
            start      = c(208248387),
            width      = 2,
            chromosome = chr2
        )
        plotTracks(
            list(itrack, gtrack, ht),
            chromosome = chr2,
            from       = 208248370,
            to         = 208248405,
            type       = "pileup",
            cex        = 0.9
        )
    }
}, error = function(e) {
    message("IDH1 plot failed: ", conditionMessage(e))
    plot.new()
    text(0.5, 0.5, paste("IDH1 plot error:\n", conditionMessage(e)), cex = 0.8)
})
dev.off()

# ── IDH2 p.R172 ──────────────────────────────────────────────────────────────
pdf(idh2_output, width = 10, height = 6)
tryCatch({
    chr15  <- bam_chr("15")
    itrack <- createCustomIdeogram("15")
    gtrack <- GenomeAxisTrack()

    idh2_plotted <- FALSE
    if (has_bsgenome) {
        tryCatch({
            sTrack       <- SequenceTrack(Hsapiens, chromosome = ucsc_chr("15"))
            Sample_track <- AlignmentsTrack(bam_file, name = "IDH2 p.R172",
                                            reverseStacking = TRUE,
                                            showMismatch    = FALSE)
            ht <- HighlightTrack(
                trackList  = list(sTrack, Sample_track),
                start      = c(90088605),
                width      = 2,
                chromosome = chr15
            )
            plotTracks(
                list(itrack, gtrack, ht),
                chromosome = chr15,
                from       = 90088587,
                to         = 90088622,
                type       = "pileup",
                cex        = 0.9
            )
            idh2_plotted <- TRUE
        }, error = function(e) {
            message("IDH2 SequenceTrack attempt failed, using pileup fallback: ",
                    conditionMessage(e))
        })
    }
    if (!idh2_plotted) {
        Sample_track <- AlignmentsTrack(bam_file, name = "IDH2 p.R172",
                                        reverseStacking = TRUE)
        ht <- HighlightTrack(
            trackList  = Sample_track,
            start      = c(90088605),
            width      = 2,
            chromosome = chr15
        )
        plotTracks(
            list(itrack, gtrack, ht),
            chromosome = chr15,
            from       = 90088587,
            to         = 90088622,
            type       = "pileup",
            cex        = 0.9
        )
    }
}, error = function(e) {
    message("IDH2 plot failed: ", conditionMessage(e))
    plot.new()
    text(0.5, 0.5, paste("IDH2 plot error:\n", conditionMessage(e)), cex = 0.8)
})
dev.off()

# ── TERTp ─────────────────────────────────────────────────────────────────────
pdf(tertp_output, width = 10, height = 6)
tryCatch({
    chr5   <- bam_chr("5")
    itrack <- createCustomIdeogram("5")
    gtrack <- GenomeAxisTrack()

    tertp_plotted <- FALSE
    if (has_bsgenome) {
        tryCatch({
            sTrack       <- SequenceTrack(Hsapiens, chromosome = ucsc_chr("5"))
            Sample_track <- AlignmentsTrack(bam_file, name = "TERTp",
                                            reverseStacking = TRUE,
                                            showMismatch    = FALSE)
            ht <- HighlightTrack(
                trackList  = list(sTrack, Sample_track),
                start      = c(1295113, 1295135),
                width      = 0,
                chromosome = chr5,
                name       = "TERTp"
            )
            plotTracks(
                list(itrack, gtrack, ht),
                chromosome = chr5,
                from       = 1295103,
                to         = 1295145,
                type       = "pileup",
                cex        = 0.9
            )
            tertp_plotted <- TRUE
        }, error = function(e) {
            message("TERTp SequenceTrack attempt failed, using pileup fallback: ",
                    conditionMessage(e))
        })
    }
    if (!tertp_plotted) {
        Sample_track <- AlignmentsTrack(bam_file, name = "TERTp",
                                        reverseStacking = TRUE)
        ht <- HighlightTrack(
            trackList  = Sample_track,
            start      = c(1295113, 1295135),
            width      = 0,
            chromosome = chr5,
            name       = "TERTp"
        )
        plotTracks(
            list(itrack, gtrack, ht),
            chromosome = chr5,
            from       = 1295103,
            to         = 1295145,
            type       = "pileup",
            cex        = 0.9
        )
    }
}, error = function(e) {
    message("TERTp plot failed: ", conditionMessage(e))
    plot.new()
    text(0.5, 0.5, paste("TERTp plot error:\n", conditionMessage(e)), cex = 0.8)
})
dev.off()

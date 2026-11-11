#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GenomicFeatures); library(GenomicRanges); library(GenomicAlignments)
  library(Rsamtools); library(BiocParallel); library(rtracklayer); library(optparse)
})
option_list <- list(
  make_option("--txdb_sqlite", type="character"),
  make_option("--rmsker", type="character"),
  make_option("--output", type="character", default="repeat_annotations.rds"),
  make_option("--cpus", type="integer", default=1)
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$txdb_sqlite) || is.null(opt$rmsker)) stop("txdb_sqlite and rmsker required")
cpus <- max(1L, as.integer(opt$cpus))
txdb <- loadDb(opt$txdb_sqlite)
seq_infos=seqinfo(txdb)

read_repeatmasker <- function(path, cpus) {
    if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) {
        stop("RepeatMasker file does not exist: ", path)
    }

    if (grepl("\\.gz$", path, ignore.case=TRUE)) {
        cmd <- paste("gzip -dc --", shQuote(path))
        x <- data.table::fread(cmd=cmd, sep=" ", nThread=cpus, skip=3, header=FALSE, fill=TRUE)
    } else {
        x <- data.table::fread(path, sep=" ", nThread=cpus, skip=3, header=FALSE, fill=TRUE)
    }

    if (nrow(x) == 0L) {
        stop("RepeatMasker input contains no annotation rows after the three-line header: ", path)
    }

    # Standard RepeatMasker .out has 15 fields. Some valid rows carry a trailing
    # '*' overlap marker, yielding a 16th field. Preserve the 15 annotation fields
    # and accept the optional marker rather than rejecting valid RepeatMasker output.
    if (ncol(x) < 15L || ncol(x) > 16L) {
        stop("RepeatMasker input must contain 15 fields plus at most the optional trailing '*' marker; observed ",
             ncol(x), " columns in: ", path)
    }
    if (ncol(x) == 16L) {
        marker <- trimws(as.character(x[[16L]]))
        bad_marker <- !is.na(marker) & nzchar(marker) & marker != "*"
        if (any(bad_marker)) {
            stop("Unexpected non-'*' values in RepeatMasker optional column 16: ", path)
        }
        x <- x[, seq_len(15L), with=FALSE]
    }

    x
}

rmsk = read_repeatmasker(opt$rmsker, cpus)
colnames(rmsk) <- c( 
    "swScore", "milliDiv", "milliDel", "milliIns", "genoName", "genoStart", "genoEnd", "genoLeft", "strand", "repName", "repClassFamily","repStart", "repEnd", "repLeft", "id"
    )

# Fail early on malformed query coordinates instead of silently constructing invalid
# genomic intervals. RepeatMasker query coordinates are 1-based inclusive.
rmsk$genoStart <- suppressWarnings(as.integer(rmsk$genoStart))
rmsk$genoEnd <- suppressWarnings(as.integer(rmsk$genoEnd))
if (anyNA(rmsk$genoStart) || anyNA(rmsk$genoEnd) ||
    any(rmsk$genoStart < 1L) || any(rmsk$genoEnd < rmsk$genoStart)) {
    stop("RepeatMasker input contains invalid query coordinates: ", opt$rmsker)
}
rmsk$strand[rmsk$strand == "C"] <- "-"
rmksGR = makeGRangesFromDataFrame(  rmsk,
                                    seqnames.field="genoName",
                                    start.field="genoStart",
                                    end.field="genoEnd",
                                    strand.field="strand",
                                    keep.extra.columns=TRUE
                                    )
# Keep only RepeatMasker intervals whose chromosome naming matches the TxDb.
# A complete mismatch (for example chr1 versus 1) must fail loudly: otherwise the
# pipeline would report zero repeat signal for a reference-annotation mismatch.
n_rmsk_input <- length(rmksGR)
keep_seqlevel <- as.character(seqnames(rmksGR)) %in% seqlevels(seq_infos)
n_rmsk_kept <- sum(keep_seqlevel)
if (n_rmsk_kept == 0L) {
    stop("No RepeatMasker intervals match TxDb seqlevels. Check genome build and chromosome naming.")
}
if (n_rmsk_kept < n_rmsk_input) {
    warning(sprintf("Dropping %d/%d RepeatMasker intervals on seqlevels absent from the TxDb",
                    n_rmsk_input - n_rmsk_kept, n_rmsk_input))
}
rmksGR <- rmksGR[keep_seqlevel]
seqlevels(rmksGR) <- seqlevels(seq_infos)
seqinfo(rmksGR) = seq_infos
rmksGR = trim(rmksGR)

names(rmksGR) = paste0(seqnames(rmksGR),":",start(rmksGR),"-",end(rmksGR),"_",strand(rmksGR),"_",rmksGR$repClassFamily,"_",rmksGR$repName)

rmksGR$id_row = names(rmksGR)
## UP/DOWN: relative to the elements
# we get the 

# Symmetric genomic flank exposure used only as a local-background diagnostic.
# It is deliberately NOT interpreted as transcript-direction upstream/downstream.
# Each side targets 30% of repeat width, is non-overlapping with the repeat body,
# and is clipped to chromosome boundaries.
w_width = pmax(1L, as.integer(round(width(rmksGR) * 0.3)))

left_end <- start(rmksGR) - 1L
left_keep <- left_end >= 1L
rmksGR_upstream <- rmksGR[left_keep]
if (length(rmksGR_upstream) > 0L) {
    left_width <- w_width[left_keep]
    left_end_kept <- left_end[left_keep]
    ranges(rmksGR_upstream) <- IRanges::IRanges(
        start=pmax(1L, left_end_kept - left_width + 1L),
        end=left_end_kept
    )
}

chrom_len <- seqlengths(seq_infos)[as.character(seqnames(rmksGR))]
right_start <- end(rmksGR) + 1L
right_keep <- is.na(chrom_len) | right_start <= chrom_len
rmksGR_downstream <- rmksGR[right_keep]
if (length(rmksGR_downstream) > 0L) {
    right_start_kept <- right_start[right_keep]
    right_end <- right_start_kept + w_width[right_keep] - 1L
    known_len <- !is.na(chrom_len[right_keep])
    right_end[known_len] <- pmin(right_end[known_len], chrom_len[right_keep][known_len])
    ranges(rmksGR_downstream) <- IRanges::IRanges(
        start=right_start_kept,
        end=right_end
    )
}

seqlevels(rmksGR_upstream) <- seqlevels(seq_infos)
seqinfo(rmksGR_upstream) = seq_infos
seqlevels(rmksGR_downstream) <- seqlevels(seq_infos)
seqinfo(rmksGR_downstream) = seq_infos

# Combine up downs
rmksGR_updown = c(rmksGR_upstream, rmksGR_downstream)
rmksGR_updown = rmksGR_updown[order(rmksGR_updown$id_row)]

# convert to a list
rmksGR_updown = split(rmksGR_updown,rmksGR_updown$id_row)
rmksGR_updown = GRangesList(rmksGR_updown)



saveRDS(list(rmksGR=rmksGR, rmksGR_updown=rmksGR_updown), file=opt$output)

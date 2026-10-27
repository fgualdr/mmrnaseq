#!/usr/bin/env Rscript

library("GenomicFeatures")
library("GenomicAlignments")
library(GenomicRanges)
library("Rsamtools")
library("BiocParallel")
library("openxlsx")
library("tidyverse")
library("yaml")
library(biomaRt)
library(rtracklayer)
library(optparse)

# Cusom functions:
`%!in%` <- Negate(`%in%`)

# R function to parse parameters using optparse.
# This will take in: 
# 1) GTF file path 
# 2) output file path
# 3) import n cores cpus

# Aim is to workout the GTF file to combine by Tx and Gene and generate Exon GRagnes, All_features GRanges, and Intron+UTRs GRanges
# The output will be a GRanges object with the following columns:
option_list <- list(
    make_option(c("-q", "--txdb_sqlite"), type="character", default=NULL, help="sqlite txdb file path", metavar="character"),
    make_option(c("-r", "--rmsker"), type="character", default=NULL, help="rmsker out file path", metavar="character"),
    make_option(c("-m", "--meta"), type="character", default=NULL, help="meta file txdb", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, help="Output file path", metavar="character"),
    make_option(c("-i", "--id"), type="character", default=NULL, help="sample_id", metavar="character"),
    make_option(c("-b", "--bam"), type="character", default=NULL, help="BAM file path", metavar="character"),
    make_option(c("-s", "--strand"), type="character", default=NULL, help="Sample strand", metavar="character"),
    make_option(c("-p", "--paired_end"), type="character", default=NULL, help="single / paired", metavar="character"),
    make_option(c("-c", "--cpus"), type="integer", default=1, help="Import n cores cpus", metavar="integer")
)

opt <- parse_args(OptionParser(option_list=option_list))

# opt <- list(
#     txdb_sqlite = "./TxDB.sqlite",
#     rmsker = "./hg38.fa.out.gz",
#     meta = "./metadata.txt",
#     output = "./",
#     id = "DONOR3_24h_SKAP2_KO_R2",
#     bam = "./DONOR3_24h_SKAP2_KO_R2.em.mark.filt.sorted.bam",
#     strand = "reverse",
#     paired_end = "paired",
#     cpus = 10
# )

cpus = as.numeric(opt$cpus)

# Load the txdb_sqlite file
txdb <- loadDb(opt$txdb_sqlite)
Meta <- read.delim(opt$meta, header=TRUE, sep="\t")

TxbyGene <- transcriptsBy(txdb,by="gene")
ExonByGene <- exonsBy(txdb,by="gene")
#  The introns are the difference between the exons and the transcripts
IntronByGene <- intronsByTranscript(txdb, use.names=TRUE)
UTR5ByTranscript <- fiveUTRsByTranscript(txdb, use.names=TRUE)
UTR3ByTranscript <- threeUTRsByTranscript(txdb, use.names=TRUE)

# Now we count the reads in the exons, introns and transcripts - consider the type of strandness

bamfile <- BamFile(opt$bam, yieldSize=1000000)
os <- Sys.info()["sysname"]
if (os == "Linux") {
    print("The system is Linux.")
    BiocParallel::register(MulticoreParam(workers=cpus, 
                                 tasks=0, 
                                 stop.on.error=TRUE, 
                                 progressbar=TRUE))
} else if (os == "Darwin") {
    print("The system is macOS (Darwin).")
    BiocParallel::register(MulticoreParam(workers=cpus, 
                                 tasks=0, 
                                 stop.on.error=TRUE, 
                                 progressbar=TRUE))
} else if (os == "Windows") {
    print("The system is Windows.")
    # For Windows, we use SnowParam with type "SOCK"
    BiocParallel::register(SnowParam(workers=cpus, 
                                 tasks=0, 
                                 stop.on.error=TRUE, 
                                 progressbar=TRUE, 
                                 type="SOCK"))
} else {
    print(paste("Unknown system:", os))
    BiocParallel::register(SnowParam(workers=cpus, 
                                 tasks=0, 
                                 stop.on.error=TRUE, 
                                 progressbar=TRUE, 
                                 type="SOCK"))
}

# In summarizeOverlaps, we set preprocess.reads to invertStrand, which will reverse the strand of each read before counting it when:
# protocols that lead to the strand of the first read being the opposite to that of the original RNA (1st family): dUTP, NSR, NNSR, Illumina stranded TruSeq PE protocol
# Else protocols that lead to the strand of the first read being that of the original RNA (2nd family): Directional Illumina (Ligation), Standard SOLiD 

# Define the strandness and single/paired end
if(opt$paired_end == "single"){
    singleEnd <- TRUE
    ignorestrand <- FALSE
}else{
    singleEnd <- FALSE
    ignorestrand <- FALSE
}

if(opt$strand == "reverse"){
    strandfun <- invertStrand
}else{
    if(opt$strand == "non_specific"){
        strandfun <- NULL
        ignorestrand <- TRUE
    }else{
        if(opt$strand == "forward"){
            strandfun <- NULL
    }
    }
}


# Count reads in exons, introns and transcripts
se_all <- summarizeOverlaps(features=TxbyGene,
    reads=bamfile,
    mode="Union",
    singleEnd=singleEnd,
    ignore.strand=ignorestrand, 
    preprocess.reads=strandfun,
    inter.feature=FALSE)

se_int <- summarizeOverlaps(features=IntronByGene,
    reads=bamfile,
    mode="IntersectionNotEmpty",
    singleEnd=singleEnd,
    ignore.strand=ignorestrand, 
    preprocess.reads=strandfun,
    inter.feature=FALSE)

se_exonic <- summarizeOverlaps(features=ExonByGene,
    reads=bamfile,
    mode="IntersectionStrict",
    singleEnd=singleEnd,
    ignore.strand=ignorestrand, 
    preprocess.reads=strandfun,
    inter.feature=FALSE)

se_5_utr <- summarizeOverlaps(features=UTR5ByTranscript,
    reads=bamfile,
    mode="Union",
    singleEnd=singleEnd,
    ignore.strand=ignorestrand, 
    preprocess.reads=strandfun,
    inter.feature=FALSE)

se_3_utr <- summarizeOverlaps(features=UTR3ByTranscript,
    reads=bamfile,
    mode="Union",
    singleEnd=singleEnd,
    ignore.strand=ignorestrand, 
    preprocess.reads=strandfun,
    inter.feature=FALSE)

# Now we can extract the counts from the summarized experiment

dat_all <- as.data.frame(assay(se_all)) %>%
    rownames_to_column( var = "gene_id" )

dat_int <- as.data.frame(assay(se_int)) %>%
    rownames_to_column( var = "gene_id" )

dat_ex <- as.data.frame(assay(se_exonic)) %>%
    rownames_to_column( var = "gene_id" )

dat_5_utr <- as.data.frame(assay(se_5_utr)) %>%
    rownames_to_column( var = "gene_id" )

dat_3_utr <- as.data.frame(assay(se_3_utr)) %>%
    rownames_to_column( var = "gene_id" )

colnames(dat_all) <- colnames(dat_int) <- colnames(dat_ex) <- colnames(dat_5_utr) <- colnames(dat_3_utr) <- c("gene_id", opt$id)
rownames(dat_all) <- dat_all$gene_id
rownames(dat_int) <- dat_int$gene_id
rownames(dat_ex) <- dat_ex$gene_id
rownames(dat_5_utr) <- dat_5_utr$gene_id
rownames(dat_3_utr) <- dat_3_utr$gene_id

# We save the individual counts
dir.create(opt$output, recursive=TRUE, showWarnings=FALSE)
# Export a named list to an rdata file
res_list = list(dat_all=dat_all, dat_int=dat_int, dat_ex=dat_ex, dat_5_utr=dat_5_utr, dat_3_utr=dat_3_utr)
saveRDS(res_list, file=paste(opt$output, opt$id, ".gene.rds", sep=""))



## We now read in the rmsker file and add the counts to the metadata
# skip the first 3 lines of the file

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


## if the data is stranded we count the reads in the sense and antisense strand and add this information
if(!ignorestrand){
        
        as_fun <- if(opt$strand == "reverse") {
            NULL
        } else {
            invertStrand
        }

        sense_reps <- summarizeOverlaps(features = rmksGR, 
                                    reads = bamfile, 
                                    mode="IntersectionStrict",
                                    ignore.strand = ignorestrand, 
                                    singleEnd = singleEnd, 
                                    inter.feature = FALSE, 
                                    preprocess.reads = strandfun)

        asense_reps <- summarizeOverlaps(features = rmksGR, 
                                    reads = bamfile, 
                                    mode="IntersectionStrict",
                                    ignore.strand = ignorestrand, 
                                    singleEnd = singleEnd, 
                                    inter.feature = FALSE, 
                                    preprocess.reads= as_fun)

        ## same for the flanking regions
        sense_reps_updown <- summarizeOverlaps(features = rmksGR_updown, 
                                    reads = bamfile, 
                                    mode="IntersectionStrict",
                                    ignore.strand = ignorestrand, 
                                    singleEnd = singleEnd, 
                                    inter.feature = FALSE, 
                                    preprocess.reads = strandfun)

        asense_reps_updown <- summarizeOverlaps(features = rmksGR_updown,
                                    reads = bamfile, 
                                    mode="IntersectionStrict",
                                    ignore.strand = ignorestrand, 
                                    singleEnd = singleEnd, 
                                    inter.feature = FALSE, 
                                    preprocess.reads= as_fun)
        
        sense_reps <- as.data.frame(assay(sense_reps)) 
        asense_reps <- as.data.frame(assay(asense_reps)) 
        sense_reps_updown <- as.data.frame(assay(sense_reps_updown)) 
        asense_reps_updown <- as.data.frame(assay(asense_reps_updown)) 

        colnames(sense_reps) <- paste0(opt$id,"_sense")
        colnames(asense_reps) <- paste0(opt$id,"_asense")
        colnames(sense_reps_updown) <- paste0(opt$id,"_sense_flank")
        colnames(asense_reps_updown) <- paste0(opt$id,"_asense_flank")
        
        # add width information
        sense_reps$width_ss = width(rmksGR[rownames(sense_reps)])
        asense_reps$width_as = width(rmksGR[rownames(asense_reps)])
        
        # Use the actual available flank exposure after chromosome-boundary clipping.
        flank_width <- vapply(
            rmksGR_updown[rownames(sense_reps_updown)],
            function(x) sum(width(x)),
            numeric(1)
        )
        sense_reps_updown$width_flank_ss = flank_width
        asense_reps_updown$width_flank_as = flank_width

        dat_rmsk <- list(
            sense_count_reps = sense_reps,
            asense__count_reps = asense_reps,
            sense_count_flank = sense_reps_updown,
            asense_count_flank = asense_reps_updown
        )
}else{

    se_all_reps <- summarizeOverlaps(features=rmksGR,
                                    reads=bamfile,
                                    mode="Union",
                                    singleEnd=singleEnd,
                                    ignore.strand=ignorestrand, 
                                    preprocess.reads=strandfun,
                                    inter.feature=FALSE)

    ## same for the flanking regions
    se_all_reps_updown <- summarizeOverlaps(features=rmksGR_updown,
                                            reads=bamfile,
                                            mode="Union",
                                            singleEnd=singleEnd,
                                            ignore.strand=ignorestrand,
                                            preprocess.reads=strandfun,
                                            inter.feature=FALSE)

    # here we get the counts and make the table but will be "counts", "counts_updown"
    dat_rmsk <- as.data.frame(assay(se_all_reps))
    dat_rmsk_updown <- as.data.frame(assay(se_all_reps_updown))

    colnames(dat_rmsk) <- paste0(opt$id,"_counts")
    colnames(dat_rmsk_updown) <- paste0(opt$id,"_counts_flank")
    dat_rmsk$width = width(rmksGR[rownames(dat_rmsk)])
    dat_rmsk_updown$width_flank = vapply(
        rmksGR_updown[rownames(dat_rmsk_updown)],
        function(x) sum(width(x)),
        numeric(1)
    )
    dat_rmsk <- list(
        count_reps = dat_rmsk,
        count_flank = dat_rmsk_updown
    )

}

saveRDS(dat_rmsk, file=paste(opt$output, opt$id, ".rmsk.rds", sep=""))
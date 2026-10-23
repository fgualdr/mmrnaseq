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
dir.create(opt$output, recursive=TRUE)
# Export a named list to an rdata file
res_list = list(dat_all=dat_all, dat_int=dat_int, dat_ex=dat_ex, dat_5_utr=dat_5_utr, dat_3_utr=dat_3_utr)
saveRDS(res_list, file=paste(opt$output, opt$id, ".gene.rds", sep=""))



## We now read in the rmsker file and add the counts to the metadata
# skip the first 3 lines of the file

seq_infos=seqinfo(txdb)
rmsk = data.table::fread(cmd=paste("zcat", opt$rmsker), sep=" ",nThread=cpus,skip=3,header=FALSE)
colnames(rmsk) <- c( 
    "swScore", "milliDiv", "milliDel", "milliIns", "genoName", "genoStart", "genoEnd", "genoLeft", "strand", "repName", "repClassFamily","repStart", "repEnd", "repLeft", "id"
    )
rmsk$strand[rmsk$strand == "C"] <- "-"
rmksGR = makeGRangesFromDataFrame(  rmsk,
                                    seqnames.field="genoName",
                                    start.field="genoStart",
                                    end.field="genoEnd",
                                    strand.field="strand",
                                    keep.extra.columns=TRUE
                                    )
# Keep seqlevels only in the txdb
rmksGR <- rmksGR[seqnames(rmksGR) %in% seqlevels(seq_infos)]
# Drop unused seqlevels
seqlevels(rmksGR) <- seqlevels(seq_infos)
seqinfo(rmksGR) = seq_infos
rmksGR = trim(rmksGR)

names(rmksGR) = paste0(seqnames(rmksGR),":",start(rmksGR),"-",end(rmksGR),"_",strand(rmksGR),"_",rmksGR$repClassFamily,"_",rmksGR$repName)

rmksGR$id_row = names(rmksGR)
## UP/DOWN: relative to the elements
# we get the 

w_width = width(rmksGR) # consider 30% of the width of the element so either upstream of it or downstream of it
w_width = round(w_width * 0.3)

rmksGR_upstream = rmksGR
start(rmksGR_upstream) = start(rmksGR_upstream) - w_width
end(rmksGR_upstream) = start(rmksGR_upstream) 

rmksGR_downstream = rmksGR
end(rmksGR_downstream) = end(rmksGR_downstream) + w_width
start(rmksGR_downstream) = end(rmksGR_downstream)

seqlevels(rmksGR_upstream) <- seqlevels(seq_infos)
seqinfo(rmksGR_upstream) = seq_infos
rmksGR_upstream = trim(rmksGR_upstream)

seqlevels(rmksGR_downstream) <- seqlevels(seq_infos)
seqinfo(rmksGR_downstream) = seq_infos
rmksGR_downstream = trim(rmksGR_downstream)

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
        
        # for the flanking should be the length of the elements * 0.3 *2
        w <- width(rmksGR[rownames(sense_reps_updown)])
        sense_reps_updown$width_flank_ss = (w*0.3*2)
        asense_reps_updown$width_flank_as = (w*0.3*2)

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
    dat_rmsk_updown$width_flank = (dat_rmsk$width*0.3*2)
    dat_rmsk <- list(
        count_reps = dat_rmsk,
        count_flank = dat_rmsk_updown
    )

}

saveRDS(dat_rmsk, file=paste(opt$output, opt$id, ".rmsk.rds", sep=""))
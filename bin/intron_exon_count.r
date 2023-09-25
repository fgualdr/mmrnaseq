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
    make_option(c("-m", "--meta"), type="character", default=NULL, help="meta file txdb", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, help="Output file path", metavar="character"),
    make_option(c("-i", "--id"), type="character", default=NULL, help="sample_id", metavar="character"),
    make_option(c("-b", "--bam"), type="character", default=NULL, help="BAM file path", metavar="character"),
    make_option(c("-s", "--strand"), type="character", default=NULL, help="Sample strand", metavar="character"),
    make_option(c("-p", "--paired_end"), type="character", default=NULL, help="single / paired", metavar="character"),
    make_option(c("-c", "--cpus"), type="integer", default=1, help="Import n cores cpus", metavar="integer")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Load the txdb_sqlite file
txdb <- loadDb(opt$txdb_sqlite)
Meta <- read.delim(opt$meta, header=TRUE, sep="\t")

TxbyGene <- transcriptsBy(txdb,by="gene")
ExonByGene <- exonsBy(txdb,by="gene")
#  The introns are the difference between the exons and the transcripts
IntronByGene <- intronsByTranscript(txdb, use.names=TRUE)

# Now we count the reads in the exons, introns and transcripts - consider the type of strandness

bamfile <- BamFile(opt$bam, yieldSize=1000000)

cpus=opt$cores
if(!is.null(cpus)){
    param <- BiocParallel::SnowParam(workers=cpus,tasks=0,stop.on.error=TRUE,progressbar = TRUE,type="SOCK")
}else{
    param <- NULL
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
    inter.feature=FALSE,
    param=param)

se_int <- summarizeOverlaps(features=IntronByGene,
    reads=bamfile,
    mode="IntersectionNotEmpty",
    singleEnd=singleEnd,
    ignore.strand=ignorestrand, 
    preprocess.reads=strandfun,
    inter.feature=FALSE,
    param=param)

se_exonic <- summarizeOverlaps(features=ExonByGene,
    reads=bamfile,
    mode="IntersectionStrict",
    singleEnd=singleEnd,
    ignore.strand=ignorestrand, 
    preprocess.reads=strandfun,
    inter.feature=FALSE,
    param=param)

# Now we can extract the counts from the summarized experiment

dat_all <- as.data.frame(assay(se_all)) %>%
    rownames_to_column( var = "gene_id" )

dat_int <- as.data.frame(assay(se_int)) %>%
    rownames_to_column( var = "gene_id" )

dat_ex <- as.data.frame(assay(se_exonic)) %>%
    rownames_to_column( var = "gene_id" )

colnames(dat_all) <- colnames(dat_int) <- colnames(dat_ex) <- c("gene_id", opt$id)
rownames(dat_all) <- dat_all$gene_id
rownames(dat_int) <- dat_int$gene_id
rownames(dat_ex) <- dat_ex$gene_id

# We save the individual counts
dir.create(opt$output, recursive=TRUE)
# Export a named list to an rdata file
res_list = list(dat_all=dat_all, dat_int=dat_int, dat_ex=dat_ex)
saveRDS(res_list, file=paste(opt$output, opt$id, ".rds", sep=""))

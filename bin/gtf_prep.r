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

# R function to parse parameters using optparse.
# This will take in: 
# 1) GTF file path 
# 2) output file path
# 3) import n cores cpus

# qsub -I -q nocg_interact -l select=1:ncpus=5:mem=80g
# singularity shell -B /hpcnfs docker://fgualdr/envr_norm_gr
# Aim is to workout the GTF file to combine by Tx and Gene and generate Exon GRagnes, All_features GRanges, and Intron+UTRs GRanges
# The output will be a GRanges object with the following columns:
option_list <- list(
    make_option(c("-g", "--gtf"), type="character", default=NULL, help="GTF file path", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, help="Output file path", metavar="character"),
    make_option(c("-s", "--chrom_size"), type="character", default=NULL, help="Chromosome size file path", metavar="character"),
    make_option(c("-c", "--cpus"), type="integer", default=1, help="Import n cores cpus", metavar="integer")
)
opt <- parse_args(OptionParser(option_list=option_list))

# To use this code from command line: 
# Rscript gtf_prep.r 
# opt <- list()
# opt$gtf = "/hpcnfs/data/GN2/fgualdrini/UsefullData/DB/Assembly/ncbi_dataset/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf"
# opt$out = "./"
# opt$chrom_size = "/hpcnfs/scratch/temporary/fgualdr_ccle_test/work/f7/e714dd3b767835c4c5f0ab1ce46ef2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fa.sizes"
# opt$cpus = 5
# -g /hpcnfs/data/GN2/fgualdrini/UsefullData/DB/Assembly/ncbi_dataset/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf 
# -o './' 
# -s /hpcnfs/scratch/temporary/fgualdr_ccle_test/work/f7/e714dd3b767835c4c5f0ab1ce46ef2/GCF_009914755.1_T2T-CHM13v2.0_genomic.fa.sizes
# -c 5

# Cusom functions:
`%!in%` <- Negate(`%in%`)

# Worout the name of the mitochondrial chromosome from the GTF file
# So we first identify by exclusion autosomes using two criterias:
# The value in the first column of the GTF (excluding the lines beginning with "#") should not be:
# 1) they should not match these patterns exactMatch = c("y", "w", "z", "m", "mt", "mtdna", "pt", "mtr", "2-micron", "ebv", "nc_000087.8","nc_005089.1")
# 2) they should not partially match with fuzzyMatch = c("chrUn", "random", "v1", "v2", "kn707", "jtfh","nw_","nt_")
# Having selected the autosomes we can remove them and the mitochondrial chromosome should be the one that are not partially matching the fuzzyMatch

exactMatch = tolower(c("chrM","y", "w", "z", "m", "mt", "mtdna", "pt", "mtr", "2-micron", "ebv", "nc_000087.8","nc_005089.1"))
fuzzyMatch = tolower(c("chrUn", "random", "v1", "v2", "kn707", "jtfh","nw_","nt_"))
sex = tolower(c("chrX", "chrY","x","y","nc_060947.1","nc_060948.1"))

# Read the GTF file
gtf <- read.table(opt$gtf, sep="\t", header=FALSE, comment.char="#", stringsAsFactors=FALSE)
chr = unique(gtf[,1])
# Remove matching
autosome <- chr[!(tolower(chr) %in% exactMatch)]
autosome <- autosome[!grepl(paste(fuzzyMatch, collapse="|"), tolower(autosome))]

mito <- chr[!(tolower(chr) %in% tolower(autosome))]
mito <- mito[!grepl(paste(fuzzyMatch, collapse="|"), tolower(mito))]
# Remove sex chromosomes
mito <- mito[!grepl(paste(sex, collapse="|"), tolower(mito))]

# Import the GTF file using the TxDb package
chrom_channel = read.delim(opt$chrom_size,header=FALSE,sep="\t",stringsAsFactors=FALSE)
colnames(chrom_channel) = c("chrom","length")
gtf = gtf[gtf[,1] %in% chrom_channel[,1],]
write.table(gtf, file = "./gtf.gtf", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

if(length(mito) > 0){
    txdb <- makeTxDbFromGFF("./gtf.gtf", format="auto", circ_seqs = mito,chrominfo=chrom_channel)
}else{
    txdb <- makeTxDbFromGFF("./gtf.gtf", format="auto",chrominfo=chrom_channel)
}

# From the GTF get the Metadata stored in the column 9 strplit by ";"
metaDat <- strsplit(gtf[,9], "; ")
# lapply and split each by the space " ".
# split by the space " " and make the first element the name of the vector
# Perform parallelization
cpus = opt$cpus
param <- BiocParallel::SnowParam(workers=cpus,tasks=0,stop.on.error=TRUE,progressbar = TRUE,type="SOCK")

metaDatV <- BiocParallel::bplapply(metaDat, function(x) {
    # remove trailing space and ;
    x <- gsub(" $", "",gsub(";$", "", x))
    # Split x by the space " "
    x <- do.call(rbind,strsplit(x, " "))
    n = x[,1]
    y = x[,2]
    names(y) = n
    return(y)
}, BPPARAM=param)

nn <- BiocParallel::bplapply(metaDatV,names, BPPARAM=param)
nn = unique(unlist(nn))
metaDatV <- BiocParallel::bplapply(metaDatV, function(x) {
    x = x[nn]
    return(x)
}, BPPARAM=param)

# Convert the list to a data frame
metaDatVM <- do.call(rbind, metaDatV)
colnames(metaDatVM) = nn

# combine the first 8 columns with the metadata - the first 8 are mandatory while the 9th are the attributes
metaDatVM = cbind(gtf[,1:8], metaDatVM)
colnames(metaDatVM)[1:8] = c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame")
# Convert to GRanges keep metadata

##############################################################################################################
# Reduce/split Transcripts by gene

Transcripts_by <- GenomicFeatures::transcriptsBy(txdb, by="gene")
Exons_by <- GenomicFeatures::exonsBy(txdb, by="gene")
CDS_by <- GenomicFeatures::cdsBy(txdb, by="gene")
Five_UTR_by <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names=TRUE)
Three_UTR_by <- GenomicFeatures::threeUTRsByTranscript(txdb, use.names=TRUE)

Transcripts_by_red <- GenomicRanges::reduce(Transcripts_by)
# Get the number of transcripts per gene
multi_tx = lapply(Transcripts_by_red,length)
w = which(unlist(multi_tx) > 1)
# add the suffix to those genes with multiple non-overlapping transcripts
Multi_tx = Transcripts_by_red[w]
Multi_tx = unlist(Multi_tx)

unique_names <- unique(names(Multi_tx))

# Separate transcripts, exons, cds and UTRs belonging to those genes without multiple non-overlapping transcripts
Transcripts_by_single =  Transcripts_by[which(names(Transcripts_by) %!in% unique_names)]
Exons_by_single =  Exons_by[which(names(Exons_by) %!in% unique_names)]
CDS_by_single =  CDS_by[which(names(CDS_by) %!in% unique_names)]

# find the column in metaDatVM that contains the gene_id elements in unique_names - metaDatVM is a GRanges so we need to convert it to a data frame
col_gene_id = apply(metaDatVM, 2, function(x) {
    z = which(x %in% unique_names)
    return(length(z))
})
col_gene_id = unlist(col_gene_id)
w_g = which.max(col_gene_id)
metaDatVM_single = metaDatVM[metaDatVM[,w_g] %!in% unique_names,]
sel_gene_col = grep("gene_",colnames(metaDatVM_single))
sel_gene_col = c(1:8,sel_gene_col)
metaDatVM_single_gene = metaDatVM_single[,sel_gene_col]
metaDatVM_single_gene = metaDatVM_single_gene[metaDatVM_single_gene$feature == "gene",]


# Add the gene_id to the transcripts
Transcripts_by_single = unlist(Transcripts_by_single)
Transcripts_by_single$gene_id = names(Transcripts_by_single)
Transcripts_by_single = split(Transcripts_by_single, Transcripts_by_single$gene_id)

# Add the gene_id to the exons
Exons_by_single = unlist(Exons_by_single)
Exons_by_single$gene_id = names(Exons_by_single)
Exons_by_single = split(Exons_by_single, Exons_by_single$gene_id)

# Add the gene_id to the cds
CDS_by_single = unlist(CDS_by_single)
CDS_by_single$gene_id = names(CDS_by_single)
CDS_by_single = split(CDS_by_single, CDS_by_single$gene_id)

# Only fix the gene name:
for (i in unique_names) {   

    cat("Processing gene: ", i,"\n")

    # select the involved tx - names will relate to the gene_id
    dd = Multi_tx[names(Multi_tx) == i]
    dd$gene_id = paste0(names(dd),"::",1:length(dd))
    names(dd) = dd$gene_id

    cat(" fix the metadata\n")
    # Overlapping MetaData:
    overlapping_meta = metaDatVM[which(metaDatVM[,w_g] == i ) ,]
    overlapping_meta = overlapping_meta[ overlapping_meta$feature == "gene",]
    gr_overlapping_meta = GRanges(seqnames = overlapping_meta$seqnames, ranges = IRanges(start = overlapping_meta$start, end = overlapping_meta$end), strand = overlapping_meta$strand)
    fov = findOverlaps(gr_overlapping_meta,dd)
    overlapping_meta = overlapping_meta[queryHits(fov),]
    rownames(overlapping_meta) = names(dd)[subjectHits(fov)]
    overlapping_meta$gene_id = names(dd)[subjectHits(fov)]
    # append to the metaDatVM_single
    overlapping_meta = overlapping_meta[,sel_gene_col]
    metaDatVM_single_gene = rbind(metaDatVM_single_gene, overlapping_meta)

    cat(" fix the transcripts\n")
    # Overlappoing tx:
    overlapping_tx = unlist(Transcripts_by[i])
    fov = findOverlaps(overlapping_tx, dd)
    names(overlapping_tx) = names(dd)[subjectHits(fov)]
    overlapping_tx$gene_id = names(dd)[subjectHits(fov)]
    overlapping_tx =  split(overlapping_tx, overlapping_tx$gene_id)
    # Append to the Transcripts_by_single
    Transcripts_by_single = c(Transcripts_by_single, overlapping_tx)

    cat(" fix the exons\n")
    # Assign the new gene_id to the overlapping exons, cds and UTRs
    if(i %in% names(Exons_by)){
        overlpping_exons = unlist(Exons_by[i])
        fov = findOverlaps(overlpping_exons, dd)
        names(overlpping_exons) = names(dd)[subjectHits(fov)]
        overlpping_exons$gene_id = names(dd)[subjectHits(fov)]
        overlpping_exons =  split(overlpping_exons,overlpping_exons$gene_id)
        # Append to the Exons_by_single
        Exons_by_single = c(Exons_by_single, overlpping_exons)
    }

    cat(" fix the cds\n")
    # Assign the new gene_id to the overlapping exons, cds
    if(i %in% names(CDS_by )){
        overlapping_cds = unlist(CDS_by[i])
        fov = findOverlaps(overlapping_cds, dd)
        names(overlapping_cds) = names(dd)[subjectHits(fov)]
        overlapping_cds$gene_id = names(dd)[subjectHits(fov)]
        overlapping_cds =  split(overlapping_cds, overlapping_cds$gene_id)
        # Append to the CDS_by_single
        CDS_by_single = c(CDS_by_single, overlapping_cds)
    }

}

# Now when reducing the GRanges features we will have the correct transcripts:
# So one Transcript per gene_id
# By reducing the Exons we will get the pseudo-exons per transcriot/gene
# 
# Add back the gene_id:
Transcripts <- unlist(GenomicRanges::reduce(Transcripts_by_single))
Transcripts$gene_id = names(Transcripts)
Transcripts = split(Transcripts, Transcripts$gene_id)

# Add back the gene_id:
Exons <- unlist(GenomicRanges::reduce(Exons_by_single))
Exons$gene_id = names(Exons)
Exons = split(Exons, Exons$gene_id)

# Add back the gene_id:
CDS <- unlist(GenomicRanges::reduce(CDS_by_single))
CDS$gene_id = names(CDS)
CDS = split(CDS, CDS$gene_id)
CDS_Rage = range(CDS_by_single)


transcripts <- as.data.frame(Transcripts)
colnames(transcripts) <- c("tx_id","tx_name","tx_chrom","tx_start","tx_end","width","tx_strand","gene_id")
transcripts$tx_id <- 1:dim(transcripts)[1]
transcripts$gene_id <- transcripts$tx_name
transcripts$tx_id <- as.integer(transcripts$tx_id)
transcripts$tx_start <- as.integer(transcripts$tx_start)
transcripts$tx_end <- as.integer(transcripts$tx_end)
transcripts$tx_name <- as.character(transcripts$tx_name)
transcripts$tx_chrom <- as.character(transcripts$tx_chrom)
transcripts$gene_id <- as.character(transcripts$gene_id)
transcripts$tx_strand <- as.character(transcripts$tx_strand)

## Make a ditionary of Tx_id and Tx_name
Dictio <- cbind(transcripts$tx_id,transcripts$tx_name)
colnames(Dictio) <- c("tx_id","tx_name")
rownames(Dictio) <- transcripts$gene_id
Dictio <- as.data.frame(Dictio)
DictioT <- Dictio$tx_id
names(DictioT) <- Dictio$tx_name

# Exons
ExonByGene <- unlist(Exons)
names(ExonByGene) = NULL
ExonByGene$tx_id <- as.numeric(as.character(DictioT[ExonByGene$gene_id]))

tx_id_all <- split(ExonByGene$tx_id,ExonByGene$tx_id)
# so we generate a progression number for each element in exon_rank
# so if exon rank has: 1,1,1,1,2,2,2,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4
# it has to turn into  1,2,3,4,1,2,3,1,2,3,4,5,6,7,8,9,1,2,3,4,5,6
exon_rank <- lapply(tx_id_all, function(x) {
    y = 1:length(x)
    return(y)
})
exon_rank <- unlist(exon_rank)
exon_rank = as.data.frame(cbind(ExonByGene$tx_id,exon_rank))
ExonByGene$exon_rank <- exon_rank$exon_rank

ExonByGene$cds_start <- NA
ExonByGene$cds_end <- NA

# Add CDS for loop with CDS information
combo <- c(ExonByGene,unlist(CDS_Rage))
combo <- disjoin(combo)
hits.1 <- findOverlaps(combo,CDS_Rage,type="within")
combo <- combo[queryHits(hits.1)]
hits.1 <- findOverlaps(combo,ExonByGene,type="within")
ExonByGene[subjectHits(hits.1)]$cds_start <- start(combo)[queryHits(hits.1)]
ExonByGene[subjectHits(hits.1)]$cds_end <- end(combo)[queryHits(hits.1)]
# Remake as GRangeList
ExonByGene_split <- split(ExonByGene, ExonByGene$tx_id)


splicings_db <- as.data.frame(ExonByGene_split)
splicings_db$exon_id <- 1:nrow(splicings_db)
splicings <- splicings_db[,c("tx_id","exon_rank","exon_id","seqnames","strand","start","end","cds_start","cds_end")]
colnames(splicings) <- c("tx_id","exon_rank","exon_id","exon_chrom","exon_strand","exon_start","exon_end","cds_start","cds_end")

splicings$tx_id <- as.integer(as.character(splicings$tx_id))
splicings$exon_rank <- as.integer(as.character(splicings$exon_rank))
splicings$exon_id <- as.integer(as.character(splicings$exon_id))
splicings$exon_chrom <- as.character(splicings$exon_chrom)
splicings$exon_strand <- as.character(splicings$exon_strand)
splicings$exon_start <- as.integer(as.character(splicings$exon_start))
splicings$exon_end <- as.integer(as.character(splicings$exon_end))


# Chrom chrominfo(mm10Txdb)
# Need to fix the chrominfo
chrominfo_fix <- as.data.frame(seqlevels(Transcripts_by))

chrominfo_fix$length <- NA
chrominfo_fix$is_circular <- FALSE
colnames(chrominfo_fix) <- c("chrom","length","is_circular")
# Find the mito chromosome by partial match as gtf can be from ucsc ensembl or refseq
chrominfo_fix$is_circular[grepl(paste0("chrm$|chrmt$|chrM$|MT$|mito"),chrominfo_fix$chrom)] <- TRUE
chrom_channel = read.delim(opt$chrom_size,header=FALSE,sep="\t",row.names=1)
chrominfo_fix$length <- chrom_channel[as.character(chrominfo_fix$chrom),1]

# Generate the final Txdb
txdb_fix <- makeTxDb(transcripts, splicings,chrominfo_fix , genes=NULL)
# zsave as sqlite
saveDb(txdb_fix, file=paste0(opt$out,"TxDB.sqlite"))
# Write the metadata as a tab delimited file
write.table(metaDatVM_single_gene,file=paste0(opt$out,"metadata.txt"),sep="\t",quote=FALSE,row.names=FALSE)

# rtracklayer::export.gff3(txdb_fix,con=paste0(opt$out,"TxDB.gff"),version="3")
# test <- import(paste0(opt$out,"TxDB.gff"))
# rtracklayer::export(test,paste0(opt$out,"TxDB.gtf"),"gtf")

# TxbyGene <- transcriptsBy(txdb_fix,by="gene")

# # Ge the Conding and NonCoding
# CDS_by <- GenomicFeatures::cdsBy(txdb_fix, by="gene")
# ExonByGene <- exonsBy(txdb_fix,by="gene")
# IntronByGene <- psetdiff(unlist(TxbyGene),CDS_by,ignore.strand=FALSE)
# # Get the UTRs 5' and 3'
# Five_UTR_by <- GenomicFeatures::fiveUTRsByTranscript(txdb_fix, use.names=TRUE)
# Three_UTR_by <- GenomicFeatures::threeUTRsByTranscript(txdb_fix, use.names=TRUE)
# # Convert to lists the UTRs. Then merge by element name and strand

# UTRs <- unlist(c(Five_UTR_by,Three_UTR_by))
# UTRs <- split(UTRs, names(UTRs))
# UTRs <- GRangesList(UTRs)


# # Save the TxbyGene, ExonByGene and NonCodingByGene as RData
# save(TxbyGene,file=paste0(opt$out,"TxbyGene.RData"))
# save(ExonByGene,file=paste0(opt$out,"ExonByGene.RData"))
# save(IntronByGene,file=paste0(opt$out,"IntronByGene.RData"))

# # Export as BED files
# rtracklayer::export(unlist(TxbyGene),con=paste0(opt$out,"TxbyGene.bed"))
# rtracklayer::export(unlist(ExonByGene),con=paste0(opt$out,"ExonByGene.bed"))
# rtracklayer::export(unlist(IntronByGene),con=paste0(opt$out,"IntronByGene.bed"))


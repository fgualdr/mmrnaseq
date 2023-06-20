#!/usr/bin/env Rscript

################################################
################################################
## REQUIREMENTS                               ##
################################################
################################################

## The Script expect in the working directory all the .rds files related to the individual samples
## test: /hpcnfs/scratch/temporary/fgualdr_smarca_rna/work/test
# singularity shell -B /hpcnfs docker://fgualdr/envnorm_fix


################################################
################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(DESeq2)
library(dplyr)
library(tidyr)
library(GeneralNormalizer)
library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)

suppressPackageStartupMessages(library(ComplexHeatmap))

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################

option_list <- list(
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory."                                                                      ),
    make_option(c("-p", "--outprefix"     ), type="character", default='custom_norm_deseq2', metavar="string" , help="Output prefix."                                                                         ),
    make_option(c("-v", "--vst"           ), type="logical"  , default=FALSE   , metavar="boolean", help="Run vst transform instead of rlog."                                                     ),
    make_option(c("-n", "--norm"          ), type="logical"  , default=FALSE   , metavar="boolean", help="Normalize with the General Normalizer Library."                                         ),
    make_option(c("-s", "--sigma_times"   ), type="integer"  , default=1       , metavar="integer", help="Number of sigma times."                                                                 ),
    make_option(c("-z", "--n_pop"          ), type="integer"  , default=1       , metavar="integer", help="Number of populations."                                                                 ),
    make_option(c("-c", "--cores"         ), type="integer"  , default=1       , metavar="integer", help="Number of cores."                                                                       )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

################################################
rdat_files <- list.files(pattern="*.rds")

if (length(rdat_files) == 0) {
    print_help(opt_parser)
    stop("No RDS files.", call.=FALSE)
}

################################################
################################################
## READ IN COUNTS FILE                        ##
################################################
################################################

# Import the rds files load the R objects:
count.table <- list()
for (i in 1:length(rdat_files)) {
    count.table[[i]] <- readRDS(rdat_files[i])
}
names(count.table) = gsub(".rds","",rdat_files)

# The imported list is composed of a list of 3 data.frames i.e. "dat_all" "dat_int" "dat_ex"
# We generate the dat_all, dat_int, dat_ex tables by merging the corresponding tables from the list
# We extract the tables named "dat_all" "dat_int" "dat_ex" from the list - there are as many tables as samples
dat_all = list()
dat_int = list()
dat_ex  = list()
for (i in 1:length(count.table)) {
    dat_all[[i]] = count.table[[i]]$dat_all
    dat_int[[i]] = count.table[[i]]$dat_int
    dat_ex[[i]]  = count.table[[i]]$dat_ex
}
# We reduce the list to a single data.frame by merging the tables
dat_all = Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), dat_all)
dat_int = Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), dat_int)
dat_ex  = Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), dat_ex)
# We set the rownames to the gene_id
rownames(dat_all) = dat_all$gene_id
rownames(dat_int) = dat_int$gene_id
rownames(dat_ex)  = dat_ex$gene_id
# Remove the gene_id column
dat_all = dat_all[,!colnames(dat_all) %in% c("gene_id")]
dat_int = dat_int[,!colnames(dat_int) %in% c("gene_id")]
dat_ex  = dat_ex[,!colnames(dat_ex) %in% c("gene_id")]

# Ok now the dat_all will be used to compute the normalization factors so we start with that


count.table           <- dat_all
colnames(count.table) <- gsub(pattern='\\.$', replacement='', colnames(count.table))


################################################
################################################
## RUN General NORM                           ##
################################################
################################################

cpus=opt$cores
if(!is.null(cpus)){
    param <- BiocParallel::SnowParam(workers=cpus,tasks=0,stop.on.error=TRUE,progressbar = TRUE,type="SOCK")
}else{
    param <- NULL
}

if(opt$norm == TRUE){

    save_folder_norm = paste0(opt$outdir,"/normalization/")
    dir.create(save_folder_norm, recursive=TRUE)

    # Assemble the DEG_DESIGN using name Layers by "_"
    # Names w/o _R[0-9]
    cn =  colnames(count.table)
    # Now the conditions we need to strtify - 
    name_components <- strsplit(cn, "_")
    n_components    <- length(unlist(name_components))
    decompose       <- n_components!=1 && all(sapply(name_components, length)==n_components)
    coldata         <- data.frame(cn, Sample_ID=cn, row.names=1)
    if (decompose) {
        groupings        <- as.data.frame(lapply(1:n_components, function(i) sapply(name_components, "[[", i)))
        names(groupings) <- paste0("Group", 1:n_components)
        n_distinct       <- sapply(groupings, function(grp) length(unique(grp)))
        groupings        <- groupings[n_distinct!=1 & n_distinct!=length(cn)]
        if (ncol(groupings)!=0) {
            coldata <- cbind(coldata, groupings)
        } else {
            decompose <- FALSE
        }
    }else{
        coldata$Group1 = gsub("_R[0-9].*","",cn)
    }

    # Take only the R[0-9] part so remove text before any _R[0-9]
    coldata$Sample_Replicate = paste0("R",gsub(".*_R","",colnames(count.table)))

    design = coldata

    x = count.table
    x_depth_p = x/colSums(x)*1000000
    # remove raws top 0.1% and bottom 0.1% based on rowMeans
    x = x[order(rowMeans(x_depth_p), decreasing = T),]
    x = x[round(nrow(x)*0.001):round(nrow(x)*0.999),]

    # Compute the normalization
    sigma_times = opt$sigma_times
    n_pop =    opt$n_pop
    Result = RunNorm(   x,
                        design,
                        sigma_times = sigma_times,
                        fix_reference="random",
                        row_name_index=rownames(x),
                        saving_path=save_folder_norm,
                        n_pop=n_pop,
                        n_pop_reference=n_pop,
                        BiocParam=param
                        )
    
    colData = Result$scaling_factors
    mat = Result$norm_mat

}

################################################
################################################
## RUN DESEQ2                                 ##
################################################
################################################

if (file.exists(opt$outdir) == FALSE) {
    dir.create(opt$outdir, recursive=TRUE)
}
setwd(opt$outdir)

samples.vec     <- colnames(count.table)
name_components <- strsplit(samples.vec, "_")
n_components    <- length(name_components[[1]])
decompose       <- n_components!=1 && all(sapply(name_components, length)==n_components)
coldata         <- data.frame(samples.vec, sample=samples.vec, row.names=1)
if (decompose) {
    groupings        <- as.data.frame(lapply(1:n_components, function(i) sapply(name_components, "[[", i)))
    n_distinct       <- sapply(groupings, function(grp) length(unique(grp)))
    groupings        <- groupings[n_distinct!=1 & n_distinct!=length(samples.vec)]
    if (ncol(groupings)!=0) {
        names(groupings) <- paste0("Group", 1:ncol(groupings))
        coldata <- cbind(coldata, groupings)
    } else {
        decompose <- FALSE
    }
}

DDSFile <- paste(opt$outprefix,".dds.RData",sep="")

counts  <- count.table[,samples.vec,drop=FALSE]
dds     <- DESeqDataSetFromMatrix(countData=round(counts), colData=coldata, design=~ 1)
dds     <- estimateSizeFactors(dds)

if(opt$norm == TRUE){
    norm = Result$scaling_factors
    sizeFactors(dds) = norm[names(sizeFactors(dds)),"scaling"]
}

if (min(dim(count.table))<=1)  { # No point if only one sample, or one gene
    save(dds,file=DDSFile)
    saveRDS(dds, file=sub("\\.dds\\.RData$", ".rds", DDSFile))
    warning("Not enough samples or genes in counts file for PCA.", call.=FALSE)
    quit(save = "no", status = 0, runLast = FALSE)
}
if (!opt$vst) {
    vst_name <- "rlog"
    rld      <- rlog(dds)
} else {
    vst_name <- "vst"
    rld      <- varianceStabilizingTransformation(dds)
}

assay(dds, vst_name) <- assay(rld)
save(dds,file=DDSFile)
saveRDS(dds, file=sub("\\.dds\\.RData$", ".rds", DDSFile))

################################################
################################################
## PLOT QC                                    ##
################################################
################################################

##' PCA pre-processeor
##'
##' Generate all the necessary information to plot PCA from a DESeq2 object
##' in which an assay containing a variance-stabilised matrix of counts is
##' stored. Copied from DESeq2::plotPCA, but with additional ability to
##' say which assay to run the PCA on.
##'
##' @param object The DESeq2DataSet object.
##' @param ntop number of top genes to use for principla components, selected by highest row variance.
##' @param assay the name or index of the assay that stores the variance-stabilised data.
##' @return A data.frame containing the projected data alongside the grouping columns.
##' A 'percentVar' attribute is set which includes the percentage of variation each PC explains,
##' and additionally how much the variation within that PC is explained by the grouping variable.
##' @author Gavin Kelly
plotPCA_vst <- function (object,  ntop = 500, assay=length(assays(object))) {
    rv         <- rowVars(assay(object, assay))
    select     <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca        <- prcomp(t(assay(object, assay)[select, ]), center=TRUE, scale=FALSE)
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    df         <- cbind( as.data.frame(colData(object)), pca$x)
    #Order points so extreme samples are more likely to get label
    ord        <- order(abs(rank(df$PC1)-median(df$PC1)), abs(rank(df$PC2)-median(df$PC2)))
    df         <- df[ord,]
    attr(df, "percentVar") <- data.frame(PC=seq(along=percentVar), percentVar=100*percentVar)
    return(df)
}

PlotFile <- paste(opt$outprefix,".plots.pdf",sep="")

pdf(file=PlotFile, onefile=TRUE, width=7, height=7)
## PCA
ntop <- c(500, Inf)
for (n_top_var in ntop) {
    pca.data      <- plotPCA_vst(dds, assay=vst_name, ntop=n_top_var)
    percentVar    <- round(attr(pca.data, "percentVar")$percentVar)
    plot_subtitle <- ifelse(n_top_var==Inf, "All genes", paste("Top", n_top_var, "genes"))
    pl <- ggplot(pca.data, aes(PC1, PC2, label=paste0(" ", sample, " "))) +
        geom_point() +
        geom_text(check_overlap=TRUE, vjust=0.5, hjust="inward") +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        labs(title = paste0("First PCs on ", vst_name, "-transformed data"), subtitle = plot_subtitle) +
        theme(legend.position="top",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    print(pl)

    if (decompose) {
        pc_names <- paste0("PC", attr(pca.data, "percentVar")$PC)
        long_pc <- reshape(pca.data, varying=pc_names, direction="long", sep="", timevar="component", idvar="pcrow")
        long_pc <- subset(long_pc, component<=5)
        long_pc_grp <- reshape(long_pc, varying=names(groupings), direction="long", sep="", timevar="grouper")
        long_pc_grp <- subset(long_pc_grp, grouper<=5)
        long_pc_grp$component <- paste("PC", long_pc_grp$component)
        long_pc_grp$grouper <- paste0(long_pc_grp$grouper, c("st","nd","rd","th","th")[long_pc_grp$grouper], " prefix")
        pl <- ggplot(long_pc_grp, aes(x=Group, y=PC)) +
            geom_point() +
            stat_summary(fun=mean, geom="line", aes(group = 1)) +
            labs(x=NULL, y=NULL, subtitle = plot_subtitle, title="PCs split by sample-name prefixes") +
            facet_grid(component~grouper, scales="free_x") +
            scale_x_discrete(guide = guide_axis(n.dodge = 3))
        print(pl)
    }
} # at end of loop, we'll be using the user-defined ntop if any, else all genes

## WRITE PC1 vs PC2 VALUES TO FILE
pca.vals           <- pca.data[,c("PC1","PC2")]
colnames(pca.vals) <- paste0(colnames(pca.vals), ": ", percentVar[1:2], '% variance')
pca.vals           <- cbind(sample = rownames(pca.vals), pca.vals)
write.table(pca.vals, file = paste(opt$outprefix, ".pca.vals.txt", sep=""),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = TRUE)

## SAMPLE CORRELATION HEATMAP
sampleDists      <- dist(t(assay(dds, vst_name)))
sampleDistMatrix <- as.matrix(sampleDists)
colors           <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(
    sampleDistMatrix,
    clustering_distance_rows=sampleDists,
    clustering_distance_cols=sampleDists,
    col=colors,
    main=paste("Euclidean distance between", vst_name, "of samples")
)

## WRITE SAMPLE DISTANCES TO FILE
write.table(cbind(sample = rownames(sampleDistMatrix), sampleDistMatrix),file=paste(opt$outprefix, ".sample.dists.txt", sep=""),
            row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
dev.off()

################################################
################################################
## SAVE SIZE FACTORS                          ##
################################################
################################################

SizeFactorsDir <- "size_factors/"
if (file.exists(SizeFactorsDir) == FALSE) {
    dir.create(SizeFactorsDir, recursive=TRUE)
}

NormFactorsFile <- paste(SizeFactorsDir,opt$outprefix, ".size_factors.RData", sep="")

normFactors <- sizeFactors(dds)
save(normFactors, file=NormFactorsFile)

for (name in names(sizeFactors(dds))) {
    sizeFactorFile <- paste(SizeFactorsDir,name, ".txt", sep="")
    write(as.numeric(sizeFactors(dds)[name]), file=sizeFactorFile)
}


# colData must have "Sample_ID" and "scaling" columns
# sizeFactors(dds) are with the custom normalization or the DESEQ2 normalization
# So we assemble the colData with the sizeFactors(dds) and write it to file
colData <- data.frame(Sample_ID=colnames(assay(dds, vst_name)), scaling=as.numeric(sizeFactors(dds)))
rownames(colData) <- colData$Sample_ID
write.table(colData,file=paste(opt$outdir, "/scaling_dat.txt", sep=""),
        row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)


norm_mat <- dat_all[, colData$Sample_ID]
norm_mat <- sweep(norm_mat[, colData$Sample_ID], 2, colData[, "scaling"], "*")
write.table(norm_mat,file=paste(opt$outdir, "/ALL_READS_norm.txt", sep=""),
        row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

norm_mat <- dat_int[, colData$Sample_ID]
norm_mat <- sweep(norm_mat[, colData$Sample_ID], 2, colData[, "scaling"], "*")
write.table(norm_mat,file=paste(opt$outdir, "/INTRON_READS_norm.txt", sep=""),
        row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)

norm_mat <- dat_ex[, colData$Sample_ID]
norm_mat <- sweep(norm_mat[, colData$Sample_ID], 2, colData[, "scaling"], "*")
write.table(norm_mat,file=paste(opt$outdir, "/EXON_READS_norm.txt", sep=""),
        row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)


################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

RLogFile <- "R_sessionInfo.log"

sink(RLogFile)
a <- sessionInfo()
print(a)
sink()

################################################
################################################
################################################
################################################

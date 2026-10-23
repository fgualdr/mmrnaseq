#!/usr/bin/env Rscript

# ============================================================
# Omics Normalization & QC Pipeline with DESeq2 + OmniNorm
# Compatible with Linux SLURM HPC via MulticoreParam
# Author: [Your Name]
# ============================================================

# ------------------------------------------------------------
# LOAD LIBRARIES
# ------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(OmniNorm)
  library(ComplexHeatmap)
  library(ggplot2)
  library(RColorBrewer)
  library(BiocParallel)
})

set.seed(12345)

# ------------------------------------------------------------
# COMMAND-LINE OPTIONS
# ------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--count_file"),     type="character", help="Path to count matrix file."),
  make_option(c("-d", "--id_col"),         type="integer",   default=1, help="Column with gene/feature IDs."),
  make_option(c("-f", "--count_col"),      type="integer",   default=3, help="Index of first sample column."),
  make_option(c("-r", "--sample_suffix"),  type="character", default="", help="Suffix to strip from sample names."),
  make_option(c("-o", "--outdir"),         type="character", default="./", help="Output directory."),
  make_option(c("-p", "--outprefix"),      type="character", default="deseq2", help="Output filename prefix."),
  make_option(c("-v", "--vst"),            type="logical",   default=FALSE, help="Use vst (TRUE) or rlog (FALSE)."),
  make_option(c("-n", "--norm"),           type="logical",   default=FALSE, help="Use OmniNorm normalization."),
  make_option(c("-s", "--sigma_times"),    type="integer",   default=1, help="Sigma multiplier for invariant filtering."),
  make_option(c("-z", "--n_pop"),          type="integer",   default=1, help="Number of populations to fit."),
  make_option(c("-c", "--cores"),          type="integer",   default=1, help="Number of parallel cores.")
)

opt <- parse_args(OptionParser(option_list=option_list))


if (is.null(opt$count_file)) {
  stop("Error: --count_file is required.", call.=FALSE)
}

# ------------------------------------------------------------
# INPUT & SAMPLE PROCESSING
# ------------------------------------------------------------
count_table <- read.delim(opt$count_file, header=TRUE, row.names=NULL)
rownames(count_table) <- count_table[[opt$id_col]]
count_table <- count_table[, opt$count_col:ncol(count_table), drop=FALSE]

# Clean up sample names
colnames(count_table) <- gsub(opt$sample_suffix, "", colnames(count_table))
colnames(count_table) <- gsub("\\.$", "", colnames(count_table))

sample_names <- colnames(count_table)

if (any(duplicated(colnames(count_table)))) {
    stop("Sample names (column names) must be unique after stripping suffix.")
}

# ------------------------------------------------------------
# MULTICORE CONFIGURATION (SLURM-compatible)
# ------------------------------------------------------------
cpus  <- opt$cores
bpp   <- if (cpus > 1) MulticoreParam(cpus,progressbar=TRUE) else SerialParam()

# ------------------------------------------------------------
# RUN NORMALIZATION (OmniNorm)
# ------------------------------------------------------------
if (opt$norm) {
    norm_outdir <- file.path(opt$outdir, "normalization")
    dir.create(norm_outdir, recursive=TRUE, showWarnings=FALSE)

    # Construct design with required columns
    Sample_IDs <- colnames(count_table)

    Sample_Condition <- gsub("_R[0-9]+.*", "", Sample_IDs)
    Sample_Replicate <- gsub(".*_R", "R", Sample_IDs)

    design_df <- data.frame(
    Sample_ID        = Sample_IDs,
    Sample_Condition = Sample_Condition,
    Sample_Replicate = Sample_Replicate,
    stringsAsFactors = FALSE
    )
    rownames(design_df) <- design_df$Sample_ID

    if (any(!grepl("^R[0-9]+$", design_df$Sample_Replicate))) {
        warning("Some Sample_Replicate entries do not match expected R[0-9]+ pattern.")
    }

    # Normalize only well-behaved rows (depth-normalized and trimmed)
    x <- count_table
    x_depth <- x / colSums(x) * 1e6
    x <- x[order(rowMeans(x_depth), decreasing=TRUE), ]
    x <- x[round(nrow(x)*0.001):round(nrow(x)*0.999), ]

    # Call RunNorm with updated API
    Result <- RunNorm(
        mat_path = x,
        design_path = design_df[colnames(x), , drop=FALSE],
        fix_reference = "random",
        row_name_index = NULL,
        n_pop = opt$n_pop,
        n_pop_reference = opt$n_pop,
        saving_path = norm_outdir,
        sigma_times = opt$sigma_times,
        dist_family = "Skew.normal",
        Norm_plot = TRUE,
        Save_results = TRUE,
        return_full = FALSE,
        BiocParam = bpp,
        cpus = cpus
    )
    sc <- Result$scaling_factors
    sc[is.na(sc$scaling) | sc$scaling <= 0, "scaling"] <- 1
    sc[is.na(sc$mean_g) | sc$mean_g <= 0, "mean_g"] <- 1
    
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

samples.vec     <- colnames(count_table)
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

if(opt$norm == TRUE){
    # ingest pre-norm table:
    counts  <- count_table
    dds     <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~ 1)
    sizeFactors(dds) = 1/sc[colnames(counts), "mean_g"]
} else {
    counts  <- count_table[,samples.vec,drop=FALSE]
    dds     <- DESeqDataSetFromMatrix(countData=round(counts), colData=coldata, design=~ 1)
    dds     <- estimateSizeFactors(dds)
}

if (min(dim(counts))<=1)  { # No point if only one sample, or one gene
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
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
        )
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

hc <- hclust(sampleDists, method = "complete")
ht <- Heatmap(
  sampleDistMatrix,
  name = "Distance",
  col = colors,
  cluster_rows = as.dendrogram(hc),
  cluster_columns = as.dendrogram(hc),
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  column_title = paste("Euclidean distance between", vst_name, "of samples"),
  row_names_side = "left"
)
draw(ht)

# Close the PDF device opened above
if (grDevices::dev.cur() > 1) grDevices::dev.off()
# End of script
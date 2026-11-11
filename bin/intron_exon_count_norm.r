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
  library(matrixStats)  # <-- add this
})

set.seed(12345)


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
rdat_files <- list.files(pattern="\\.rds$")  # proper regex
if (length(rdat_files) == 0) {
    print_help(opt_parser)
    stop("No RDS files.", call.=FALSE)
}
# read
count.table <- lapply(rdat_files, readRDS)
names(count.table) <- gsub("\\.rds$", "", rdat_files)  # escape dot, anchor

# extract the three tables per sample
dat_all <- lapply(count.table, `[[`, "dat_all")
dat_int <- lapply(count.table, `[[`, "dat_int")
dat_ex  <- lapply(count.table, `[[`, "dat_ex")

# merge by gene_id (outer join), then set rownames and drop gene_id col
merge_by_gene <- function(lst){
  out <- Reduce(function(x,y) merge(x,y, by="gene_id", all=TRUE), lst)
  if (anyDuplicated(out$gene_id)) stop("Duplicate gene IDs in count input")
  rownames(out) <- out$gene_id
  out <- out[, setdiff(colnames(out), "gene_id"), drop=FALSE]
  # replace NA with 0 (DESeq2 cannot take NA)
  out[is.na(out)] <- 0
  # enforce numeric
  out[] <- lapply(out, function(z) { z <- as.numeric(z); ifelse(is.na(z), 0, z) })
  as.data.frame(out, check.names=FALSE)
}

dat_all <- merge_by_gene(dat_all)
dat_int <- merge_by_gene(dat_int)
dat_ex  <- merge_by_gene(dat_ex)

count.table <- dat_all
colnames(count.table) <- gsub("\\.$", "", colnames(count.table))

################################################
################################################
## RUN General NORM                           ##
################################################
################################################

cpus=opt$cores
if(!is.null(cpus)){
    param <- BiocParallel::MulticoreParam(workers=cpus, progressbar = TRUE)
}else{
    param <- BiocParallel::SerialParam(progressbar = TRUE)
}

raw_counts <- as.matrix(count.table)
if(opt$norm == TRUE){

    save_folder_norm <- file.path(opt$outdir, "normalization")
    unlink(save_folder_norm, recursive=TRUE, force=TRUE)
    dir.create(save_folder_norm, recursive=TRUE, showWarnings=FALSE)

    Sample_IDs      <- colnames(raw_counts)
    Sample_Condition <- gsub("_R[0-9]+.*", "", Sample_IDs)
    Sample_Replicate <- gsub(".*_R", "R", Sample_IDs)

    if (anyDuplicated(Sample_IDs)) stop("Sample IDs are not unique after suffix parsing.")

    design_df <- data.frame(
      Sample_ID        = Sample_IDs,
      Sample_Condition = Sample_Condition,
      Sample_Replicate = Sample_Replicate,
      stringsAsFactors = FALSE,
      row.names        = Sample_IDs
    )

    x <- raw_counts

    # guard zero-sum columns in depth scaling
    cs <- colSums(x)
    cs_safe <- pmax(cs, 1)
    x_depth_p <- sweep(x, 2, cs_safe, "/") * 1e6

    # robust 0.1% trimming
    n <- nrow(x_depth_p)
    start <- max(1, floor(n*0.001) + 1)
    end   <- max(start, ceiling(n*0.999))
    x <- x[order(rowMeans(x_depth_p), decreasing=TRUE), , drop=FALSE]
    x <- x[start:end, , drop=FALSE]

    Result <- RunNorm(
      mat_path         = x,
      design_path      = design_df,
      fix_reference    = "random",
      row_name_index   = NULL,
      n_pop            = opt$n_pop,
      n_pop_reference  = opt$n_pop,
      saving_path      = save_folder_norm,
      sigma_times      = opt$sigma_times,
      dist_family      = "Skew.normal",
      Norm_plot        = TRUE,
      Save_results     = TRUE,
      return_full      = FALSE,
      BiocParam        = param,
      cpus             = cpus
    )

    omni_scaling <- Result$scaling_factors
    omni_norm_mat <- Result$norm_mat  # matrix

    # also write a tidy OmniNorm-normalized matrix with gene IDs
    # Coerce to base matrix → data.frame so cbind.Matrix is not used
    omni_norm_df <- data.frame(
      gene_id = rownames(omni_norm_mat),
      as.matrix(omni_norm_mat),
      check.names = FALSE
    )

    write.table(
      omni_norm_df,
      file = file.path(save_folder_norm, "OmniNorm_normalized_counts.txt"),
      sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
    )
}

################################################
################################################
## RUN DESEQ2                                 ##
################################################
################################################
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive=TRUE)
setwd(opt$outdir)

samples.vec     <- colnames(raw_counts)
name_components <- strsplit(samples.vec, "_")
n_components    <- length(name_components[[1]])
decompose       <- n_components != 1 && all(sapply(name_components, length) == n_components)
coldata         <- data.frame(samples.vec, sample=samples.vec, row.names=1)

if (decompose) {
  groupings  <- as.data.frame(lapply(1:n_components, function(i) sapply(name_components, "[[", i)))
  n_distinct <- sapply(groupings, function(grp) length(unique(grp)))
  groupings  <- groupings[n_distinct != 1 & n_distinct != length(samples.vec)]
  if (ncol(groupings) != 0) {
    names(groupings) <- paste0("Group", 1:ncol(groupings))
    coldata <- cbind(coldata, groupings)
  } else {
    decompose <- FALSE
  }
}

DDSFile <- paste0(opt$outprefix, ".dds.RData")

dds <- DESeqDataSetFromMatrix(countData = round(raw_counts[, samples.vec, drop=FALSE]),
                              colData   = coldata,
                              design    = ~ 1)

dds <- estimateSizeFactors(dds)

if (isTRUE(opt$norm)) {
  # Use OmniNorm's scaling: DESeq2 offsets should be inverse of scaling
  stopifnot(exists("omni_scaling"))
  sf <- 1 / omni_scaling[ names(sizeFactors(dds)), "scaling" ]
  sizeFactors(dds) <- sf
}

# Apply the single set of sample size factors to every gene region.
# sizeFactors are DESeq2 denominators; the file scaling factors are their inverse.
sf_used <- 1 / sizeFactors(dds)
stopifnot(identical(names(sf_used), colnames(raw_counts)),
          all(is.finite(sf_used)), all(sf_used > 0))

write_region <- function(mat, label) {
  if (!setequal(colnames(mat), names(sf_used)))
    stop("Sample mismatch in ", label, " count matrix")
  if (anyDuplicated(rownames(mat))) stop("Duplicate gene IDs in ", label)
  mat <- as.matrix(mat[, names(sf_used), drop=FALSE])
  storage.mode(mat) <- "double"
  if (anyNA(mat) || any(mat < 0)) stop("Invalid raw counts in ", label)
  normalized <- sweep(mat, 2L, sf_used, `*`)
  write.table(data.frame(gene_id=rownames(mat), mat, check.names=FALSE),
              file=paste0(label,"_READS_raw.txt"), sep="\t", row.names=FALSE, quote=FALSE)
  write.table(data.frame(gene_id=rownames(normalized), normalized, check.names=FALSE),
              file=paste0(label,"_READS_norm.txt"), sep="\t", row.names=FALSE, quote=FALSE)
  invisible(normalized)
}
write_region(dat_all, "ALL")
write_region(dat_int, "INTRON")
write_region(dat_ex, "EXON")

# Keep backward-compatible downstream filename for the ALL matrix.
file.copy("ALL_READS_norm.txt", "Final_norm_count_mat.txt", overwrite=TRUE)
write.table(data.frame(Sample_ID=names(sf_used), scaling=unname(sf_used)),
            file="scaling_dat.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
dir.create("size_factors", showWarnings=FALSE)
for (sample_id in names(sf_used))
  write(sf_used[[sample_id]], file=file.path("size_factors",paste0(sample_id,".txt")))
normFactors <- sizeFactors(dds)
save(normFactors, file=file.path("size_factors", paste0(opt$outprefix, ".size_factors.RData")))

if (min(dim(raw_counts)) <= 1) {
  save(dds, file=DDSFile); saveRDS(dds, file=sub("\\.dds\\.RData$", ".rds", DDSFile))
  warning("Not enough samples or genes in counts file for PCA.", call.=FALSE)
  quit(save="no", status=0, runLast=FALSE)
}

vst_name <- if (!isTRUE(opt$vst)) "rlog" else "vst"
rld <- if (vst_name == "rlog") rlog(dds) else varianceStabilizingTransformation(dds)
assay(dds, vst_name) <- assay(rld)

save(dds, file=DDSFile)
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

PlotFile <- paste0(opt$outprefix, ".plots.pdf")
pdf(file=PlotFile, onefile=TRUE, width=7, height=7)
on.exit({ if (grDevices::dev.cur() > 1) grDevices::dev.off() }, add=TRUE)
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
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
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
library(ComplexHeatmap)

sampleDists      <- dist(t(assay(dds, vst_name)))
sampleDistMatrix <- as.matrix(sampleDists)
colors           <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

ht <- Heatmap(
    sampleDistMatrix,
    name = "Distance",
    col = colors,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_title = paste("Euclidean distance between", vst_name, "of samples"),
    show_row_names = TRUE,
    show_column_names = TRUE
)

draw(ht)

# Close the PDF device opened above
if (grDevices::dev.cur() > 1) grDevices::dev.off()

## WRITE SAMPLE DISTANCES TO FILE
cn_add = c("sample", colnames(sampleDistMatrix))
write.table(cbind(sample = rownames(sampleDistMatrix), sampleDistMatrix),file=paste(opt$outprefix, ".sample.dists.txt", sep=""),
            row.names=FALSE, sep="\t", quote=FALSE)


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

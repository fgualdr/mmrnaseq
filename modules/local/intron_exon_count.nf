process INTRON_EXON_COUNT {
    tag "$meta.id"
    label 'process_medium'

    container 'docker://fgualdr/envr_norm_gr'

    input:
    tuple val(meta), path(bam), path(bai)
    path txdb_sqlite
    path metadata
    path rmsker

    output:
    path '*.gene.rds', emit: count_rds
    path '*.rmsk.rds', emit: rmsk_rds

    script:
    
    def strandedness = 'non_specific'
    if (meta.strandedness == 'forward') {
        strandedness = 'forward'
    } else if (meta.strandedness == 'reverse') {
        strandedness = 'reverse'
    }
    // define single/paired-end
    def paired_end = meta.single_end ? 'single' :  'paired'
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    intron_exon_count.r \\
            -q $txdb_sqlite \\
            -r $rmsker \\
            -o './' \\
            -m $metadata \\
            -b $bam \\
            -i $prefix \\
            -s $strandedness \\
            -p $paired_end \\
            -c $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-GenomicFeatures: \$(Rscript -e "library(GenomicFeatures); cat(as.character(packageVersion('GenomicFeatures')))")
        bioconductor-GenomicAlignments: \$(Rscript -e "library(GenomicAlignments); cat(as.character(packageVersion('GenomicAlignments')))")
        bioconductor-GenomicRanges: \$(Rscript -e "library(GenomicRanges); cat(as.character(packageVersion('GenomicRanges')))")
    END_VERSIONS
    """
}

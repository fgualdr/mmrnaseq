process PREP_GTF {
    tag "$gtf"
    label 'process_medium'

    container 'docker://fgualdr/envr_norm_gr'

    input:
    path gtf
    path chrom_size

    output:
    path 'TxDB.sqlite'       , emit: txdb_sqlite
    path 'metadata.txt'      , emit: metadata
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    """
    gtf_prep.r \\
            -g ${gtf} \\
            -o './' \\
            -s ${chrom_size} \\
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

// Reference-wide RepeatMasker GRanges + flanks. Cached once per matching inputs.
process PREP_REPEAT_ANNOTATION {
    label 'process_high'
    container 'docker://fgualdr/envr_norm_gr'

    input:
    path txdb_sqlite
    path rmsker

    output:
    path 'repeat_annotations.rds', emit: annotation
    path 'versions.yml', emit: versions

    script:
    """
    prepare_repeat_annotation.r \\
        --txdb_sqlite $txdb_sqlite \\
        --rmsker $rmsker \\
        --output repeat_annotations.rds \\
        --cpus $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript --version 2>&1 | head -n1 | sed 's/"/\\"/g')
    END_VERSIONS
    """
}

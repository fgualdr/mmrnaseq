process MM_BAYES_RESOLVER {
    tag "$meta.id"
    label 'process_high_memory'

    // Reuse the pipeline's established MultiQC/Python biocontainer. It already
    // runs correctly under Nextflow 25.10 task-metric collection and provides
    // Python + NumPy, while avoiding a separate C++ build/runtime image.
    conda "bioconda::multiqc=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(candidates), path(unique_fragments), path(unique_anchors)
    path resolver_script
    path chrom_sizes

    output:
    tuple val(meta), path("*.posterior_alignment.tsv.gz"), emit: posterior
    tuple val(meta), path("*.resolution.tsv"), emit: resolution
    tuple val(meta), path("*.resolver_qc.tsv"), emit: qc
    tuple val(meta), path("*.density_model.tsv"), emit: density_model
    tuple val(meta), path("*.density_cv.tsv"), emit: density_cv
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    def fdr = params.mm_resolution_fdr ?: 0.05
    """
    python3 $resolver_script \
        --candidates $candidates \
        --unique-fragments $unique_fragments \
        --anchors $unique_anchors \
        --chrom-sizes $chrom_sizes \
        --resolution-fdr $fdr \
        --output-prefix ${prefix}

    gzip -n -f ${prefix}.posterior_alignment.tsv

    pyver=\$(python3 -c 'import platform; print(platform.python_version())')
    npver=\$(python3 -c 'import numpy; print(numpy.__version__)')
    printf '"%s":\n    resolver: "STAR-candidate Bayesian resolver Python v1 + PATCH05 density-CV diagnostics"\n    python: "%s"\n    numpy: "%s"\n' \
        '${task.process}' "\$pyver" "\$npver" > versions.yml
    """
}

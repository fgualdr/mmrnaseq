process MM_QC_MERGE {
    tag "$meta.id"
    label 'process_low'

    // Reuse the already validated Python/NumPy runtime used by the resolver.
    // Unlike python:* slim images this container has worked with Nextflow 25.10
    // task-metric collection on the target cluster.
    conda "bioconda::multiqc=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(contract_qc), path(resolver_qc), path(realization_qc), path(density_model), path(density_cv)

    output:
    tuple val(meta), path('*.mm.qc.tsv'), emit: qc
    tuple val(meta), path('*.mmresolver_mqc.tsv'), emit: multiqc
    tuple val(meta), path('*.mmdensity_mqc.tsv'), emit: density_multiqc
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    python3 ${projectDir}/bin/mm_merge_qc.py \
        --sample '${meta.id}' \
        --contract $contract_qc \
        --resolver $resolver_qc \
        --realization $realization_qc \
        --density $density_model \
        --density-cv $density_cv \
        --output-prefix ${prefix}

    pyver=\$(python3 -c 'import platform; print(platform.python_version())')
    printf '"%s":\n    python: "%s"\n    qc_contract: "PATCH05 cross-stage fragment accounting + density-k diagnostics"\n' \
        '${task.process}' "\$pyver" > versions.yml
    """
}

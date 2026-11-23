process MM_LOCUS_EM {
    tag "$meta.id"
    label 'process_high_memory'

    // Supply a tested image containing Python 3.11, NumPy and samtools for
    // Singularity/Apptainer; alternatively enable the pinned Conda environment.
    conda "${projectDir}/envs/mm_graph_em.yml"
    container "${params.mm_experimental_container ?: ''}"

    input:
    tuple val(meta), path(candidates), path(scores), path(unique_anchors)
    path loci

    output:
    tuple val(meta), path('*.em_locus_abundance.tsv'), emit: abundance
    tuple val(meta), path('*.em_candidate_posterior.tsv'), emit: posterior
    tuple val(meta), path('*.em_qc.tsv'), emit: qc
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    python3 ${projectDir}/bin/mm_locus_em.py \
        --candidates $candidates \
        --scores $scores \
        --loci $loci \
        --output-prefix ${prefix} \
        --max-edges ${params.mm_em_max_edges}

    printf '"%s":\n    inference: "experimental_conditional_EM_not_used_for_FINAL_BAM"\n    numpy: "%s"\n' \
        '${task.process}' "\$(python3 -c 'import numpy; print(numpy.__version__)')" > versions.yml
    """
}

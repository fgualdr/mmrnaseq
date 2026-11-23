process MM_QUALITY_SCORE {
    tag "$meta.id"
    label 'process_high'

    // Supply a tested image containing Python 3.11, NumPy and samtools for
    // Singularity/Apptainer; alternatively enable the pinned Conda environment.
    conda "${projectDir}/envs/mm_graph_em.yml"
    container "${params.mm_experimental_container ?: ''}"

    input:
    tuple val(meta), path(candidates), path(namesorted_bam)
    path fasta

    output:
    tuple val(meta), path('*.quality_scores.tsv'), emit: scores
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    samtools faidx --fai-idx ${prefix}.fai $fasta
    samtools view -h $namesorted_bam | \
        python3 ${projectDir}/bin/mm_quality_scores.py \
            --candidates $candidates --sam - --fasta $fasta --fai ${prefix}.fai \
            --output ${prefix}.quality_scores.tsv

    printf '"%s":\n    samtools: "%s"\n    quality_scoring: "experimental_phred_indel_pseudolikelihood"\n' \
        '${task.process}' "\$(samtools --version | head -1 | awk '{print \$2}')" > versions.yml
    """
}

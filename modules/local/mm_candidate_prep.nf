process MM_CANDIDATE_PREP {
    tag "$meta.id"
    label 'process_samtools_sort'

    conda (params.enable_conda ? "bioconda::samtools=1.17 coreutils=9.3 gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-9d3a458f6420e5712103ae2af82c94d26d63f059:60b54b43045e8cf39ba307fd683c69d4c57240ce-0' :
        'quay.io/biocontainers/mulled-v2-9d3a458f6420e5712103ae2af82c94d26d63f059:60b54b43045e8cf39ba307fd683c69d4c57240ce-0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.mm.candidates.tsv.gz"), emit: candidates
    tuple val(meta), path("*.mm.unique_fragments.tsv.gz"), emit: unique_fragments
    tuple val(meta), path("*.mm.unique_anchors.tsv.gz"), emit: unique_anchors
    tuple val(meta), path("*.mm.excluded.tsv.gz"), emit: excluded
    tuple val(meta), path("*.mm.contract_qc.tsv"), emit: qc
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? task.memory.mega : false
    def sort_memory = avail_mem ? Math.max(256, ((avail_mem - 1000) / task.cpus).intValue()) : 2000
    def singleEnd = meta.single_end ? 1 : 0
    def strand = meta.strandedness ?: 'unstranded'
    """
    samtools sort \
        -@ ${task.cpus} \
        -m ${sort_memory}M \
        -n \
        -T ${prefix}.mm.nsort \
        -o ${prefix}.mm.nsort.bam \
        $bam

    samtools view ${prefix}.mm.nsort.bam | \
        awk \
            -v single_end=${singleEnd} \
            -v strandedness='${strand}' \
            -v candidate_out=${prefix}.mm.candidates.tsv \
            -v unique_out=${prefix}.mm.unique_fragments.tsv \
            -v anchor_out=${prefix}.mm.unique_anchors.tsv \
            -v excluded_out=${prefix}.mm.excluded.tsv \
            -v qc_out=${prefix}.mm.contract_qc.tsv \
            -f ${projectDir}/bin/mm_build_candidates.awk

    gzip -n -f ${prefix}.mm.candidates.tsv
    gzip -n -f ${prefix}.mm.unique_fragments.tsv
    gzip -n -f ${prefix}.mm.unique_anchors.tsv
    gzip -n -f ${prefix}.mm.excluded.tsv

    printf '"%s":\n    samtools: "%s"\n    awk: "%s"\n' \
        '${task.process}' \
        "\$(samtools --version | awk 'NR==1{print \$2}')" \
        "\$(awk -W version 2>&1 | head -1 | sed 's/"/\\"/g')" \
        > versions.yml
    """
}

process MM_BAM_REALIZE {
    tag "$meta.id"
    label 'process_samtools_sort'

    conda (params.enable_conda ? "bioconda::samtools=1.17 coreutils=9.3 gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-9d3a458f6420e5712103ae2af82c94d26d63f059:60b54b43045e8cf39ba307fd683c69d4c57240ce-0' :
        'quay.io/biocontainers/mulled-v2-9d3a458f6420e5712103ae2af82c94d26d63f059:60b54b43045e8cf39ba307fd683c69d4c57240ce-0' }"

    input:
    tuple val(meta), path(bam), path(resolution), path(excluded) // bam is query-name sorted by MM_CANDIDATE_PREP

    output:
    tuple val(meta), path("*.mm.full.unsorted.bam"),         optional: true, emit: full_bam
    tuple val(meta), path("*.mm.quantitative.unsorted.bam"), emit: quantitative_bam
    tuple val(meta), path("*.mm.highconf.unsorted.bam"),     optional: true, emit: highconf_bam
    tuple val(meta), path("*.mm.realization_qc.tsv"),         emit: qc
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    def fdr = (params.mm_resolution_fdr ?: 0.05) as double
    def highconfMin = 1.0d - fdr
    // Preserve the AWK accounting in every run, but materialize auxiliary BAM tiers
    // only when requested. /dev/null avoids writing huge diagnostic SAM/BAM files.
    def fullRequested = params.save_align_intermeds
    def highconfRequested = params.save_align_intermeds && !params.skip_em
    def fullSam = fullRequested ? "${prefix}.mm.full.sam" : '/dev/null'
    def highconfSam = highconfRequested ? "${prefix}.mm.highconf.sam" : '/dev/null'
    def fullConvert = fullRequested ? "samtools view -bS ${prefix}.mm.full.sam > ${prefix}.mm.full.unsorted.bam" : ''
    def highconfConvert = highconfRequested ? "samtools view -bS ${prefix}.mm.highconf.sam > ${prefix}.mm.highconf.unsorted.bam" : '' 
    """
    {
        cat $resolution
        gzip -dc $excluded | awk 'BEGIN{FS=OFS="\\t"} NR>1{print \$1,"EXCLUDED",".",0,0,0,".",0,1,0,0,0,0,0,0,0}'
    } > ${prefix}.mm.decisions.tsv

    samtools view -h $bam | \
        awk \
            -v full_out=${fullSam} \
            -v quantitative_out=${prefix}.mm.quantitative.sam \
            -v highconf_out=${highconfSam} \
            -v highconf_min=${highconfMin} \
            -v qc_out=${prefix}.mm.realization_qc.tsv \
            -f ${projectDir}/bin/mm_realize_bam.awk \
            ${prefix}.mm.decisions.tsv -

    samtools view -bS ${prefix}.mm.quantitative.sam > ${prefix}.mm.quantitative.unsorted.bam
    ${fullConvert}
    ${highconfConvert}

    printf '"%s":\n    samtools: "%s"\n    awk: "%s"\n    highconf_min_pmax: "%s"\n' \
        '${task.process}' \
        "\$(samtools --version | awk 'NR==1{print \$2}')" \
        "\$(awk -W version 2>&1 | head -1 | sed 's/"/\\"/g')" \
        '${highconfMin}' \
        > versions.yml
    """
}

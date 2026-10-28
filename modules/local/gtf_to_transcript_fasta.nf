process GTF_TO_TRANSCRIPT_FASTA {
    tag "$gtf"
    label 'process_low'

    conda "bioconda::gffread=0.12.7"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffread:0.12.7--hdcf5f25_4' :
        'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4' }"

    input:
    path fasta
    path gtf

    output:
    path "*.transcripts.fa", emit: fasta
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: fasta.baseName
    """
    gffread \
        $gtf \
        -g $fasta \
        -w ${prefix}.transcripts.fa

    printf '%s\n' \
      '"${task.process}":' \
      '    gffread: "'"\$(gffread --version 2>&1 | head -1 | sed 's/"/\\"/g')"'"' \
      > versions.yml
    """
}

/*
 * Download SRA files from NCBI
 */
process SRX_DOWNLOAD {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::sra-tools>=3.0.0 bioconda::parallel-fastq-dump" : null)
    container 'docker://fgualdr/parallel_fastq_dump:latest'

    input:
    val(meta)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: fastq
    path "versions.yml", emit: versions

    script:

    def prefix           = task.ext.prefix ?: "${meta.id}"
    def srr             = meta.run_accession

    """
    parallel-fastq-dump \\
            --sra-id $srr \\
            --threads $task.cpus \\
            --outdir ./ \\
            --split-files \\
            --gzip \\
            --tmpdir ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parallel_fastq_dump: \$(parallel-fastq-dump --version | grep 'parallel-fastq-dump : ' | sed -e "s/parallel-fastq-dump : //g")
        fastq-dump: \$(parallel-fastq-dump --version | grep 'fastq-dump' | grep -v 'parallel-fastq-dump : ' | sed -e "s/fastq-dump : //g")
    END_VERSIONS
    """
}
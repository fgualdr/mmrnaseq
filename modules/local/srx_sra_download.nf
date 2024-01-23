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
    // add an if else in case "*._1.fastq.gz" and "*._2.fastq.gz" are not present
    // in case both are present 
    // in case only _1 is present therefore _2 is null and single_end must be turned on
    if(path("*._1.fastq.gz").isFile() && path("*._2.fastq.gz").isFile()){
        tuple val(meta), path("*._1.fastq.gz"), path("*._2.fastq.gz"), emit: fastq
    } else if(path("*._1.fastq.gz").isFile() && !path("*._2.fastq.gz").isFile()){
        tuple val(meta), path("*._1.fastq.gz"), null, emit: fastq
    } else {
        tuple val(meta), null, path("*._2.fastq.gz"), emit: fastq
    }
    path "versions.yml"           , emit: versions

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
        parallel-fastq-dump: \$(echo \$(parallel-fastq-dump --version 2>&1) | sed 's/^.*parallel-fastq-dump //; s/Using.*\$//')
    END_VERSIONS
    """
}


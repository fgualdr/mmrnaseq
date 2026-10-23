process FIX_SRA_FASTQ {

    label 'process_high'

    conda "bioconda::bbmap=39.01"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap:39.01--h5c4e2a8_0' :
        'quay.io/biocontainers/bbmap:39.01--h5c4e2a8_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def r1 = fastq.find { it.name =~ /_1\.fastq\.gz$/ }
    def r2 = fastq.find { it.name =~ /_2\.fastq\.gz$/ }

    if (!r1 || !r2) {
        throw new IllegalArgumentException("Paired-end FASTQ files not found in: ${fastq*.name}")
    }

    def prefix = meta.id

    """
    set -euo pipefail

    echo "[INFO] Repairing FASTQ pairing for ${meta.id}"
    repair.sh \\
        in1=${r1} \\
        in2=${r2} \\
        out1=${prefix}_repaired_1.fastq.gz \\
        out2=${prefix}_repaired_2.fastq.gz \\
        outs=${prefix}_unpaired.fastq.gz \\
        overwrite=true

    echo "[INFO] Writing version file..."
    echo \\
    'FIX_SRA_FASTQ:
        bbmap: $(repair.sh version | head -n 1)' > versions.yml
    """
}
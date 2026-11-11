process GTF2BED {
    tag "$gtf"
    label 'process_low'

    conda "conda-forge::perl=5.26.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'quay.io/biocontainers/perl:5.26.2' }"

    input:
    path gtf
    path fasta

    output:
    path '*.bed'       , emit: bed
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = gtf.baseName
    """
    set -euo pipefail

    # The compatibility check/filter is implemented inside gtf2bed itself so
    # this process only relies on Perl, which is guaranteed by the container.
    # Do not add coreutils-only helpers (comm/sort/cut/wc/head/tr) here.
    gtf2bed \
        --fasta $fasta \
        $gtf \
        > ${prefix}.bed

    test -s ${prefix}.bed

    printf '%s\n' \
        '"${task.process}":' \
        "    perl: \$(perl -e 'print \$]')" \
        > versions.yml
    """
}

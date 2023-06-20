/*
 * Filter BAM file
 */
process BAM_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0':
        'quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.filter2.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:

    def prefix           = task.ext.prefix ?: "${meta.id}"
    def filter_params    = meta.single_end ? '-F 0x004' : '-F 0x004 -F 0x0008 -f 0x001 -f 0x002' // proper pair selection!!

    """
    # we set the maximum allowed insert_size to a fixed 50kb
    max_frag=50000
    
    # General filtering of the bam
    samtools view \\
        $filter_params \\
        -b $bam > ${prefix}.filter1.bam

    # Filter those pairs that have insert size compatible with the fragment length:
    # For the RNAseq this is set to 5 times the mimum gettin 90% of the data ~ 10kb *5
    # This will be possible to set it accordingly to the longer interon or the tapestation data
    
    samtools view -h ${prefix}.filter1.bam | \\
    awk -v var="\$max_frag" '{if(substr(\$0,1,1)=="@" || ((\$9>=0?\$9:-\$9)<=var)) print \$0}' | \\
    samtools view -b > ${prefix}.filter2.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

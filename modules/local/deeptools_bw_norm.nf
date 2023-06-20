process DEEPTOOLS_BIGWIG_NORM {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::deeptools=3.5.1 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:62d1ebe2d3a2a9d1a7ad31e0b902983fa7c25fa7-0':
        'quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:62d1ebe2d3a2a9d1a7ad31e0b902983fa7c25fa7-0' }"

    input:
    tuple val(meta), path(bam), path(bai), val(scaling)

    output:
    path "*.unstranded.norm.bw" , optional:true, emit: unstranded_bw
    path "*.fwd.norm.bw"        , optional:true, emit: fw_bw
    path "*.rev.norm.bw"        , optional:true, emit: rev_bw
    path "versions.yml"         , emit: versions
 
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def pe     = meta.single_end ? '' : '-pc'
    def extend = (meta.single_end && params.fragment_size > 0) ? "-fs ${params.fragment_size}" : ''

    def strandedness = 'non_stranded'
    if (meta.strandedness == 'forward') {
        strandedness = 'forward'
    } else if (meta.strandedness == 'reverse') {
        strandedness = 'reverse'
    }

    if(strandedness == 'non_stranded'){
        """
        echo $scaling
        echo $prefix
        
        bamCoverage \\
                --numberOfProcessors $task.cpus \\
                --binSize 1 \\
                --scaleFactor $scaling \\
                --bam $bam \\
                -o ${prefix}.unstranded.norm.bw

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            deeptools: \$(computeMatrix --version | sed -e "s/computeMatrix //g")
        END_VERSIONS
        """
    } else {
        """
        echo $scaling
        echo $prefix

        # include reads that are 2nd in a pair (128);
        # exclude reads that are mapped to the reverse strand (16)
        samtools view -b -f 128 -F 16 $bam > ${prefix}.fwd1.bam

        # exclude reads that are mapped to the reverse strand (16) and
        # first in a pair (64): 64 + 16 = 80
        samtools view -b -f 80 $bam > ${prefix}.fwd2.bam

        # combine the temporary files
        samtools merge -f ${prefix}.fwd.bam ${prefix}.fwd1.bam ${prefix}.fwd2.bam
        rm -rf ${prefix}.fwd1.bam ${prefix}.fwd2.bam

        # index the filtered BAM file
        samtools index ${prefix}.fwd.bam

        # run bamCoverage
        bamCoverage \\
            --bam ${prefix}.fwd.bam \\
            --numberOfProcessors $task.cpus \\
            --binSize 1 \\
            --scaleFactor $scaling \\
            -o ${prefix}.fwd.norm.bw
        
        rm -rf ${prefix}.fwd.bam ${prefix}.fwd.bam.bai

        # include reads that map to the reverse strand (128)
        # and are second in a pair (16): 128 + 16 = 144
        samtools view -b -f 144 $bam > ${prefix}.rev1.bam

        # include reads that are first in a pair (64), but
        # exclude those ones that map to the reverse strand (16)
        samtools view -b -f 64 -F 16 $bam > ${prefix}.rev2.bam

        # merge the temporary files
        samtools merge -f ${prefix}.rev.bam ${prefix}.rev1.bam ${prefix}.rev2.bam
        rm -rf ${prefix}.rev1.bam ${prefix}.rev2.bam

        # index the merged, filtered BAM file
        samtools index ${prefix}.rev.bam

        # run bamCoverage
        bamCoverage \\
            --bam ${prefix}.rev.bam \\
            --numberOfProcessors $task.cpus \\
            --binSize 1 \\
            --scaleFactor $scaling \\
            -o ${prefix}.rev.norm.bw

        rm -rf ${prefix}.rev.bam ${prefix}.rev.bam.bai

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            deeptools: \$(computeMatrix --version | sed -e "s/computeMatrix //g")
        END_VERSIONS
        """
    }
        
}

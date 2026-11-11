/*
 * Basic BAM filtering before multimapper resolution.
 *
 * RNA-seq note:
 * We intentionally do NOT filter paired-end RNA alignments by |TLEN|. TLEN is the
 * genomic outer span of a pair and therefore includes introns; a fixed 50 kb TLEN
 * threshold removes biologically valid long-intron fragments. Spurious splice
 * junctions are instead controlled in the STAR alignment step (BySJout and splice
 * junction filters; see conf/modules.config).
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
    path "versions.yml"                    , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Keep all mapped alignments for SE. For PE, require a mapped, properly paired
    // fragment with both mates mapped. Secondary alignments are deliberately retained
    // because they are candidate loci for the multimapper resolver.
    def filter_params = meta.single_end ? '-F 0x004' : '-f 0x003 -F 0x00C'

    """
    samtools view \\
        $filter_params \\
        -b $bam > ${prefix}.filter2.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

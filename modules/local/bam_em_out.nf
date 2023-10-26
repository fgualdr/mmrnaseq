/*
 * Filter BAM file
 */
process BAM_EM_OUT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0 bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0':
        'quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0' }"

    input:
    tuple val(meta), path(bam_label), path(final_read_target)
    path chrom_sizes

    output:
    tuple val(meta), path("*.filtered.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:
    def prefix           = task.ext.prefix ?: "${meta.id}"

    """
    # We save the header:
    samtools view -h $bam_label > samp.sam

    # Select the reads in the samp.sam havin the read name in the final_read_target or NH:i:1
    awk 'BEGIN{OFS="\\t"} FNR==NR{ids[\$1];next} {if(\$1~/^@/){print}else{if(\$1 in ids || \$0~/NH:i:1/){print}}}' $final_read_target samp.sam > filter.sam
    # awk 'BEGIN{OFS="\\t"} FNR==NR{ids[\$1];next} {if(\$1~/^@/){print}else{if(\$1 in ids){print}}}' $final_read_target samp.sam > filter.sam
    samtools view -bS filter.sam > ${prefix}.filtered.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

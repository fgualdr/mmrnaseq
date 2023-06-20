/*
 * SQUIRE_COUNT
 */
process SQUIRE_COUNT {
    tag "$meta.id"
    label 'process_medium'

    container 'docker://lizatym/squire'

    input:
    tuple val(meta), path(bam), path(bai)
    path ch_gtf
    path clean_folder
    val ch_read_length
 
    output:
    path("*_abund.txt") , emit: squire_abund_count
    path("*refGenecounts.txt") , emit: squire_refgenecount
    path("*subFcounts.txt") , emit: squire_subf_count
    path("*TEcounts.txt") , emit: squire_te_count
    path  "versions.yml" , emit: versions
    
    script:
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def strandedness = 0
    if (meta.strandedness == 'forward') {
        strandedness = 1
    } else if (meta.strandedness == 'reverse') {
        strandedness = 2
    }
    def iter = params.squire_iter ? params.squire_iter : 1000

    """   
    # Make the folder 
    squire Count \\
        -c $clean_folder \\
        -m ./ \\
        -r $ch_read_length \\
        -s $strandedness \\
        -p $task.cpus \\
        -e $iter \\
        -n $prefix \\
        -f './' \\
        -t './' \\
        -v True \\
        -o './'
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        squire: \$(echo 'SQuIRE_v0.9.9.9a-beta')
    END_VERSIONS
    """
    
}




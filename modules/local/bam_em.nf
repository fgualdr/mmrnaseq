/*
 * Filter BAM file
 */
process BAM_EM {
    tag "$meta.id"
    label 'process_high_memory'

    conda (params.enable_conda ? "" : null)
    container 'docker://fgualdr/empy'

    input:
    tuple val(meta), path(labeled_bedpe)

    output:
    tuple val(meta), path("*Final.bedpe"), emit: final_bedpe
    tuple val(meta), path("*posterior_target_probabilities.txt"), emit: posterior_txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    // added to control memory usage on HPC
    def prefix = task.ext.prefix ?: "${meta.id}"    
    def eps = params.em_eps ? params.em_eps : 1e-6
    def iter = params.em_iter ? params.em_iter : 1000
    def processors = task.cpus
    """
 
    em_algorithm_bedpe_sm_bis.py \\
        -i $labeled_bedpe \\
        -o ./ \\
        -m $iter \\
        -c $eps


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

/*
 * SQUIRE_CLEAN
 */
process SQUIRE_CLEAN {

    label 'process_medium'

    container 'docker://lizatym/squire'

    input:
    path ch_gtf
    path ch_rmsker
 
    output:
    path "squire_clean" , emit: clean_folder
    path  "versions.yml"                      , emit: versions

    script:

    """   
    # Make the folder 
    mkdir squire_clean

    # Unzip ch_rmsker if it ends with .gz rename the file as basename of ch_rmsker without .gz give it to a variable ch_rmsker_fixed
    if [[ $ch_rmsker == *.gz ]]; then
        gunzip -c $ch_rmsker > \$(basename $ch_rmsker .gz)
        ch_rmsker_fixed=\$(basename $ch_rmsker .gz)
    else
        ch_rmsker_fixed=$ch_rmsker
    fi
    squire Clean \\
        -r \$ch_rmsker_fixed \\
        -o squire_clean
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        squire: \$(echo 'SQuIRE_v0.9.9.9a-beta')
    END_VERSIONS
    """
    
}



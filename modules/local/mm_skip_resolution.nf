process MM_SKIP_RESOLUTION {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::gzip=1.13 coreutils=9.3 gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-9d3a458f6420e5712103ae2af82c94d26d63f059:60b54b43045e8cf39ba307fd683c69d4c57240ce-0' :
        'quay.io/biocontainers/mulled-v2-9d3a458f6420e5712103ae2af82c94d26d63f059:60b54b43045e8cf39ba307fd683c69d4c57240ce-0' }"

    input:
    tuple val(meta), path(candidates)

    output:
    tuple val(meta), path("*.mm.skip_resolution.tsv"), emit: resolution
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    gzip -dc $candidates | \
        awk 'BEGIN{FS=OFS="\\t"}
        NR==1 {
            for (i=1; i<=NF; i++) h[\$i]=i
            print "fragment_id","status","selected_candidate_id","selected_hi","nh","candidate_count","orientation_group","pmax","posterior_error","entropy","best_as","second_as","as_delta","tied_best","normalized_as","absolute_score_pass"
            next
        }
        function emit_group(    second,delta,entropy) {
            if (frag=="") return
            second=(tied>1 ? best_as : (have_second ? second_as : best_as))
            delta=best_as-second
            entropy=(n>1 ? log(n) : 0)
            print frag,"UNRESOLVED",".",0,nh,n,ori,0,1,entropy,best_as,second,delta,tied,best_norm,0
        }
        {
            f=\$(h["fragment_id"])
            if (frag!="" && f!=frag) emit_group()
            if (f!=frag) {
                frag=f
                nh=\$(h["nh"])+0
                ori=\$(h["orientation_group"])
                n=0
                best_as=-1e300
                second_as=-1e300
                have_second=0
                tied=0
                best_norm=0
            }
            n++
            asv=\$(h["as"])+0
            norm=\$(h["normalized_as"])+0
            if (asv>best_as) {
                if (best_as>-1e299) { second_as=best_as; have_second=1 }
                best_as=asv
                best_norm=norm
                tied=1
            } else if (asv==best_as) {
                tied++
            } else if (!have_second || asv>second_as) {
                second_as=asv
                have_second=1
            }
        }
        END { emit_group() }' \
        > ${prefix}.mm.skip_resolution.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -W version 2>&1 | head -1 | sed 's/"/\\"/g')
        mode: "skip_bayesian_resolution_all_multimappers_unresolved"
    END_VERSIONS
    """
}

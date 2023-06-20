process NORMALIZE_COUNTS_INTRON_EXON {
    label "process_medium"

    // (Bio)conda packages have intentionally not been pinned to a specific version
    // This was to avoid the pipeline failing due to package conflicts whilst creating the environment when using -profile conda
    container 'docker://fgualdr/envnorm_fix'

    input:
    path rds
    path pca_header_multiqc
    path clustering_header_multiqc

    output:
    path "normalization"        , optional:true, emit: noamlization
    path "scaling_dat.txt"      , optional:true, emit: noamlization_txt
    path "ALL_READS_norm.txt"                 , optional:true, emit:final_all_scaled_mat
    path "INTRON_READS_norm.txt"                 , optional:true, emit:final_intron_scaled_mat
    path "EXON_READS_norm.txt"                 , optional:true, emit:final_exon_scaled_mat
    path "*.pdf"                , optional:true, emit: pdf
    path "*.RData"              , optional:true, emit: rdata
    path "*pca.vals.txt"        , optional:true, emit: pca_txt
    path "*pca.vals_mqc.tsv"    , optional:true, emit: pca_multiqc
    path "*sample.dists.txt"    , optional:true, emit: dists_txt
    path "*sample.dists_mqc.tsv", optional:true, emit: dists_multiqc
    path "*.log"                , optional:true, emit: log
    path "size_factors"         , optional:true, emit: size_factors
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args  ?: ''
    def args2 = task.ext.args2 ?: ''
    def label_lower = args2.toLowerCase()
    def label_upper = args2.toUpperCase()

    def norm = params.normalize ? 'TRUE' : 'FALSE'
    def sigma_times = params.sigma_times ? params.sigma_times : '1'
    def n_pop = params.n_pop ? params.n_pop : '1'

    """
    intron_exon_count_norm.r \\
        --outdir ./ \\
        --norm $norm \\
        --sigma_times $sigma_times \\
        --n_pop $n_pop \\
        --cores $task.cpus \\
        $args

    if [ -f "R_sessionInfo.log" ]; then
        sed "s/deseq2_pca/${label_lower}_deseq2_pca/g" <$pca_header_multiqc >tmp.txt
        sed -i -e "s/DESeq2 PCA/${label_upper} DESeq2 PCA/g" tmp.txt
        cat tmp.txt *.pca.vals.txt > ${label_lower}.pca.vals_mqc.tsv

        sed "s/deseq2_clustering/${label_lower}_deseq2_clustering/g" <$clustering_header_multiqc >tmp.txt
        sed -i -e "s/DESeq2 sample/${label_upper} DESeq2 sample/g" tmp.txt
        cat tmp.txt *.sample.dists.txt > ${label_lower}.sample.dists_mqc.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}

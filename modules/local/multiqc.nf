process MULTIQC {
    label 'process_medium'

    conda "bioconda::multiqc=1.35"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.35--pyhdfd78af_2' :
        'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_2' }"

    input:
    path multiqc_config
    path multiqc_custom_config
    path software_versions
    path workflow_summary
    path methods_description
    path logo
    path fail_trimming_summary
    path fail_mapping_summary
    path strand_qc
    path fragment_usability_summary
    path mmresolver_qc
    path mmdensity_qc
    path ('fastqc/raw/*')
    path ('fastqc/trim/*')
    path ('trim_log/*')
    path ('sortmerna/*')
    path ('star/*')
    path ('salmon/*')
    path ('samtools/stats/*')
    path ('samtools/flagstat/*')
    path ('samtools/idxstats/*')
    path ('picard/markduplicates/*')
    path ('featurecounts/*')
    path ('deseq2/aligner_pca/*')
    path ('deseq2/aligner_clustering/*')
    path ('preseq/*')
    path ('qualimap/*')
    path ('dupradar/*')
    path ('rseqc/bam_stat/*')
    path ('rseqc/infer_experiment/*')
    path ('rseqc/inner_distance/*')
    path ('rseqc/junction_annotation/*')
    path ('rseqc/junction_saturation/*')
    path ('rseqc/read_distribution/*')
    path ('rseqc/read_duplication/*')
    path ('rseqc/tin/*')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    set -euo pipefail

    # MultiQC 1.35 Biocontainer uses Python 3.14. On Linux, Python 3.14
    # defaults multiprocessing to forkserver; on NFS this can break / hang
    # Kaleido static plot export. Force the portable spawn start method and
    # put multiprocessing/Kaleido temporary state on node-local /tmp.
    mqc_tmp=\$(mktemp -d /tmp/mmrnaseq_multiqc.XXXXXX)
    trap 'rm -rf "\$mqc_tmp"' EXIT
    printf '%s\n' 'import multiprocessing as _mp' 'try:' '    _mp.set_start_method("spawn", force=True)' 'except RuntimeError:' '    pass' > "\$mqc_tmp/sitecustomize.py"
    export PYTHONPATH="\$mqc_tmp\${PYTHONPATH:+:\$PYTHONPATH}"
    export TMPDIR="\$mqc_tmp"

    multiqc -f $args $custom_config .

    # PATCH07 contract: if static export is requested, publishing empty
    # png/svg/pdf directories is a QC failure, not a successful report.
    plots_dir=\$(find . -maxdepth 1 -type d -name '*_plots' -print -quit)
    if [ -z "\$plots_dir" ]; then
        echo "ERROR: MultiQC static plot directory was not created." >&2
        exit 1
    fi
    for ext in png svg pdf; do
        n=\$(find "\$plots_dir" -type f -name "*.\$ext" | wc -l | tr -d ' ')
        if [ "\$n" -eq 0 ]; then
            echo "ERROR: MultiQC static plot export produced zero .\$ext files. Check Kaleido, TMPDIR and export timeout." >&2
            exit 1
        fi
    done

    printf '%s\n' '"${task.process}":' '    multiqc: "'"\$(multiqc --version 2>&1 | sed -e 's/multiqc, version //g' -e 's/"/\\"/g')"'"' > versions.yml
    """
}

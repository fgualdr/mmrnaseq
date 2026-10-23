process SRX_DOWNLOAD {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::sra-tools>=3.0.0 bioconda::parallel-fastq-dump bioconda::bbmap" : null)
    container 'docker://fgualdr/parallel_fastq_dump_bbmap:latest'

    input:
    val(meta)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: fastq
    path "versions.yml", emit: versions

    script:
    def srr = meta.run_accession
    def prefix = task.ext.prefix ?: meta.id
    def cpu = task.cpus

    """
    set -euo pipefail

    srr="${srr}"
    prefix="${prefix}"
    cpu="${cpu}"

    echo "[INFO] Downloading SRA accession \$srr with prefetch..."
    prefetch \$srr --max-size u
    mkdir -p ./tmp
    echo "[INFO] Converting to FASTQ with fasterq-dump..."
    fasterq-dump \$srr --split-files -e \$cpu -O ./ --temp ./tmp

    r1="${srr}_1.fastq"
    r2="${srr}_2.fastq"
    r_single="${srr}.fastq"

    if [[ -f "\$r1" && -f "\$r2" ]]; then
        echo "[INFO] Detected paired-end data for \$srr"

        reads_r1=\$(wc -l < "\$r1")
        reads_r2=\$(wc -l < "\$r2")

        echo "[gzip] the files"
        gzip "\$r1" "\$r2"
        r1="\$r1.gz"
        r2="\$r2.gz"
        
        rm -f *.fastq

        if [[ \$((reads_r1/4)) -ne \$((reads_r2/4)) ]]; then
            echo "[WARNING] Paired-end read counts mismatch: R1=\$((reads_r1/4)) reads, R2=\$((reads_r2/4)) reads"
            echo "[INFO] Running repair.sh to fix..."

            repair.sh \\
                in="\$r1" \\
                in2="\$r2" \\
                out="\${srr}_repaired_1.fastq.gz" \\
                out2="\${srr}_repaired_2.fastq.gz" \\
                outs="\${srr}_unpaired.fastq.gz" \\
                overwrite=true

            rm -f "\${srr}_1.fastq.gz" "\${srr}_2.fastq.gz"
            mv "\${srr}_repaired_1.fastq.gz" "\${srr}_1.fastq.gz"
            mv "\${srr}_repaired_2.fastq.gz" "\${srr}_2.fastq.gz"
            echo "[INFO] Repair completed."
        else
            echo "[INFO] Paired-end FASTQs already synchronized."
        fi
    elif [[ -f "\$r_single" ]]; then
        echo "[INFO] Detected single-end data for \$srr"
        reads=\$(wc -l < "\$r_single")
        if [[ \$((reads % 4)) -ne 0 ]]; then
            echo "[ERROR] Incomplete single-end FASTQ: lines=\$reads"
            rm "\$r_single"
            exit 1
        else
            echo "[INFO] FASTQ file is complete."
            echo "[INFO] GZIP the file"
            # rename as _1.fastq.gz
            mv "\$r_single" "\${srr}_1.fastq"
            r_single="\${srr}_1.fastq"
            # gzip the file
            gzip "\$r_single"
            rm -f *.fastq
        fi
    else
        echo "[[DEBUG]] Check if \$r1 is the only present so single end named _1.fastq"
        if [[ -f "\$r1" ]]; then
            echo "[INFO] Detected single-end data for \$srr"
            reads=\$(wc -l < "\$r1")
            if [[ \$((reads % 4)) -ne 0 ]]; then
                echo "[ERROR] Incomplete single-end FASTQ: lines=\$reads"
                rm "\$r1"
                exit 1
            else
                echo "[INFO] FASTQ file is complete."
                echo "[INFO] GZIP the file"
                # in this case is already named _1.fastq.gz
                r_single="\${srr}_1.fastq"
                # gzip the file
                gzip "\$r_single"
                rm -f *.fastq
            fi
        else
            echo "[ERROR] No FASTQ files found for \$srr"
            exit 1
        fi
    fi

    echo "[DEBUG] Output files:"
    ls -lh *.fastq.gz

    echo "[INFO] Writing version file..."
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasterq_dump: "\$(fasterq-dump --version | sed 's/.*: //')"
        prefetch: "\$(prefetch --version |  sed 's/.*: //')"
        bbmap_repair: "repair.sh from BBMap"
    END_VERSIONS
    """
}
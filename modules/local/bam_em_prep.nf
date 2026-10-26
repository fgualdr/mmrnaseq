/*
 * Filter BAM file
 */
process BAM_EM_PREP {
    tag "$meta.id"
    label 'process_samtools_sort'

    conda (params.enable_conda ? "bedtools=2.31.0 samtools=1.17 coreutils=9.3 gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-9d3a458f6420e5712103ae2af82c94d26d63f059:60b54b43045e8cf39ba307fd683c69d4c57240ce-0':
        'quay.io/biocontainers//mulled-v2-9d3a458f6420e5712103ae2af82c94d26d63f059:60b54b43045e8cf39ba307fd683c69d4c57240ce-0' }"

    input:
    tuple val(meta), path(bam)
    path chrom_sizes

    output:
    tuple val(meta), path("*.read_target_match.txt"), emit: read_target_match
    tuple val(meta), path("*.candidate_table.tsv"), emit: candidate_table
    tuple val(meta), path("*.labeled.nsort.bam"), emit: bam_label
    tuple val(meta), path("*.multi.hotspots.bed"), emit: mm_hotspots
    path "versions.yml"           , emit: versions

    script:
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def overlap = params.label_overlap ? params.label_overlap : 0.25
    // convert to megabytes the memory available
    def avail_mem = task.memory ? task.memory.mega : false
    // sort memory is the available memory minus 1GB divided by the number of cpus
    def sort_memory = avail_mem ? Math.max(256, ((avail_mem - 1000) / task.cpus).intValue()) : 2000
    def single_end = meta.single_end
    def strandedness = 0
    if (meta.strandedness == 'forward') {
        strandedness = 1
    } else if (meta.strandedness == 'reverse') {
        strandedness = 2
    }

    if(strandedness == 0){
        if(single_end){
            """
            dist_r=${params.em_hotspot_distance ?: 2000}

            # Add the HI:i bit to each
            samtools view -h $bam | awk 'BEGIN{OFS="\\t"} \$1~/^@/{print;next} {for(i=12;i<=NF;i++){if(\$i~/^HI:i:/){ \$1=\$1"_"\$i}};print}' | samtools view -bS - > ${prefix}.labeled.bam
            samtools sort -@ $task.cpus -m ${sort_memory}M -n -o ${prefix}.labeled.nsort.bam -T ${prefix}.labeled.nsort ${prefix}.labeled.bam
            bedtools bamtobed -i ${prefix}.labeled.nsort.bam > ${prefix}.labeled.nsort.bed
            
            #################################################################################
            #   GENERATE THE APPROPRIATE FILES FOR EM                                       #
            #################################################################################
        
            # Generate the files required for the EM step:
            # Separate multi from single - this is done on the BAM first then on the BEDPE
            samtools view -h ${prefix}.labeled.nsort.bam | awk 'BEGIN{OFS="\t"} /^@/{print;next} {nh=0; for(i=12;i<=NF;i++){if(\$i~/^NH:i:/){split(\$i,a,":"); nh=a[3]+0}} if(nh>1) print}' | samtools view -bS - > ${prefix}.labeled.multi.bam
            bedtools bamtobed -i ${prefix}.labeled.multi.bam > ${prefix}.labeled.nsort.multi.bed

            # Hotspots
            cut -f 1,2,3 ${prefix}.labeled.nsort.multi.bed > ${prefix}.multi.bed
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.multi.bed > ${prefix}.multi.sort.bed
            bedtools merge -d $dist_r -i ${prefix}.multi.sort.bed > ${prefix}.multi.sort.merged.bed
            awk 'BEGIN {FS="\\t"; OFS="\\t"} {print \$0"\\t""id_"NR}' ${prefix}.multi.sort.merged.bed > ${prefix}.multi.hotspots.temp.bed
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.multi.hotspots.temp.bed > ${prefix}.multi.hotspots.bed
            
            # find intersections - at least 50% of the read has to be within default!
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.labeled.nsort.bed > ${prefix}.labeled.sort.bed

            bedtools intersect \\
                    -a ${prefix}.labeled.sort.bed \\
                    -b ${prefix}.multi.hotspots.bed \\
                    -sorted \\
                    -f ${overlap} \\
                    -wo -loj > ${prefix}.nsort.match.bed
            
            # Replace the no match "." with "id_nomatch"

            awk 'BEGIN{FS=OFS="\\t"} \$10=="."{\$10="id_nomatch"}1' ${prefix}.nsort.match.bed > ${prefix}.nsort.match.fix.bed

            #################################
            #   OUTPUT DATA                 #
            #################################

            # Generate the read_target table for the EM which is composed of read_full_id, read_id_ori, target_id (add id_nomatch)
            awk -F'\\t' 'BEGIN{OFS="\\t"} {id=\$4;sub("_HI:.*", "", \$4); print id, \$4, \$10}' ${prefix}.nsort.match.fix.bed > ${prefix}.read_target.txt
            awk -F'\\t' 'BEGIN{OFS="\\t"} \$3!="id_nomatch"{print \$1, \$2, \$3}' ${prefix}.read_target.txt > ${prefix}.read_target_match.txt
            awk -F'\t' 'BEGIN{OFS="\t"} {id=\$4; orig=\$4; sub("_HI:.*","",orig); if(\$10!="id_nomatch") print id,orig,\$10,\$1,\$2,\$3}' ${prefix}.nsort.match.fix.bed > ${prefix}.candidate_geometry.tsv
            # Build one score record per STAR alignment (HI).
            # The BAM is name-sorted after HI has been appended to QNAME, therefore the
            # two PE mate records for one STAR alignment are adjacent. Collapse them
            # with O(1) memory instead of storing every QNAME in an awk hash.
            samtools view ${prefix}.labeled.nsort.bam | awk 'BEGIN{OFS="\\t"; prev=""}
                {
                    q=\$1; hi=""; nh=""; as=""; nm="";
                    for(i=12;i<=NF;i++){
                        if(\$i~/^HI:i:/){split(\$i,a,":");hi=a[3]}
                        else if(\$i~/^NH:i:/){split(\$i,a,":");nh=a[3]}
                        else if(\$i~/^AS:i:/){split(\$i,a,":");as=a[3]}
                        else if(\$i~/^nM:i:/){split(\$i,a,":");nm=a[3]}
                    }
                    if(q!=prev){print q,hi,nh,as,nm; prev=q}
                }' > ${prefix}.alignment_scores.tsv

            # Join geometry and STAR scores by the HI-labelled alignment ID.
            # sort/join keeps this step streaming and avoids an O(number_of_alignments)
            # in-memory associative array.
            LC_ALL=C sort -k1,1 ${prefix}.alignment_scores.tsv > ${prefix}.alignment_scores.sorted.tsv
            LC_ALL=C sort -k1,1 ${prefix}.candidate_geometry.tsv > ${prefix}.candidate_geometry.sorted.tsv

            LC_ALL=C join -t "\$(printf '\\t')" -v 2 -1 1 -2 1 \\
                ${prefix}.alignment_scores.sorted.tsv \\
                ${prefix}.candidate_geometry.sorted.tsv \\
                > ${prefix}.candidate_missing_scores.tsv

            if [ -s ${prefix}.candidate_missing_scores.tsv ]; then
                echo "ERROR: candidate geometry rows could not be matched to STAR AS/NH/HI/nM records" >&2
                head -20 ${prefix}.candidate_missing_scores.tsv >&2
                exit 2
            fi

            LC_ALL=C join -t "\$(printf '\\t')" -1 1 -2 1 \\
                -o 1.1,2.2,2.3,2.4,2.5,2.6,1.2,1.3,1.4,1.5 \\
                ${prefix}.alignment_scores.sorted.tsv \\
                ${prefix}.candidate_geometry.sorted.tsv \\
                > ${prefix}.candidate_table.body.tsv

            printf 'alignment_id\tfragment_id\ttarget_id\tchrom\tstart\tend\thi\tnh\tas\tnM\n' > ${prefix}.candidate_table.tsv
            cat ${prefix}.candidate_table.body.tsv >> ${prefix}.candidate_table.tsv

            # Fail early if the STAR scoring contract is not present.
            awk -F'\\t' 'NR>1 && (NF!=10 || \$7=="" || \$8=="" || \$9=="" || \$10==""){bad++}
                END{if(bad>0){print "ERROR: " bad " malformed candidate-table rows" > "/dev/stderr"; exit 3}}' \\
                ${prefix}.candidate_table.tsv
    

            samtools_version=\$(samtools --version 2>&1 | awk 'NR==1{print \$2; exit}')
            bedtools_version=\$(bedtools --version 2>&1 | awk 'NR==1{print \$NF; exit}')
            {
                echo "\"${task.process}\":"
                echo "    samtools: \"\$samtools_version\""
                echo "    bedtools: \"\$bedtools_version\""
            } > versions.yml
            """
        }else{
            """
            dist_r=${params.em_hotspot_distance ?: 2000}

            # Add the HI:i bit to each
            samtools view -h $bam | awk 'BEGIN{OFS="\\t"} \$1~/^@/{print;next} {for(i=12;i<=NF;i++){if(\$i~/^HI:i:/){ \$1=\$1"_"\$i}};print}' | samtools view -bS - > ${prefix}.labeled.bam
            samtools sort -@ $task.cpus -m ${sort_memory}M -n -o ${prefix}.labeled.nsort.bam -T ${prefix}.labeled.nsort ${prefix}.labeled.bam
            bedtools bamtobed -i ${prefix}.labeled.nsort.bam -bedpe > ${prefix}.labeled.nsort.bedpe 
            
            #################################################################################
            #   GENERATE THE APPROPRIATE FILES FOR EM                                       #
            #################################################################################
        
            # Generate the files required for the EM step:
            # Separate multi from single - this is done on the BAM first then on the BEDPE
            samtools view -h ${prefix}.labeled.nsort.bam | awk 'BEGIN{OFS="\t"} /^@/{print;next} {nh=0; for(i=12;i<=NF;i++){if(\$i~/^NH:i:/){split(\$i,a,":"); nh=a[3]+0}} if(nh>1) print}' | samtools view -bS - > ${prefix}.labeled.multi.bam
            bedtools bamtobed -i ${prefix}.labeled.multi.bam -bedpe > ${prefix}.labeled.nsort.multi.bedpe 

            # Hotspots
            cut -f 1,2,6 ${prefix}.labeled.nsort.multi.bedpe  > ${prefix}.multi.bed
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.multi.bed > ${prefix}.multi.sort.bed
            bedtools merge -d \$dist_r -i ${prefix}.multi.sort.bed > ${prefix}.multi.sort.merged.bed
            awk 'BEGIN {FS="\\t"; OFS="\t"} {print \$0"\\t""id_"NR}' ${prefix}.multi.sort.merged.bed > ${prefix}.multi.hotspots.temp.bed
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.multi.hotspots.temp.bed > ${prefix}.multi.hotspots.bed
            
            # find intersections - at least 50% of the read has to be within default!

            cut -f 1,2,6,7 ${prefix}.labeled.nsort.bedpe  > ${prefix}.labeled.nsort.bed
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.labeled.nsort.bed > ${prefix}.labeled.sort.bed

            bedtools intersect \\
                    -a ${prefix}.labeled.sort.bed \\
                    -b ${prefix}.multi.hotspots.bed \\
                    -sorted \\
                    -f ${overlap} \\
                    -wo -loj > ${prefix}.nsort.match.bed
            
            # Replace the no match "." with "id_nomatch"

            awk 'BEGIN{FS=OFS="\\t"} \$8=="."{\$8="id_nomatch"}1' ${prefix}.nsort.match.bed > ${prefix}.nsort.match.fix.bed

            #################################
            #   OUTPUT DATA                 #
            #################################

            # Generate the read_target table for the EM which is composed of read_full_id, read_id_ori, target_id (add id_nomatch)
            awk -F'\\t' 'BEGIN{OFS="\\t"} {id=\$4;sub("_HI:.*", "", \$4); print id, \$4, \$8}' ${prefix}.nsort.match.fix.bed > ${prefix}.read_target.txt
            awk -F'\\t' 'BEGIN{OFS="\\t"} \$3!="id_nomatch"{print \$1, \$2, \$3}' ${prefix}.read_target.txt > ${prefix}.read_target_match.txt
            awk -F'\t' 'BEGIN{OFS="\t"} {id=\$4; orig=\$4; sub("_HI:.*","",orig); if(\$8!="id_nomatch") print id,orig,\$8,\$1,\$2,\$3}' ${prefix}.nsort.match.fix.bed > ${prefix}.candidate_geometry.tsv
            # Build one score record per STAR alignment (HI).
            # The BAM is name-sorted after HI has been appended to QNAME, therefore the
            # two PE mate records for one STAR alignment are adjacent. Collapse them
            # with O(1) memory instead of storing every QNAME in an awk hash.
            samtools view ${prefix}.labeled.nsort.bam | awk 'BEGIN{OFS="\\t"; prev=""}
                {
                    q=\$1; hi=""; nh=""; as=""; nm="";
                    for(i=12;i<=NF;i++){
                        if(\$i~/^HI:i:/){split(\$i,a,":");hi=a[3]}
                        else if(\$i~/^NH:i:/){split(\$i,a,":");nh=a[3]}
                        else if(\$i~/^AS:i:/){split(\$i,a,":");as=a[3]}
                        else if(\$i~/^nM:i:/){split(\$i,a,":");nm=a[3]}
                    }
                    if(q!=prev){print q,hi,nh,as,nm; prev=q}
                }' > ${prefix}.alignment_scores.tsv

            # Join geometry and STAR scores by the HI-labelled alignment ID.
            # sort/join keeps this step streaming and avoids an O(number_of_alignments)
            # in-memory associative array.
            LC_ALL=C sort -k1,1 ${prefix}.alignment_scores.tsv > ${prefix}.alignment_scores.sorted.tsv
            LC_ALL=C sort -k1,1 ${prefix}.candidate_geometry.tsv > ${prefix}.candidate_geometry.sorted.tsv

            LC_ALL=C join -t "\$(printf '\\t')" -v 2 -1 1 -2 1 \\
                ${prefix}.alignment_scores.sorted.tsv \\
                ${prefix}.candidate_geometry.sorted.tsv \\
                > ${prefix}.candidate_missing_scores.tsv

            if [ -s ${prefix}.candidate_missing_scores.tsv ]; then
                echo "ERROR: candidate geometry rows could not be matched to STAR AS/NH/HI/nM records" >&2
                head -20 ${prefix}.candidate_missing_scores.tsv >&2
                exit 2
            fi

            LC_ALL=C join -t "\$(printf '\\t')" -1 1 -2 1 \\
                -o 1.1,2.2,2.3,2.4,2.5,2.6,1.2,1.3,1.4,1.5 \\
                ${prefix}.alignment_scores.sorted.tsv \\
                ${prefix}.candidate_geometry.sorted.tsv \\
                > ${prefix}.candidate_table.body.tsv

            printf 'alignment_id\tfragment_id\ttarget_id\tchrom\tstart\tend\thi\tnh\tas\tnM\n' > ${prefix}.candidate_table.tsv
            cat ${prefix}.candidate_table.body.tsv >> ${prefix}.candidate_table.tsv

            # Fail early if the STAR scoring contract is not present.
            awk -F'\\t' 'NR>1 && (NF!=10 || \$7=="" || \$8=="" || \$9=="" || \$10==""){bad++}
                END{if(bad>0){print "ERROR: " bad " malformed candidate-table rows" > "/dev/stderr"; exit 3}}' \\
                ${prefix}.candidate_table.tsv
    

            samtools_version=\$(samtools --version 2>&1 | awk 'NR==1{print \$2; exit}')
            bedtools_version=\$(bedtools --version 2>&1 | awk 'NR==1{print \$NF; exit}')
            {
                echo "\"${task.process}\":"
                echo "    samtools: \"\$samtools_version\""
                echo "    bedtools: \"\$bedtools_version\""
            } > versions.yml
            """
        }
    } else {
        if(single_end){
            """
            ### Check properly !!!!!!!!!!!
            dist_r=${params.em_hotspot_distance ?: 2000}

            # Add the HI:i bit to each
            samtools view -h $bam | awk 'BEGIN{OFS="\\t"} \$1~/^@/{print;next} {for(i=12;i<=NF;i++){if(\$i~/^HI:i:/){ \$1=\$1"_"\$i}};print}' | samtools view -bS - > ${prefix}.labeled.bam
            samtools sort -@ $task.cpus -m ${sort_memory}M -n -o ${prefix}.labeled.nsort.bam -T ${prefix}.labeled.nsort ${prefix}.labeled.bam

            #################################################################################
            # Add the if! so if unstranded - should not split into FW and RW!!!:            #
            #################################################################################
            # Generate the files required for the EM step:
            # Forward strand.
            # 1. alignments mapping to the FW strand (we are single ended: and RNAseq strand oriented)
            samtools view -b -F 16 ${prefix}.labeled.bam > ${prefix}.fwd.bam
            samtools index ${prefix}.fwd.bam

            samtools sort -@ $task.cpus -m ${sort_memory}M -n -o ${prefix}.fwd.nsort.bam -T ${prefix}.fwd.nsort ${prefix}.fwd.bam
            bedtools bamtobed -i ${prefix}.fwd.nsort.bam > ${prefix}.fwd.nsort.bedpe 
            cut -f 1,2,3,4 ${prefix}.fwd.nsort.bedpe > ${prefix}.fwd.nsort.bed

            # Separate multi from single:
            samtools view -h ${prefix}.fwd.nsort.bam | awk 'BEGIN{OFS="\t"} /^@/{print;next} {nh=0; for(i=12;i<=NF;i++){if(\$i~/^NH:i:/){split(\$i,a,":"); nh=a[3]+0}} if(nh>1) print}' | samtools view -bS - > ${prefix}.labeled.multi.fwd.bam
            bedtools bamtobed -i ${prefix}.labeled.multi.fwd.bam > ${prefix}.fwd.nsort.multi.bedpe 

            cut -f 1,2,3,4 ${prefix}.fwd.nsort.multi.bedpe > ${prefix}.fwd.multi.bed
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.fwd.multi.bed > ${prefix}.fwd.multi.sort.bed
            bedtools merge -d \$dist_r -i ${prefix}.fwd.multi.sort.bed > ${prefix}.fwd.multi.sort.merged.bed
            awk 'BEGIN {FS="\\t"; OFS="\\t"} {print \$0"\\t""id_fw_"NR}' ${prefix}.fwd.multi.sort.merged.bed > ${prefix}.fwd.multi.sort.merged.named.tmp.bed
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.fwd.multi.sort.merged.named.tmp.bed > ${prefix}.fwd.multi.sort.merged.named.bed

            # find intersections - at least 50% of the read has to be within default!
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.fwd.nsort.bed > ${prefix}.fwd.sort.bed

            bedtools intersect \\
                    -a ${prefix}.fwd.sort.bed \\
                    -b ${prefix}.fwd.multi.sort.merged.named.bed \\
                    -sorted \\
                    -f ${overlap} \\
                    -wo -loj > ${prefix}.fwd.nsort.match.bed

            # Reverse strand
            # 1. alignments mapping to the RW strand (we are single ended: and RNAseq strand oriented)
            samtools view -b -f 16 ${prefix}.labeled.bam > ${prefix}.rev.bam
            samtools index ${prefix}.rev.bam

            samtools sort -@ $task.cpus -m ${sort_memory}M -n -o ${prefix}.rev.nsort.bam -T ${prefix}.rev.nsort ${prefix}.rev.bam
            bedtools bamtobed -i ${prefix}.rev.nsort.bam > ${prefix}.rev.nsort.bedpe 
            cut -f 1,2,3,4 ${prefix}.rev.nsort.bedpe > ${prefix}.rev.nsort.bed

            # Separate multi from single:
            samtools view -h ${prefix}.rev.nsort.bam | awk 'BEGIN{OFS="\t"} /^@/{print;next} {nh=0; for(i=12;i<=NF;i++){if(\$i~/^NH:i:/){split(\$i,a,":"); nh=a[3]+0}} if(nh>1) print}' | samtools view -bS - > ${prefix}.labeled.multi.rev.bam
            bedtools bamtobed -i ${prefix}.labeled.multi.rev.bam > ${prefix}.rev.nsort.multi.bedpe 

            cut -f 1,2,3,4 ${prefix}.rev.nsort.multi.bedpe > ${prefix}.rev.multi.bed
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.rev.multi.bed > ${prefix}.rev.multi.sort.bed
            bedtools merge -d \$dist_r -i ${prefix}.rev.multi.sort.bed > ${prefix}.rev.multi.sort.merged.bed
            awk 'BEGIN {FS="\\t"; OFS="\\t"} {print \$0"\\t""id_rw_"NR}' ${prefix}.rev.multi.sort.merged.bed > ${prefix}.rev.multi.sort.merged.named.tmp.bed
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.rev.multi.sort.merged.named.tmp.bed > ${prefix}.rev.multi.sort.merged.named.bed

            # find intersections - at least 50% of the read has to be within default!
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.rev.nsort.bed > ${prefix}.rev.sort.bed
            bedtools intersect \\
                    -a ${prefix}.rev.sort.bed \\
                    -b ${prefix}.rev.multi.sort.merged.named.bed \\
                    -sorted \\
                    -f ${overlap} \\
                    -wo -loj > ${prefix}.rev.nsort.match.bed

            #################################
            #   OUTPUT DATA                 #
            #################################
            # Combine hotspots
            cat ${prefix}.fwd.multi.sort.merged.named.bed ${prefix}.rev.multi.sort.merged.named.bed > ${prefix}.multi.hotspots.bed
            # Generate the read_target table for the EM which is composed of read_full_id, read_id_ori, target_id
            cat ${prefix}.fwd.nsort.match.bed ${prefix}.rev.nsort.match.bed > ${prefix}.read_target.bed
            # Replace the no match "." with "id_nomatch"
            awk 'BEGIN{FS=OFS="\\t"} \$8=="."{\$8="id_nomatch"}1' ${prefix}.read_target.bed > ${prefix}.read_target.fix.bed
            awk -F'\\t' 'BEGIN{OFS="\\t"} {id=\$4;sub("_HI:.*", "", \$4); print id, \$4, \$8}' ${prefix}.read_target.fix.bed > ${prefix}.read_target.txt
            awk -F'\\t' 'BEGIN{OFS="\\t"} \$3!="id_nomatch"{print \$1, \$2, \$3}' ${prefix}.read_target.txt > ${prefix}.read_target_match.txt
            awk -F'\t' 'BEGIN{OFS="\t"} {id=\$4; orig=\$4; sub("_HI:.*","",orig); if(\$8!="id_nomatch") print id,orig,\$8,\$1,\$2,\$3}' ${prefix}.read_target.fix.bed > ${prefix}.candidate_geometry.tsv
            # Build one score record per STAR alignment (HI).
            # The BAM is name-sorted after HI has been appended to QNAME, therefore the
            # two PE mate records for one STAR alignment are adjacent. Collapse them
            # with O(1) memory instead of storing every QNAME in an awk hash.
            samtools view ${prefix}.labeled.nsort.bam | awk 'BEGIN{OFS="\\t"; prev=""}
                {
                    q=\$1; hi=""; nh=""; as=""; nm="";
                    for(i=12;i<=NF;i++){
                        if(\$i~/^HI:i:/){split(\$i,a,":");hi=a[3]}
                        else if(\$i~/^NH:i:/){split(\$i,a,":");nh=a[3]}
                        else if(\$i~/^AS:i:/){split(\$i,a,":");as=a[3]}
                        else if(\$i~/^nM:i:/){split(\$i,a,":");nm=a[3]}
                    }
                    if(q!=prev){print q,hi,nh,as,nm; prev=q}
                }' > ${prefix}.alignment_scores.tsv

            # Join geometry and STAR scores by the HI-labelled alignment ID.
            # sort/join keeps this step streaming and avoids an O(number_of_alignments)
            # in-memory associative array.
            LC_ALL=C sort -k1,1 ${prefix}.alignment_scores.tsv > ${prefix}.alignment_scores.sorted.tsv
            LC_ALL=C sort -k1,1 ${prefix}.candidate_geometry.tsv > ${prefix}.candidate_geometry.sorted.tsv

            LC_ALL=C join -t "\$(printf '\\t')" -v 2 -1 1 -2 1 \\
                ${prefix}.alignment_scores.sorted.tsv \\
                ${prefix}.candidate_geometry.sorted.tsv \\
                > ${prefix}.candidate_missing_scores.tsv

            if [ -s ${prefix}.candidate_missing_scores.tsv ]; then
                echo "ERROR: candidate geometry rows could not be matched to STAR AS/NH/HI/nM records" >&2
                head -20 ${prefix}.candidate_missing_scores.tsv >&2
                exit 2
            fi

            LC_ALL=C join -t "\$(printf '\\t')" -1 1 -2 1 \\
                -o 1.1,2.2,2.3,2.4,2.5,2.6,1.2,1.3,1.4,1.5 \\
                ${prefix}.alignment_scores.sorted.tsv \\
                ${prefix}.candidate_geometry.sorted.tsv \\
                > ${prefix}.candidate_table.body.tsv

            printf 'alignment_id\tfragment_id\ttarget_id\tchrom\tstart\tend\thi\tnh\tas\tnM\n' > ${prefix}.candidate_table.tsv
            cat ${prefix}.candidate_table.body.tsv >> ${prefix}.candidate_table.tsv

            # Fail early if the STAR scoring contract is not present.
            awk -F'\\t' 'NR>1 && (NF!=10 || \$7=="" || \$8=="" || \$9=="" || \$10==""){bad++}
                END{if(bad>0){print "ERROR: " bad " malformed candidate-table rows" > "/dev/stderr"; exit 3}}' \\
                ${prefix}.candidate_table.tsv

            samtools_version=\$(samtools --version 2>&1 | awk 'NR==1{print \$2; exit}')
            bedtools_version=\$(bedtools --version 2>&1 | awk 'NR==1{print \$NF; exit}')
            {
                echo "\"${task.process}\":"
                echo "    samtools: \"\$samtools_version\""
                echo "    bedtools: \"\$bedtools_version\""
            } > versions.yml
            """
        }else{
            """
            
            dist_r=${params.em_hotspot_distance ?: 2000}

            # Add the HI:i bit to each
            samtools view -h $bam | awk 'BEGIN{OFS="\\t"} \$1~/^@/{print;next} {for(i=12;i<=NF;i++){if(\$i~/^HI:i:/){ \$1=\$1"_"\$i}};print}' | samtools view -bS - > ${prefix}.labeled.bam
            samtools sort -@ $task.cpus -m ${sort_memory}M -n -o ${prefix}.labeled.nsort.bam -T ${prefix}.labeled.nsort ${prefix}.labeled.bam

            #################################################################################
            # Add the if! so if unstranded - should not split into FW and RW!!!:            #
            #################################################################################
            # Generate the files required for the EM step:
            # Forward strand.
            # 1. alignments of the second in pair if they map to the forward strand
            samtools view -b -f 128 -F 16 ${prefix}.labeled.bam > fwd1.bam
            samtools index fwd1.bam
            samtools view -b -f 80 ${prefix}.labeled.bam > fwd2.bam
            samtools index fwd2.bam
            samtools merge -f ${prefix}.fwd.bam fwd1.bam fwd2.bam
            samtools index ${prefix}.fwd.bam

            samtools sort -@ $task.cpus -m ${sort_memory}M -n -o ${prefix}.fwd.nsort.bam -T ${prefix}.fwd.nsort ${prefix}.fwd.bam
            bedtools bamtobed -i ${prefix}.fwd.nsort.bam -bedpe > ${prefix}.fwd.nsort.bedpe 
            cut -f 1,2,6,7 ${prefix}.fwd.nsort.bedpe > ${prefix}.fwd.nsort.bed

            # Separate multi from single:
            samtools view -h ${prefix}.fwd.nsort.bam | awk 'BEGIN{OFS="\t"} /^@/{print;next} {nh=0; for(i=12;i<=NF;i++){if(\$i~/^NH:i:/){split(\$i,a,":"); nh=a[3]+0}} if(nh>1) print}' | samtools view -bS - > ${prefix}.labeled.multi.fwd.bam
            bedtools bamtobed -i ${prefix}.labeled.multi.fwd.bam -bedpe > ${prefix}.fwd.nsort.multi.bedpe 

            cut -f 1,2,6 ${prefix}.fwd.nsort.multi.bedpe > ${prefix}.fwd.multi.bed
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.fwd.multi.bed > ${prefix}.fwd.multi.sort.bed
            bedtools merge -d \$dist_r -i ${prefix}.fwd.multi.sort.bed > ${prefix}.fwd.multi.sort.merged.bed
            awk 'BEGIN {FS="\\t"; OFS="\\t"} {print \$0"\\t""id_fw_"NR}' ${prefix}.fwd.multi.sort.merged.bed > ${prefix}.fwd.multi.sort.merged.named.tmp.bed
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.fwd.multi.sort.merged.named.tmp.bed > ${prefix}.fwd.multi.sort.merged.named.bed

            # find intersections - at least 50% of the read has to be within default!
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.fwd.nsort.bed > ${prefix}.fwd.sort.bed

            bedtools intersect \\
                    -a ${prefix}.fwd.sort.bed \\
                    -b ${prefix}.fwd.multi.sort.merged.named.bed \\
                    -sorted \\
                    -f ${overlap} \\
                    -wo -loj > ${prefix}.fwd.nsort.match.bed

            # Reverse strand
            # 1. alignments of the second in pair if they map to the reverse strand
            # 2. alignments of the first in pair if they map to the forward strand
            samtools view -b -f 144 ${prefix}.labeled.bam > rev1.bam
            samtools index rev1.bam
            samtools view -b -f 64 -F 16 ${prefix}.labeled.bam > rev2.bam
            samtools index rev2.bam
            samtools merge -f ${prefix}.rev.bam rev1.bam rev2.bam
            samtools index ${prefix}.rev.bam

            samtools sort -@ $task.cpus -m ${sort_memory}M -n -o ${prefix}.rev.nsort.bam -T ${prefix}.rev.nsort ${prefix}.rev.bam
            bedtools bamtobed -i ${prefix}.rev.nsort.bam -bedpe > ${prefix}.rev.nsort.bedpe 
            cut -f 1,2,6,7 ${prefix}.rev.nsort.bedpe > ${prefix}.rev.nsort.bed

            # Separate multi from single:
            samtools view -h ${prefix}.rev.nsort.bam | awk 'BEGIN{OFS="\t"} /^@/{print;next} {nh=0; for(i=12;i<=NF;i++){if(\$i~/^NH:i:/){split(\$i,a,":"); nh=a[3]+0}} if(nh>1) print}' | samtools view -bS - > ${prefix}.labeled.multi.rev.bam
            bedtools bamtobed -i ${prefix}.labeled.multi.rev.bam -bedpe > ${prefix}.rev.nsort.multi.bedpe 

            cut -f 1,2,6 ${prefix}.rev.nsort.multi.bedpe > ${prefix}.rev.multi.bed
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.rev.multi.bed > ${prefix}.rev.multi.sort.bed
            bedtools merge -d \$dist_r -i ${prefix}.rev.multi.sort.bed > ${prefix}.rev.multi.sort.merged.bed
            awk 'BEGIN {FS="\\t"; OFS="\\t"} {print \$0"\\t""id_rw_"NR}' ${prefix}.rev.multi.sort.merged.bed > ${prefix}.rev.multi.sort.merged.named.tmp.bed
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.rev.multi.sort.merged.named.tmp.bed > ${prefix}.rev.multi.sort.merged.named.bed

            # find intersections - at least 50% of the read has to be within default!
            sort --parallel=$task.cpus -k1,1 -k2,2n ${prefix}.rev.nsort.bed > ${prefix}.rev.sort.bed
            bedtools intersect \\
                    -a ${prefix}.rev.sort.bed \\
                    -b ${prefix}.rev.multi.sort.merged.named.bed \\
                    -sorted \\
                    -f ${overlap} \\
                    -wo -loj > ${prefix}.rev.nsort.match.bed

            #################################
            #   OUTPUT DATA                 #
            #################################
            # Combine hotspots
            cat ${prefix}.fwd.multi.sort.merged.named.bed ${prefix}.rev.multi.sort.merged.named.bed > ${prefix}.multi.hotspots.bed
            # Generate the read_target table for the EM which is composed of read_full_id, read_id_ori, target_id
            cat ${prefix}.fwd.nsort.match.bed ${prefix}.rev.nsort.match.bed > ${prefix}.read_target.bed
            # Replace the no match "." with "id_nomatch"
            awk 'BEGIN{FS=OFS="\\t"} \$8=="."{\$8="id_nomatch"}1' ${prefix}.read_target.bed > ${prefix}.read_target.fix.bed
            awk -F'\\t' 'BEGIN{OFS="\\t"} {id=\$4;sub("_HI:.*", "", \$4); print id, \$4, \$8}' ${prefix}.read_target.fix.bed > ${prefix}.read_target.txt
            awk -F'\\t' 'BEGIN{OFS="\\t"} \$3!="id_nomatch"{print \$1, \$2, \$3}' ${prefix}.read_target.txt > ${prefix}.read_target_match.txt
            awk -F'\t' 'BEGIN{OFS="\t"} {id=\$4; orig=\$4; sub("_HI:.*","",orig); if(\$8!="id_nomatch") print id,orig,\$8,\$1,\$2,\$3}' ${prefix}.read_target.fix.bed > ${prefix}.candidate_geometry.tsv
            # Build one score record per STAR alignment (HI).
            # The BAM is name-sorted after HI has been appended to QNAME, therefore the
            # two PE mate records for one STAR alignment are adjacent. Collapse them
            # with O(1) memory instead of storing every QNAME in an awk hash.
            samtools view ${prefix}.labeled.nsort.bam | awk 'BEGIN{OFS="\\t"; prev=""}
                {
                    q=\$1; hi=""; nh=""; as=""; nm="";
                    for(i=12;i<=NF;i++){
                        if(\$i~/^HI:i:/){split(\$i,a,":");hi=a[3]}
                        else if(\$i~/^NH:i:/){split(\$i,a,":");nh=a[3]}
                        else if(\$i~/^AS:i:/){split(\$i,a,":");as=a[3]}
                        else if(\$i~/^nM:i:/){split(\$i,a,":");nm=a[3]}
                    }
                    if(q!=prev){print q,hi,nh,as,nm; prev=q}
                }' > ${prefix}.alignment_scores.tsv

            # Join geometry and STAR scores by the HI-labelled alignment ID.
            # sort/join keeps this step streaming and avoids an O(number_of_alignments)
            # in-memory associative array.
            LC_ALL=C sort -k1,1 ${prefix}.alignment_scores.tsv > ${prefix}.alignment_scores.sorted.tsv
            LC_ALL=C sort -k1,1 ${prefix}.candidate_geometry.tsv > ${prefix}.candidate_geometry.sorted.tsv

            LC_ALL=C join -t "\$(printf '\\t')" -v 2 -1 1 -2 1 \\
                ${prefix}.alignment_scores.sorted.tsv \\
                ${prefix}.candidate_geometry.sorted.tsv \\
                > ${prefix}.candidate_missing_scores.tsv

            if [ -s ${prefix}.candidate_missing_scores.tsv ]; then
                echo "ERROR: candidate geometry rows could not be matched to STAR AS/NH/HI/nM records" >&2
                head -20 ${prefix}.candidate_missing_scores.tsv >&2
                exit 2
            fi

            LC_ALL=C join -t "\$(printf '\\t')" -1 1 -2 1 \\
                -o 1.1,2.2,2.3,2.4,2.5,2.6,1.2,1.3,1.4,1.5 \\
                ${prefix}.alignment_scores.sorted.tsv \\
                ${prefix}.candidate_geometry.sorted.tsv \\
                > ${prefix}.candidate_table.body.tsv

            printf 'alignment_id\tfragment_id\ttarget_id\tchrom\tstart\tend\thi\tnh\tas\tnM\n' > ${prefix}.candidate_table.tsv
            cat ${prefix}.candidate_table.body.tsv >> ${prefix}.candidate_table.tsv

            # Fail early if the STAR scoring contract is not present.
            awk -F'\\t' 'NR>1 && (NF!=10 || \$7=="" || \$8=="" || \$9=="" || \$10==""){bad++}
                END{if(bad>0){print "ERROR: " bad " malformed candidate-table rows" > "/dev/stderr"; exit 3}}' \\
                ${prefix}.candidate_table.tsv

            samtools_version=\$(samtools --version 2>&1 | awk 'NR==1{print \$2; exit}')
            bedtools_version=\$(bedtools --version 2>&1 | awk 'NR==1{print \$NF; exit}')
            {
                echo "\"${task.process}\":"
                echo "    samtools: \"\$samtools_version\""
                echo "    bedtools: \"\$bedtools_version\""
            } > versions.yml
            """
        }
    }


}

//
// Check input samplesheet and get read channels
//

include { SRX_DOWNLOAD } from '../../modules/local/srx_sra_download'

workflow FASTQ_FROM_SRA {
    take:
    public_data_ids // file: /path/to/samplesheet.csv

    main:
    // first we modify the samplesheet to make a new channel of meta:
    Channel.from(public_data_ids)
        .splitCsv ( header:true, sep:',' )
        .map { create_meta_channel(it) }
        .set { ch_public_data_ids }

    SRX_DOWNLOAD ( ch_public_data_ids )
    SRX_DOWNLOAD.out.fastq                              // channel: tuple val(meta), path('*.fastq.gz'), emit: fastq from the SRX_DOWNLOAD module - so fastq can be of length 1 or 2
        .map { check_fastq_files(it[0], it[1]) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SRX_DOWNLOAD.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta ]
def create_meta_channel(LinkedHashMap row) {
    // create meta map

    // filed will match the header of the csv: sample_accession,experiment_accession,run_accession,strandedness,sample,replicate,study_accession
    // for now it works by SRR
    def meta = [:]
    meta.id           = row.sample + "_" + row.replicate + "_" + row.run_accession
    meta.replicate    = row.replicate
    meta.strandedness = row.strandedness
    meta.study_accession        = row.study_accession
    meta.experiment_accession   = row.experiment_accession
    meta.run_accession          = row.run_accession
    meta.sample_accession       = row.sample_accession
    return meta
}

// Function to check for files ending in '_1.fastq.gz or _2.fastq.gz' and generating: [ meta + [ single_end: true/false], [ reads_1, reads_2 ] ]

def check_fastq_files(meta, fastq) {
    def meta_new = [:]
    meta_new.id = meta.id
    meta_new.strandedness = meta.strandedness
    
    def fastq_meta = []
    println(fastq)
    fastq_string = fastq.toString()
    println(fastq_string)
    // we now check for the fastq files
    def reads_1 = fastq_string.findAll(/_1.fastq.gz/)
    def reads_2 = fastq_string.findAll(/_2.fastq.gz/)

    if (reads_1 && reads_2) {
        if (reads_1.size() == 1 && reads_2.size() == 1) {
            // Exactly one _1.fastq.gz and one _2.fastq.gz file
            reads_1 = fastq[0]
            reads_2 = fastq[1]
            fastq_meta = [ meta_new + [ single_end: false], [ reads_1, reads_2 ] ]
            println "DEBUG: ${meta_new.id} reads_1: ${reads_1} reads_2: ${reads_2}"
        } else {
            // Handle other cases or print debug information
            println "DEBUG: Unexpected size for reads_1 or reads_2"
        }
    } else if (reads_1) {
        // Only _1.fastq.gz file
        reads_1 = fastq
        fastq_meta = [ meta_new + [ single_end: true], [ reads_1 ] ]
        println "DEBUG: ${meta_new.id} reads_1: ${reads_1}"
    } else if (reads_2) {
        // Exit with ERROR: Only _2.fastq.gz file
        exit 1, "ERROR: ${meta_new.id} , ${reads_2}, has only _2.fastq.gz files"
    } else {
        // Exit with ERROR: No _1.fastq.gz or _2.fastq.gz files
        exit 1, "ERROR: ${meta_new.id} , has no _1.fastq.gz or _2.fastq.gz files"
    }
    return fastq_meta
    
}


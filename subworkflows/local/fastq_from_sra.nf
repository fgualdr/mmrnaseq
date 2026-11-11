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
// One SRX is an experiment; multiple SRRs belonging to that SRX are technical runs.
// The biological sample and replicate labels are retained, but SRXs are not
// automatically merged merely because their labels happen to match.
def create_meta_channel(LinkedHashMap row) {
    for (key in ['sample','replicate','run_accession','experiment_accession','strandedness']) {
        if (!row[key]?.toString()?.trim()) {
            throw new IllegalArgumentException("Missing ${key} in public-data samplesheet: ${row}")
        }
    }
    def groupId = "${row.sample}_${row.replicate}_${row.experiment_accession}"
    def meta = [
        id: "${groupId}_${row.run_accession}",
        group_id: groupId,
        replicate: row.replicate,
        strandedness: row.strandedness,
        experiment_accession: row.experiment_accession,
        run_accession: row.run_accession,
        sample_accession: row.sample_accession,
        study_accession: row.study_accession
    ]
    return meta
}

def check_fastq_files(meta, fastq) {
    def paths = fastq instanceof Collection ? fastq.toList() : [fastq]
    def mate1 = paths.findAll { it.toString().endsWith('_1.fastq.gz') }
    def mate2 = paths.findAll { it.toString().endsWith('_2.fastq.gz') }
    def single = paths.findAll { it.toString().endsWith('.fastq.gz') &&
        !it.toString().endsWith('_1.fastq.gz') && !it.toString().endsWith('_2.fastq.gz') }
    def reads
    boolean isSingle
    if (mate1.size() == 1 && mate2.size() == 1 && paths.size() == 2) {
        reads = [mate1[0], mate2[0]]
        isSingle = false
    } else if (single.size() == 1 && paths.size() == 1) {
        reads = single
        isSingle = true
    } else if (mate1.size() == 1 && paths.size() == 1) {
        reads = mate1
        isSingle = true
    } else {
        throw new IllegalArgumentException("Unexpected SRA FASTQ layout for ${meta.run_accession}: ${paths}")
    }
    return [meta + [single_end: isSingle], reads]
}

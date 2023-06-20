//
// SQUIRE
//

include { SQUIRE_COUNT     } from '../../modules/local/squire_count'

workflow SQUIRE {
    take:
    ch_bam_bai   // channel: [ val(meta), [ bam ] ]
    ch_gtf
    clean_folder
    ch_chrom_sizes
    ch_read_length


    main: 

    ch_versions = Channel.empty()

    // 2) Apply SQUIRE_COUNT
    SQUIRE_COUNT( 
        ch_bam_bai, 
        ch_gtf,
        clean_folder,
        ch_read_length
        )
    ch_versions = ch_versions.mix(SQUIRE_COUNT.out.versions)

    emit:
    count_out      = SQUIRE_COUNT.out.count_folder // this is the folder with the count files - should be named according to meta.id
    versions = ch_versions                     // channel: [ versions.yml ]
}


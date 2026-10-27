//
// Filter BAM via Bayesian EM algorithm 


include { BAM_FILTER } from '../../modules/local/bam_filter'
include { MM_CANDIDATE_PREP } from '../../modules/local/mm_candidate_prep'
include { BAM_EM_PREP } from '../../modules/local/bam_em_prep'
include { BAM_EM } from '../../modules/local/bam_em'
include { BAM_EM_OUT } from '../../modules/local/bam_em_out'


workflow BAM_FILTER_EM {
    take:
    ch_bam_bai   // channel: [ val(meta), [ bam ] ]
    ch_chrom_sizes


    main:

    ch_versions = Channel.empty()

    // Basic mapped/proper-pair filtering. Secondary STAR candidates are retained.
    BAM_FILTER(ch_bam_bai)
    ch_versions = ch_versions.mix(BAM_FILTER.out.versions.first())

    // Exact fragment-candidate contract used by the rebuilt resolver. PATCH 01 is
    // non-disruptive: the legacy resolver still controls the BAM until PATCH 02.
    MM_CANDIDATE_PREP(BAM_FILTER.out.bam)
    ch_versions = ch_versions.mix(MM_CANDIDATE_PREP.out.versions.first())

    if( !params.skip_em  ) {

        // Prepare for EM
        BAM_EM_PREP(BAM_FILTER.out.bam,ch_chrom_sizes)
        ch_versions = ch_versions.mix(BAM_EM_PREP.out.versions.first())

        BAM_EM(BAM_EM_PREP.out.read_target_match)
        ch_versions = ch_versions.mix(BAM_EM.out.versions.first())
        // This needs to be combined in order to avoid error!!!!!!
        BAM_EM_PREP
            .out
            .bam_label
            .join(BAM_EM.out.final_bedpe, by: [0])
            .set { ch_em_bam }

        // Make a final BAM cleaned with unique mappers
        BAM_EM_OUT( ch_em_bam,
                    ch_chrom_sizes
                    ) 

        ch_versions = ch_versions.mix(BAM_EM_OUT.out.versions.first())
        ch_filter_bam = BAM_EM_OUT.out.bam
    } else {
        ch_filter_bam = BAM_FILTER.out.bam
    }

    emit:
    bam      = ch_filter_bam
    candidate_table = params.skip_em ? Channel.empty() : BAM_EM_PREP.out.candidate_table
    mm_candidates = MM_CANDIDATE_PREP.out.candidates
    mm_unique_fragments = MM_CANDIDATE_PREP.out.unique_fragments
    mm_unique_anchors = MM_CANDIDATE_PREP.out.unique_anchors
    mm_excluded = MM_CANDIDATE_PREP.out.excluded
    mm_contract_qc = MM_CANDIDATE_PREP.out.qc
    versions = ch_versions                     // channel: [ versions.yml ]
}


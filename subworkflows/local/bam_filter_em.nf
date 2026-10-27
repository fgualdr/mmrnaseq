// STAR-candidate multimapper resolver.
// PATCH 02-PY computes posterior assignments with a streaming NumPy/Python
// implementation but keeps the filtered STAR BAM unchanged.

include { BAM_FILTER }        from '../../modules/local/bam_filter'
include { MM_CANDIDATE_PREP } from '../../modules/local/mm_candidate_prep'
include { MM_BAYES_RESOLVER } from '../../modules/local/mm_bayes_resolver'

workflow BAM_FILTER_EM {
    take:
    ch_bam_bai
    ch_chrom_sizes

    main:
    ch_versions = Channel.empty()

    BAM_FILTER(ch_bam_bai)
    ch_versions = ch_versions.mix(BAM_FILTER.out.versions.first())
    ch_filter_bam = BAM_FILTER.out.bam

    MM_CANDIDATE_PREP(BAM_FILTER.out.bam)
    ch_versions = ch_versions.mix(MM_CANDIDATE_PREP.out.versions.first())

    if (!params.skip_em) {
        MM_CANDIDATE_PREP.out.candidates
            .join(MM_CANDIDATE_PREP.out.unique_fragments, by: [0])
            .join(MM_CANDIDATE_PREP.out.unique_anchors, by: [0])
            .set { ch_resolver_inputs }

        resolver_script = file("${projectDir}/bin/mm_bayes_resolver.py", checkIfExists: true)

        MM_BAYES_RESOLVER(
            ch_resolver_inputs,
            resolver_script,
            ch_chrom_sizes
        )
        ch_versions = ch_versions.mix(MM_BAYES_RESOLVER.out.versions.first())
    }

    emit:
    bam = ch_filter_bam
    mm_candidates = MM_CANDIDATE_PREP.out.candidates
    mm_unique_fragments = MM_CANDIDATE_PREP.out.unique_fragments
    mm_unique_anchors = MM_CANDIDATE_PREP.out.unique_anchors
    mm_excluded = MM_CANDIDATE_PREP.out.excluded
    mm_contract_qc = MM_CANDIDATE_PREP.out.qc
    mm_posterior = params.skip_em ? Channel.empty() : MM_BAYES_RESOLVER.out.posterior
    mm_resolution = params.skip_em ? Channel.empty() : MM_BAYES_RESOLVER.out.resolution
    mm_resolver_qc = params.skip_em ? Channel.empty() : MM_BAYES_RESOLVER.out.qc
    mm_density_model = params.skip_em ? Channel.empty() : MM_BAYES_RESOLVER.out.density_model
    versions = ch_versions
}

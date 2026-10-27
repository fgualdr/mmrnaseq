// STAR-candidate Bayesian multimapper resolver with three explicit BAM tiers.
//
// k / density selection belongs entirely to MM_BAYES_RESOLVER and is learned
// per sample (and per observed orientation). This realization stage is therefore
// agnostic to the selected k/grid and remains compatible with later adaptive-k
// refinements as long as resolution.tsv retains its contract.

include { BAM_FILTER }        from '../../modules/local/bam_filter'
include { MM_CANDIDATE_PREP } from '../../modules/local/mm_candidate_prep'
include { MM_BAYES_RESOLVER } from '../../modules/local/mm_bayes_resolver'
include { MM_BAM_REALIZE }    from '../../modules/local/mm_bam_realize'
include { MM_QC_MERGE }       from '../../modules/local/mm_qc_merge'

workflow BAM_FILTER_EM {
    take:
    ch_bam_bai
    ch_chrom_sizes

    main:
    ch_versions = Channel.empty()

    BAM_FILTER(ch_bam_bai)
    ch_versions = ch_versions.mix(BAM_FILTER.out.versions.first())

    MM_CANDIDATE_PREP(BAM_FILTER.out.bam)
    ch_versions = ch_versions.mix(MM_CANDIDATE_PREP.out.versions.first())

    if (!params.skip_em) {
        MM_CANDIDATE_PREP.out.candidates
            .join(MM_CANDIDATE_PREP.out.unique_fragments, by: [0])
            .join(MM_CANDIDATE_PREP.out.unique_anchors, by: [0])
            .set { ch_resolver_inputs }

        resolver_script = file("${projectDir}/bin/mm_bayes_resolver.py", checkIfExists: true)
        MM_BAYES_RESOLVER(ch_resolver_inputs, resolver_script, ch_chrom_sizes)
        ch_versions = ch_versions.mix(MM_BAYES_RESOLVER.out.versions.first())

        BAM_FILTER.out.bam
            .join(MM_BAYES_RESOLVER.out.resolution, by: [0])
            .join(MM_CANDIDATE_PREP.out.excluded, by: [0])
            .set { ch_realize_inputs }

        MM_BAM_REALIZE(ch_realize_inputs)
        ch_versions = ch_versions.mix(MM_BAM_REALIZE.out.versions.first())

        ch_quantitative_bam = MM_BAM_REALIZE.out.quantitative_bam
        ch_highconf_bam     = MM_BAM_REALIZE.out.highconf_bam
        ch_full_bam         = MM_BAM_REALIZE.out.full_bam

        // Merge and validate PATCH01 contract, PATCH02P posterior-resolution,
        // PATCH03 BAM realization, and sample/orientation-specific density QC.
        MM_CANDIDATE_PREP.out.qc
            .join(MM_BAYES_RESOLVER.out.qc, by: [0])
            .join(MM_BAM_REALIZE.out.qc, by: [0])
            .join(MM_BAYES_RESOLVER.out.density_model, by: [0])
            .join(MM_BAYES_RESOLVER.out.density_cv, by: [0])
            .set { ch_mm_qc_inputs }

        MM_QC_MERGE(ch_mm_qc_inputs)
        ch_versions = ch_versions.mix(MM_QC_MERGE.out.versions.first())
    } else {
        ch_quantitative_bam = BAM_FILTER.out.bam
        ch_highconf_bam     = BAM_FILTER.out.bam
        ch_full_bam         = Channel.empty()
    }

    emit:
    // Backward-compatible generic output means the quantitative, globally
    // expected-FDR-controlled set. Do not silently substitute the stricter tier.
    bam = ch_quantitative_bam
    bam_quantitative = ch_quantitative_bam
    bam_highconf = ch_highconf_bam
    bam_confident = ch_highconf_bam
    bam_full = ch_full_bam

    mm_candidates = MM_CANDIDATE_PREP.out.candidates
    mm_unique_fragments = MM_CANDIDATE_PREP.out.unique_fragments
    mm_unique_anchors = MM_CANDIDATE_PREP.out.unique_anchors
    mm_excluded = MM_CANDIDATE_PREP.out.excluded
    mm_contract_qc = MM_CANDIDATE_PREP.out.qc
    mm_posterior = params.skip_em ? Channel.empty() : MM_BAYES_RESOLVER.out.posterior
    mm_resolution = params.skip_em ? Channel.empty() : MM_BAYES_RESOLVER.out.resolution
    mm_resolver_qc = params.skip_em ? Channel.empty() : MM_BAYES_RESOLVER.out.qc
    mm_density_model = params.skip_em ? Channel.empty() : MM_BAYES_RESOLVER.out.density_model
    mm_density_cv = params.skip_em ? Channel.empty() : MM_BAYES_RESOLVER.out.density_cv
    mm_realization_qc = params.skip_em ? Channel.empty() : MM_BAM_REALIZE.out.qc
    mm_qc = params.skip_em ? Channel.empty() : MM_QC_MERGE.out.qc
    mm_multiqc = params.skip_em ? Channel.empty() : MM_QC_MERGE.out.multiqc
    versions = ch_versions
}

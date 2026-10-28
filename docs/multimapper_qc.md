# Multimapper and fragment-level QC

mmrnaseq distinguishes **biological fragments** from SAM/BAM alignment records. For paired-end libraries, one read pair is one fragment even though a selected alignment normally contributes two BAM records.

## Quantitative fragment usability

`--min_mapped_reads` is evaluated **after multimapper handling**, using STAR input fragments as the denominator.

- With Bayesian resolution enabled (`--skip_em false`):
  `quantitative usable = trusted unique + globally expected-FDR-resolved multimappers`.
- With Bayesian resolution disabled (`--skip_em true`):
  `quantitative usable = trusted unique only`.

STAR's own unique, multiple-loci, and too-many-loci percentages are still reported, but they are aligner-level availability metrics and do not determine downstream sample acceptance by themselves.

The MultiQC **Quantitative fragment usability** table reports STAR input and mapping categories together with PATCH01 trusted-unique / eligible-multimapper accounting, retained resolved multimappers, unresolved multimappers, excluded mapped fragments, the final quantitative-usable percentage, and PASS/FAIL status.

## Bayesian resolver

When Bayesian resolution is enabled, the **Multimapper resolver** table reports fragment-level candidate multiplicity, STAR-score ambiguity, posterior confidence, expected-FDR performance, absolute-score failures, resolution rate, and the three BAM realization tiers.

The BAM tiers have different meanings:

- `*.mm.quantitative.sorted.bam`: trusted unique fragments plus all multimappers accepted by the global posterior expected-FDR procedure. This is the primary downstream quantitative BAM.
- `high_confidence/*.mm.highconf.sorted.bam`: trusted unique fragments plus the stricter per-fragment high-confidence subset (`pmax >= 1 - mm_resolution_fdr`). This is a robustness/QC product.
- `full_ambiguity/*.mm.full.sorted.bam`: trusted unique fragments plus unresolved alternatives and excluded archival records. This preserves ambiguity/provenance and is not a direct counting BAM.

## Density bandwidth diagnostics

The **Multimapper density bandwidth** table reports the sample- and observed-orientation-specific adaptive unique-anchor density fit, including selected `k`, the maximum tested `k`, boundary status, cross-validation likelihood, gain relative to the previous grid point, and physical genomic radii implied by the selected neighborhood.

A boundary hit means only that the selected `k` equals the largest value currently tested. It is a diagnostic to review; it is not automatically a failure and PATCH07 does not change the resolver's `k` grid or posterior assignments.

## Preseq and final BAM QC

Preseq and the main samtools QC are run on the **quantitative BAM**, not the pre-resolution STAR BAM, so reported library-complexity and alignment statistics correspond to the data actually used for downstream locus-level analysis.

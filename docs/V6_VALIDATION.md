# v6 HPC comparison and release requirements

Each stage should run on the **same biological samples** (at least two PE replicates, including one with substantial multimapping), reference build, STAR index, GTF, RepeatMasker file, and versions. Archive the baseline outputs before proceeding. Different samples or references cannot establish equivalence.

## Per-stage gates

- STEP 01: inspect `fastp.json` mate/fragment arithmetic, UMI-transformed headers, lane merging, unique sample IDs. Adjust `--min_trimmed_reads` threshold if previously expressed in individual mates.
- STEP 02: inspect MultiQC HTML and `multiqc_data`: custom general stats and fragment fate columns, STAR fragment conservation, report navigation order. Native flagstat files remain outside MultiQC; don't compare pair counts with alignment-record counts.
- STEP 03: for each sample and region, verify `normalized[g,s]=raw[g,s]*scaling[s]`, no missing sample IDs, and invariant ALL scaling/PCA. A change in normalization factors requires investigation.
- STEP 04: `ALL/INTRON/EXON` matrices must agree with STEP 03 exactly; `RMSK` RDS (counts, flank lengths, column metadata) must agree exactly for the same inputs. Re-enable `--count_utr true` to check old UTR values. `--count_repeats false` must skip RMSK cleanly.
- STEP 05: confirm no featureCounts/biotype QC tasks; main count and MultiQC metrics unchanged.
- STEP 06: compare QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, NH/HI/AS/posterior tags in final quantitative BAMs; compare candidate/resolver/realization TSVs and all matrices. Exercise `--skip_em true` and `false` and UMI mode if supported. Coordinate-sorted BAM may differ in compressed bytes or tie ordering: compare records, not BAM hashes.
- STEP 07: `nextflow config -flat` resolves without fetching remote nf-core configs, standard non-auto samples do not trigger transcript FASTA creation; test `strandedness=auto`, optional Salmon quantification, and explicit `--prepare_transcript_fasta true`.
- STEP 08: all tests pass, production smoke test passes, MultiQC static PNG/SVG/PDF plots are present, and no numerical regressions remain unexplained.

## General sanity

Expected template conservation:

`N_STAR_input = N_unique_STAR + N_multi_STAR + N_too_many_STAR + N_unmapped_STAR`

`N_contract = N_unique_trusted + N_multi_eligible + N_excluded`

`N_multi_eligible = N_resolved + N_unresolved`

`N_quantitative = N_unique_trusted + N_resolved`

`N_not_retained = N_STAR_input - N_quantitative`

The `N_STAR_input - N_contract` category is outside the filter contract. It is not synonymous with STAR-unmapped reads, and should not be relabeled as such. STAR PE singletons are unavailable from the current summary contract; do not fabricate a singleton count.

## Final acceptance

Stop and investigate if there is any mismatch in individual gene/repeat counts, QNAME-level multimapper decisions, fragment conservation, or newly missing MultiQC sections. Execution without failure is necessary but insufficient. Record wall time, peak RSS and disk write/size of each stage on the same HPC nodes when possible.

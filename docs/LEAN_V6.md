# mmRNAseq v6 — Lean execution and fragment accounting

This is an incremental refactor of the existing mmRNAseq DSL2 workflow, not a replacement of its Bayesian multimapper inference.

## Normal operation

- IEO input files: biological sample identifier includes replicate (`_R1`, `_R2`), and only the generated technical-lane `_T1`, `_T2`, etc. suffix is stripped for grouping.
- Public SRA input: SRRs are technical sequencing runs within their SRX; SRXs are not merged automatically.
- Fastp PE thresholds refer to surviving **fragments**; historical mate-level thresholds must be recalibrated.
- The quantitative BAM contains trusted unique plus Bayesian FDR-resolved fragments; `--skip_em` retains trusted unique only.
- `--count_utr false` avoids UTR passes; `--count_repeats true` retains repeat and flank counts with reference-wide annotation preprocessing.
- OmniNorm size factors are estimated from ALL, then applied unchanged to ALL/INTRON/EXON matrices; `scaling = 1/sizeFactor`.
- `--save_align_intermeds false` (default) suppresses auxiliary full/high-confidence BAM generation and all intermediate BAM publication; final quantitative BAM remains sorted/indexed.
- MultiQC General Statistics counts are **fragments**; native samtools diagnostics remain outside the MultiQC report and are **alignment records**, not fragments.
- Transcript FASTA creation is skipped unless needed for optional Salmon quantification, a sample with `strandedness=auto`, or `--prepare_transcript_fasta true`.
- nf-core modules retained under `modules/nf-core/` are vendored tools and do not require network access to nf-core config repositories.

## Scientific invariants to verify

1. STAR input = STAR unique + multi-locus + too-many-loci + unmapped template counts.
2. Candidate contract: unique + eligible multi + excluded = mapped candidate fragments.
3. Bayesian output: eligible multi = resolved + unresolved.
4. Final quantitative: trusted unique + resolved, or trusted unique only in `--skip_em` mode.
5. Normalized region count = raw region count × the single ALL-derived sample scaling factor.
6. No change in overlap modes, gene models, RepeatMasker coordinate conventions or flank width.

## Key flags

| Parameter | Default | Note |
|---|---|---|
| `--min_trimmed_reads` | 500000 | Fragment threshold in fastp mode |
| `--skip_quantify_salmon` | true | Salmon still used for `auto` strandedness |
| `--count_utr` | false | Avoid additional 5'/3' UTR overlap scans |
| `--count_repeats` | true | Requires compatible RepeatMasker file |
| `--save_align_intermeds` | false | Generate and publish all extra BAMs only when true |
| `--prepare_transcript_fasta` | false | Force transcript FASTA generation |

Any scientific regression must halt adoption of the corresponding step. Read the stage-specific STEP.md files in the v6 delivery archive before applying each replacement.

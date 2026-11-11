# Final BAM and intermediate-alignment output policy

This output contract is applied at the **Nextflow module and publication** levels. It does **not** modify Bayesian posterior inference, the expected-FDR threshold, MAPQ/tag realization or gene/repeat count generation.

## Single published final BAM

For each sample, downstream counting consumes the coordinate-sorted/indexed quantitative BAM. By default, only these two BAM artifacts are published:

```text
outdir/star_salmon/FINAL_BAM/SAMPLE.final.sorted.bam
outdir/star_salmon/FINAL_BAM/SAMPLE.final.sorted.bam.bai
```

With `--bam_csi_index true`, a `.csi` is emitted instead of the usual `.bai` (subject to samtools/index behavior).

Resolution modes:

| Parameters | Quantitative content | Optional diagnostics |
|---|---|---|
| `--skip_em false --save_align_intermeds false` | Trusted unique + FDR-resolved multimappers | Not materialized |
| `--skip_em true --save_align_intermeds false` | Trusted unique only | Not materialized |
| `--skip_em false --save_align_intermeds true` | Trusted unique + FDR-resolved multimappers | Full ambiguity, high-confidence and other intermediate BAMs |
| `--skip_em true --save_align_intermeds true` | Trusted unique only | Full ambiguity and other intermediate BAMs; no redundant high-confidence tier |

**Note:** “final” means final for the quantitative/counting path. It does not mean every mapped fragment is retained. The legacy name `--skip_em` refers to skipping Bayesian resolution, not skipping a conventional EM iteration (none is present in the current resolver).

## File creation versus publication

Input BAMs that are needed for STAR mapping, candidate extraction, fragment filtering, QC or counting may still be created within Nextflow `work/`. The `--save_align_intermeds` switch controls **published intermediate alignment BAMs** and generation of the two *auxiliary realized BAM tiers*. It does not disable necessary temporary data, version tracking or QC.

When the switch is false, `MM_BAM_REALIZE` writes the optional full and high-confidence SAM streams to `/dev/null` while still calculating their **fragment and alignment-record counters** in the realization QC. The process creates only the quantitative BAM. Optional Nextflow outputs allow absent tiers without altering the quantitative-channel contract.

## QC always retained

The main run retains STAR/FASTQ/fastp logs, MultiQC inputs and quantitative BAM native samtools diagnostics independently of BAM publication. Report-level mapping numbers use **fragments**, not doubled paired-end SAM records. Native samtools `.stats`, `.flagstat`, and `.idxstats` are retained under `outdir/qc/samtools_quantitative/` for debugging; they are not treated as fragment counts in MultiQC.

Merged model QC is published in `outdir/qc/multimapper/`. Per-fragment candidate, posterior and decision tables remain available in work directories; intermediate publication of verbose tables follows existing `save_align_intermeds` policy.

## FastQC naming

The current baseline already has explicit FastQC prefixes and destinations:

```text
outdir/fastqc/raw/SAMPLE.raw_1_fastqc.html
outdir/fastqc/trim/SAMPLE.trimmed_1_fastqc.html
```

The two runs are deliberate. One represents pre-trim FASTQs and one post-fastp FASTQs; reporting both is not duplicate execution.

## Regression requirements

Apply to a fresh `--outdir` and preserve the original Nextflow work folder for `-resume`:

1. With intermediate output off: one `.final.sorted.bam` and index per sample in `FINAL_BAM/`; no auxiliary BAM tiers in published output or `MM_BAM_REALIZE` work directory.
2. With intermediate output on: correct auxiliary BAMs appear, with `.full` and `.highconf` containing the documented distinct fragment sets; quantitative BAM does not change except for filename/path.
3. Run both `--skip_em true` and false. In skip mode no selected multimappers should enter final BAM and high-confidence should not be generated.
4. Compare final BAM counts, normalized region/repeat matrices and MultiQC fragment statistics across both publication modes.
5. Require conservation invariants from resolver and realization QC, and check that the same quantitative BAM reaches `INTRON_EXON_COUNT` in all cases.

This patch **does not claim** end-to-end validation on Slurm/Singularity; only the local static and AWK/NumPy known-answer tests have been run.

# mmRNAseq

**RNA-seq processing with fragment-aware counting and conservative Bayesian resolution of genomic multimappers.**

mmRNAseq is a [Nextflow DSL2](https://www.nextflow.io/) workflow that combines STAR alignment, custom multimapper resolution, gene-region and RepeatMasker counting, normalization, and MultiQC. It builds on selected [nf-core/rnaseq](https://github.com/nf-core/rnaseq) components.

## Workflow

```text
FASTQ → FastQC / fastp → STAR → fragment filtering
                                ↓
                     unique anchors + candidate alignments
                                ↓
                       Bayesian resolver (optional)
                                ↓
                        FINAL_BAM + index
                                ↓
               gene / intron / exon / repeat counting
                                ↓
                       normalization + MultiQC
```

## Multimapping: current conservative resolver

**Fixed unique-anchor density → candidate posterior → global expected-FDR selection → quantitative BAM with one selected alignment per accepted fragment.**

For each multimapped fragment, the resolver combines STAR alignment scores with local density estimated from uniquely mapped fragment anchors. It retains the best-supported genomic placement when alignment-quality and model-based expected-FDR criteria are met; otherwise the fragment is excluded from quantitative counting.

**Primary use:** high-confidence genomic visualization, fragment placement, and conventional gene/repeat counting.

This is **not iterative expectation–maximization (EM) abundance estimation**. Ambiguous fragments do not update locus abundance iteratively, and posterior error probabilities are model-dependent rather than independently calibrated. Highly similar loci may remain unresolved.

## Final BAM, counts, and QC

By default, **one coordinate-sorted, indexed `FINAL_BAM` per sample** is published and used downstream:

| Option | Quantitative BAM |
| --- | --- |
| `--skip_em false` (default) | Trusted unique fragments + Bayesian-resolved multimappers |
| `--skip_em true` | Trusted unique fragments only |
| `--save_align_intermeds true` | Also generate/publish intermediate BAMs for inspection |

The `--skip_em` parameter name is historical; the current resolver does not use iterative EM.

Gene quantification includes gene-body (`ALL`), intron, and exon counts, with optional UTR counts (`--count_utr true`). RepeatMasker locus/flank counting is enabled by default (`--count_repeats true`). Different overlap rules mean intron and exon counts need not sum to `ALL`. Normalization factors estimated from the `ALL` matrix are also applied to the intron/exon matrices.

MultiQC prioritizes **fragment-level** mapping and resolver metrics instead of treating paired-end SAM/BAM records as fragments; raw samtools diagnostics remain available separately.

## Quick start

Use matching genome FASTA, GTF, and RepeatMasker annotations, plus a samplesheet with `sample,fastq_1,fastq_2,strandedness` columns (see [`assets/samplesheet.csv`](assets/samplesheet.csv)).

```bash
export NXF_VER=25.10.9

nextflow run ./main.nf \
  --input samplesheet.csv \
  --fasta /refs/genome.fa.gz \
  --gtf /refs/genes.gtf \
  --rmsker /refs/repeatMasker.out.gz \
  --outdir /results/mmrnaseq \
  -profile singularity \
  -w /scratch/mmrnaseq/work \
  --nf_boost_cleanup false \
  -resume
```

Specify a site-specific Nextflow configuration with `-c conf/<site>.config` as needed. Keep the work directory for reliable `-resume`; automatic intermediate cleanup can invalidate cached outputs. To disable repeat counting, set `--count_repeats false`. The model-based multimapper threshold is set with `--mm_resolution_fdr` (default `0.05`).

## Main outputs

```text
outdir/
├── fastqc/raw/                       # Before trimming
├── fastqc/trim/                      # After trimming
├── fastp/                            # Trimming QC
├── star_salmon/FINAL_BAM/            # Final quantitative BAM + index
├── star_salmon/norm_deseq2_intron_exon/  # Normalized gene-region matrices
├── qc/multimapper/                   # Fragment and resolver QC
├── qc/samtools_quantitative/         # Native BAM diagnostics
└── multiqc/                          # Integrated QC report
```

**Limitation:** `FINAL_BAM` represents selected genomic placements, not fractional locus abundances or all alignment alternatives. TE loci with insufficient discriminatory evidence cannot be reliably resolved.

## References

- Dobin *et al.* (2013), [STAR](https://doi.org/10.1093/bioinformatics/bts635). *Bioinformatics*.
- Chen *et al.* (2018), [fastp](https://doi.org/10.1093/bioinformatics/bty560). *Bioinformatics*.
- Li *et al.* (2009), [SAMtools](https://doi.org/10.1093/bioinformatics/btp352). *Bioinformatics*.
- Ewels *et al.* (2016), [MultiQC](https://doi.org/10.1093/bioinformatics/btw354). *Bioinformatics*.
- Ewels *et al.* (2020), [nf-core](https://doi.org/10.1038/s41587-020-0439-x). *Nature Biotechnology*.

The Bayesian resolver is specific to mmRNAseq; these references acknowledge upstream tools rather than independently validating its inference model.

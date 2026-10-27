#!/usr/bin/env python3
"""Static/known-answer contract tests for PATCH 04.

R/Bioconductor is intentionally not required by this test so it can run on the
repository login environment. The real workflow remains the integration test.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def genomic_flanks(start: int, end: int, chrom_len: int, frac: float = 0.30):
    """Python mirror of the PATCH04 flank arithmetic for known-answer tests."""
    width = end - start + 1
    target = max(1, round(width * frac))

    left = None
    left_end = start - 1
    if left_end >= 1:
        left_start = max(1, left_end - target + 1)
        left = (left_start, left_end)

    right = None
    right_start = end + 1
    if right_start <= chrom_len:
        right_end = min(chrom_len, right_start + target - 1)
        right = (right_start, right_end)

    return left, right


def assert_nonoverlap(flank, body):
    if flank is None:
        return
    fs, fe = flank
    bs, be = body
    assert fe < bs or fs > be, (flank, body)


def main():
    r = (ROOT / "bin/intron_exon_count.r").read_text()
    module = (ROOT / "modules/local/intron_exon_count.nf").read_text()
    workflow = (ROOT / "workflows/rnaseq.nf").read_text()

    # RepeatMasker input contract: plain/gzip, empty/malformed validation, and
    # the legitimate optional trailing '*' overlap marker.
    assert 'grepl("\\\\.gz$", path, ignore.case=TRUE)' in r
    assert 'gzip -dc --' in r
    assert 'data.table::fread(path' in r
    assert 'ncol(x) < 15L || ncol(x) > 16L' in r
    assert 'marker != "*"' in r
    assert 'x <- x[, seq_len(15L), with=FALSE]' in r
    assert 'anyNA(rmsk$genoStart)' in r
    assert 'any(rmsk$genoEnd < rmsk$genoStart)' in r

    # A reference/annotation chromosome mismatch must not silently become an
    # all-zero repeat result.
    assert 'No RepeatMasker intervals match TxDb seqlevels' in r
    assert 'Check genome build and chromosome naming' in r

    # Flanks must be non-overlapping intervals with actual exposure reported.
    assert 'left_end <- start(rmksGR) - 1L' in r
    assert 'right_start <- end(rmksGR) + 1L' in r
    assert 'IRanges::IRanges' in r
    assert 'function(x) sum(width(x))' in r
    assert '(dat_rmsk$width*0.3*2)' not in r
    assert 'start(rmksGR_upstream) = start(rmksGR_upstream) - w_width' not in r

    # Known-answer flank tests, including chromosome boundaries.
    left, right = genomic_flanks(101, 200, 1000)
    assert left == (71, 100)
    assert right == (201, 230)
    assert_nonoverlap(left, (101, 200))
    assert_nonoverlap(right, (101, 200))

    left, right = genomic_flanks(1, 100, 1000)
    assert left is None
    assert right == (101, 130)
    assert_nonoverlap(right, (1, 100))

    left, right = genomic_flanks(901, 1000, 1000)
    assert left == (871, 900)
    assert right is None
    assert_nonoverlap(left, (901, 1000))

    # Small elements still get a one-base target flank on each available side.
    left, right = genomic_flanks(10, 10, 20)
    assert left == (9, 9)
    assert right == (11, 11)

    # The module now exposes its existing versions.yml for provenance.
    assert "path 'versions.yml', emit: versions" in module
    assert 'path rmsker' in module
    assert 'intron_exon_count.r' in module

    # PATCH03 contract must survive the rebase unchanged: quantitative is the
    # primary downstream BAM; high-confidence and full ambiguity remain separate.
    assert 'ch_genome_bam      = BAM_FILTER_EM.out.bam_quantitative' in workflow
    assert 'ch_mm_highconf_bam = BAM_FILTER_EM.out.bam_highconf' in workflow
    assert 'ch_mm_full_bam     = BAM_FILTER_EM.out.bam_full' in workflow
    assert 'BAM_SORT_STATS_SAMTOOLS_MM_HIGHCONF' in workflow
    assert 'BAM_SORT_STATS_SAMTOOLS_MM_FULL' in workflow
    assert 'ch_genome_bam = BAM_FILTER_EM.out.bam_confident' not in workflow

    # RepeatMasker counting consumes the main quantitative BAM after sorting and
    # stages the annotation through the reusable path channel.
    assert '.join(ch_genome_bam_index, by: [0])' in workflow
    assert 'ch_rmsker // RepeatMasker .out/.out.gz staged as a Nextflow path' in workflow
    assert 'ch_versions = ch_versions.mix(INTRON_EXON_COUNT.out.versions.first())' in workflow

    print('PASS PATCH04 RepeatMasker hardening rebased on Bayesian BAM tiers')


if __name__ == '__main__':
    main()

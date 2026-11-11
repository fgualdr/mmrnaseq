#!/usr/bin/env python3
"""Static v6 I/O invariants; biological equivalence requires matched HPC runs."""
import hashlib
import json
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]

MM_HASHES = {'bin/mm_bayes_resolver.py': '1dc844631ba662c8ca508a52b9d9978cc8aaef1704c6ba5553abaa9236a2a1dc', 'bin/mm_build_candidates.awk': 'fe5bdc222bbac48c80b2ebdc923243bce5e78647960fa42e682a80f9a9a2705e', 'bin/mm_realize_bam.awk': '0f415f7f167797c321a37cb2f6ffa45b2e37f95d49c2026d6d4a4db8bce9cbed', 'modules/local/mm_bayes_resolver.nf': 'ea9540e98d043c0df0fbaefddcc60ad42c9b0defbff4d8e845322483f4156b95', 'modules/local/mm_skip_resolution.nf': '248dc0c89f78a75a9e6a960006acfafe7ef55780419a43a8075029e92e79653a', 'bin/mm_merge_qc.py': 'e2130264b6fec9b8643d01f72c9ef95cece24717f456f2edffa7bb2dbd897f64'}


def source(p):
    return (ROOT / p).read_text()


def test_multimapper_algorithm_untouched():
    for path, expected in MM_HASHES.items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == expected, path


def test_reference_and_repeat_caching():
    wf = source('workflows/rnaseq.nf')
    r = source('bin/intron_exon_count.r')
    prep = source('bin/prepare_repeat_annotation.r')
    assert 'PREP_REPEAT_ANNOTATION(PREP_GTF.out.txdb_sqlite, ch_rmsker)' in wf
    assert 'readRDS(opt$repeat_annotations)' in r
    assert 'saveRDS(list(rmksGR=rmksGR, rmksGR_updown=rmksGR_updown)' in prep
    assert 'if (isTRUE(opt$count_utr))' in r
    assert 'if (isTRUE(opt$count_repeats))' in r


def test_normalization_outputs_and_single_scale():
    norm = source('bin/intron_exon_count_norm.r')
    nf = source('modules/local/norm_counts_intron_exon.nf')
    assert 'write_region(dat_all, "ALL")' in norm
    assert 'write_region(dat_int, "INTRON")' in norm
    assert 'write_region(dat_ex, "EXON")' in norm
    assert 'sf_used <- 1 / sizeFactors(dds)' in norm
    for tag in ('ALL', 'INTRON', 'EXON'):
        for kind in ('raw', 'norm'):
            assert f'{tag}_READS_{kind}.txt' in nf
    assert 'Final_norm_count_mat.txt' in norm


def test_bam_io_and_skip_mode_contract():
    wf = source('workflows/rnaseq.nf')
    sub = source('subworkflows/local/bam_filter_em.nf')
    module = source('modules/local/mm_candidate_prep.nf')
    bamreal = source('modules/local/mm_bam_realize.nf')
    assert 'BAM_FILTER_EM(ch_genome_bam, PREPARE_GENOME.out.chrom_sizes)' in wf
    assert 'MM_CANDIDATE_PREP.out.namesorted_bam' in sub
    assert 'samtools sort' not in bamreal
    assert 'samtools sort' in module
    assert 'MM_SKIP_RESOLUTION' in sub
    assert 'BAM_SORT_STATS_SAMTOOLS_BAM(ch_genome_bam' in wf
    assert 'params.save_align_intermeds' in wf


def test_fragment_qc_and_multiqc_order():
    wf = source('workflows/rnaseq.nf')
    helper = source('lib/WorkflowRnaseq.groovy')
    conf = source('assets/multiqc_config.yml')
    mod = source('modules/local/multiqc.nf')
    assert 'public static Map getQuantitativeFragmentUsability' in helper
    assert 'trusted_unique + resolved' in helper
    assert 'resolved + unresolved != eligible_multi' in helper
    assert 'mm_fragment_usability_mqc.tsv' in wf
    assert 'mm_fragment_generalstats_mqc.tsv' in wf
    assert 'report_section_order:' in conf
    assert 'module_order:' in conf
    assert '  - samtools' not in conf
    assert "path ('samtools/" not in mod


def test_fastq_fragment_and_reference_contract():
    fastp = source('subworkflows/nf-core/fastq_fastqc_umitools_fastp/main.nf')
    prep = source('subworkflows/local/prepare_genome.nf')
    wf = source('workflows/rnaseq.nf')
    assert 'FASTP (\n            umi_reads' in fastp
    assert 'getFastpFragmentsAfterFiltering' in fastp
    assert 'if (need_transcript_fasta)' in prep
    assert 'hasAutoStrandedness' in wf
    assert 'PREPARE_GENOME.out.transcript_fasta' in wf


def test_parameters_parse():
    schema = json.loads(source('nextflow_schema.json'))
    c = source('nextflow.config')
    for p in ('count_utr','count_repeats','save_align_intermeds','prepare_transcript_fasta'):
        assert p in c and p in json.dumps(schema)

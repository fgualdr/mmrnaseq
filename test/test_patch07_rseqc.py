#!/usr/bin/env python3
from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run(cmd, **kwargs):
    return subprocess.run(cmd, check=True, text=True, capture_output=True, **kwargs)


def read_status(path: Path) -> dict[str, str]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == 1
    return rows[0]


def test_gtf_with_gff_header() -> None:
    # This reproduces the failure mode in the user's T2T annotation: a file can
    # carry a GFF-style header while the records themselves use GTF transcript_id
    # attributes. The historical converter globally switched parsers and emitted
    # an empty RSeQC BED in that situation.
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        gtf = d / 'mixed_header.gtf'
        gtf.write_text(
            '##gff-version 3\n'
            'chr1\ttest\texon\t101\t150\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";\n'
            'chr1\ttest\texon\t201\t250\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";\n'
            'chr2\ttest\texon\t501\t550\t.\t-\t.\tgene_id "g2"; transcript_id "tx2";\n'
        )
        out = run(['perl', str(ROOT / 'bin/gtf2bed'), str(gtf)])
        lines = [x.split('\t') for x in out.stdout.strip().splitlines()]
        assert len(lines) == 2
        tx1 = next(x for x in lines if x[3] == 'tx1')
        assert len(tx1) == 12
        assert tx1[0] == 'chr1'
        assert tx1[1:3] == ['100', '250']
        assert tx1[5] == '+'
        assert tx1[9] == '2'
        assert tx1[10] == '50,50,'
        assert tx1[11] == '0,100,'
        assert 'wrote 2 BED12 transcript models' in out.stderr


def test_gff3_parent() -> None:
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        gff = d / 'a.gff3'
        gff.write_text(
            '##gff-version 3\n'
            'chr1\ttest\texon\t11\t20\t.\t-\t.\tID=ex1;Parent=transcript:txA\n'
            'chr1\ttest\texon\t31\t40\t.\t-\t.\tID=ex2;Parent=transcript:txA\n'
        )
        out = run(['perl', str(ROOT / 'bin/gtf2bed'), str(gff)])
        cols = out.stdout.strip().split('\t')
        assert len(cols) == 12
        assert cols[3] == 'txA'
        assert cols[5] == '-'
        assert cols[9:] == ['2', '10,10,', '0,20,']



def test_fasta_namespace_filter_is_perl_only() -> None:
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        gtf = d / 'a.gtf'
        fasta = d / 'ref.fa'
        gtf.write_text(
            'chr1\ttest\texon\t11\t20\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";\n'
            'chrUn\ttest\texon\t31\t40\t.\t+\t.\tgene_id "g2"; transcript_id "tx2";\n'
        )
        fasta.write_text('>chr1 primary\n' + 'A' * 100 + '\n')
        out = run([
            'perl', str(ROOT / 'bin/gtf2bed'),
            '--fasta', str(fasta), str(gtf),
        ])
        rows = [x.split('\t') for x in out.stdout.strip().splitlines()]
        assert len(rows) == 1
        assert rows[0][0] == 'chr1'
        assert rows[0][3] == 'tx1'
        assert 'reference compatibility kept 1 / 2 transcript models' in out.stderr


def test_zero_fasta_overlap_fails_loudly() -> None:
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        gtf = d / 'a.gtf'
        fasta = d / 'ref.fa'
        gtf.write_text(
            'chr1\ttest\texon\t11\t20\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";\n'
        )
        fasta.write_text('>NC_000001.11\n' + 'A' * 100 + '\n')
        proc = subprocess.run(
            ['perl', str(ROOT / 'bin/gtf2bed'), '--fasta', str(fasta), str(gtf)],
            text=True, capture_output=True,
        )
        assert proc.returncode != 0
        assert 'zero contig-name overlap' in proc.stderr
        assert 'chr1' in proc.stderr
        assert 'NC_000001.11' in proc.stderr

def test_zero_usable_is_indeterminate() -> None:
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        stdout = d / 'infer.txt'
        stderr = d / 'infer.err'
        status = d / 'status.tsv'
        stdout.write_text('')
        stderr.write_text('Total 0 usable reads were sampled\n')
        run([
            sys.executable, str(ROOT / 'bin/rseqc_infer_status.py'),
            '--sample', 'S0', '--provided', 'reverse',
            '--stdout', str(stdout), '--stderr', str(stderr),
            '--cutoff', '30', '--output', str(status),
        ])
        row = read_status(status)
        assert row['Inference status'] == 'NO_USABLE_READS'
        assert row['Check'] == 'INDETERMINATE'
        assert row['Inferred strandedness'] == 'NA'
        assert row['Sense (%)'] == 'NA'
        assert row['Antisense (%)'] == 'NA'
        assert row['Undetermined (%)'] == 'NA'
        assert row['Usable reads'] == '0'
        assert 'zero usable reads' in row['Reason'].lower()


def test_valid_reverse_pairend() -> None:
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        stdout = d / 'infer.txt'
        stderr = d / 'infer.err'
        status = d / 'status.tsv'
        stdout.write_text(
            'This is PairEnd Data\n'
            'Fraction of reads failed to determine: 0.0500\n'
            'Fraction of reads explained by "1++,1--,2+-,2-+": 0.0100\n'
            'Fraction of reads explained by "1+-,1-+,2++,2--": 0.9400\n'
        )
        stderr.write_text('Total 200000 usable reads were sampled\n')
        run([
            sys.executable, str(ROOT / 'bin/rseqc_infer_status.py'),
            '--sample', 'S1', '--provided', 'reverse',
            '--stdout', str(stdout), '--stderr', str(stderr),
            '--cutoff', '30', '--output', str(status),
        ])
        row = read_status(status)
        assert row['Inference status'] == 'VALID'
        assert row['Check'] == 'PASS'
        assert row['Inferred strandedness'] == 'reverse'
        assert abs(float(row['Sense (%)']) - 1.0) < 1e-8
        assert abs(float(row['Antisense (%)']) - 94.0) < 1e-8
        assert abs(float(row['Undetermined (%)']) - 5.0) < 1e-8
        assert row['Usable reads'] == '200000'


def test_malformed_fractions_are_indeterminate() -> None:
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        stdout = d / 'infer.txt'
        stderr = d / 'infer.err'
        status = d / 'status.tsv'
        stdout.write_text(
            'This is PairEnd Data\n'
            'Fraction of reads failed to determine: 0.1\n'
            'Fraction of reads explained by "1++,1--,2+-,2-+": 0.1\n'
        )
        stderr.write_text('Total 200000 usable reads were sampled\n')
        run([
            sys.executable, str(ROOT / 'bin/rseqc_infer_status.py'),
            '--sample', 'S2', '--provided', 'reverse',
            '--stdout', str(stdout), '--stderr', str(stderr),
            '--output', str(status),
        ])
        row = read_status(status)
        assert row['Inference status'] == 'INVALID_OUTPUT'
        assert row['Check'] == 'INDETERMINATE'
        assert row['Inferred strandedness'] == 'NA'


def test_workflow_plumbing() -> None:
    gtf_module = (ROOT / 'modules/local/gtf2bed.nf').read_text()
    prepare = (ROOT / 'subworkflows/local/prepare_genome.nf').read_text()
    infer = (ROOT / 'modules/nf-core/rseqc/inferexperiment/main.nf').read_text()
    rseqc_sub = (ROOT / 'subworkflows/nf-core/bam_rseqc/main.nf').read_text()
    workflow = (ROOT / 'workflows/rnaseq.nf').read_text()
    helper = (ROOT / 'lib/WorkflowRnaseq.groovy').read_text()

    assert 'path fasta' in gtf_module
    assert '--fasta $fasta' in gtf_module
    assert 'comm -12' not in gtf_module
    assert 'compatible.contigs' not in gtf_module
    assert 'Do not add coreutils-only helpers' in gtf_module
    assert 'GTF2BED ( ch_gtf, ch_fasta )' in prepare

    assert '*.infer_experiment.stderr.txt' in infer
    assert 'emit: diagnostic' in infer
    assert 'emit: qc_status' in infer
    assert 'emit: stderr' not in infer
    assert '*.infer_experiment.status.tsv' in infer
    assert 'rseqc_infer_status.py' in infer
    assert 'inferexperiment_diagnostic = RSEQC_INFEREXPERIMENT.out.diagnostic' in rseqc_sub
    assert 'RSEQC_INFEREXPERIMENT.out.stderr' not in rseqc_sub
    assert 'inferexperiment_status = RSEQC_INFEREXPERIMENT.out.qc_status' in rseqc_sub
    assert 'BAM_RSEQC.out.inferexperiment_status' in workflow
    assert 'getRseqcStrandStatus' in workflow
    assert "pass_strand_check[meta.id] = 'INDETERMINATE'" in workflow
    assert "strand_qc_mqc.tsv" in workflow
    assert 'getInferexperimentStrandedness' not in helper


def main() -> None:
    test_gtf_with_gff_header()
    test_gff3_parent()
    test_fasta_namespace_filter_is_perl_only()
    test_zero_fasta_overlap_fails_loudly()
    test_zero_usable_is_indeterminate()
    test_valid_reverse_pairend()
    test_malformed_fractions_are_indeterminate()
    test_workflow_plumbing()
    print('PASS PATCH07 corrected RSeQC BED12 + fail-closed strandedness QC')


if __name__ == '__main__':
    main()

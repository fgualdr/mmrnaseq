#!/usr/bin/env python3
"""Validate fragment QC stability when auxiliary BAM materialization is disabled.

The Nextflow container/runtime wiring additionally needs an HPC test.
"""
import importlib.util
import json
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load_fixtures():
    spec = importlib.util.spec_from_file_location('bam_fixtures', ROOT / 'test/test_patch03_bam_realization.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.RES, module.SAM


def test_auxiliary_output_switch_preserves_quantitative_and_qc():
    resolution, sam = load_fixtures()
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / 'decisions.tsv').write_text(resolution)
        (root / 'reads.sam').write_text(sam)
        results = {}
        for mode in ('all', 'quantitative_only'):
            d = root / mode
            d.mkdir()
            full = str(d / 'full.sam') if mode == 'all' else '/dev/null'
            highconf = str(d / 'high.sam') if mode == 'all' else '/dev/null'
            subprocess.run([
                'awk', '-v', f'full_out={full}',
                '-v', f'quantitative_out={d / "quant.sam"}',
                '-v', f'highconf_out={highconf}',
                '-v', 'highconf_min=0.95',
                '-v', f'qc_out={d / "qc.tsv"}',
                '-f', str(ROOT / 'bin/mm_realize_bam.awk'),
                str(root / 'decisions.tsv'), str(root / 'reads.sam'),
            ], check=True)
            results[mode] = ((d / 'quant.sam').read_bytes(), (d / 'qc.tsv').read_bytes())
            if mode == 'quantitative_only':
                assert not (d / 'full.sam').exists()
                assert not (d / 'high.sam').exists()
        assert results['all'] == results['quantitative_only']


def test_one_final_bam_and_optional_aux_contract():
    conf = (ROOT / 'conf/modules.config').read_text()
    module = (ROOT / 'modules/local/mm_bam_realize.nf').read_text()
    wf = (ROOT / 'workflows/rnaseq.nf').read_text()
    schema = json.loads((ROOT / 'nextflow_schema.json').read_text())
    assert '${params.outdir}/${params.aligner}/FINAL_BAM' in conf
    assert '${meta.id}.final.sorted' in conf
    assert 'fastqc/raw' in conf and 'fastqc/trim' in conf
    assert "def fullRequested = params.save_align_intermeds" in module
    assert "def highconfRequested = params.save_align_intermeds && !params.skip_em" in module
    assert 'optional: true, emit: full_bam' in module
    assert 'optional: true, emit: highconf_bam' in module
    assert 'if (params.save_align_intermeds)' in wf
    assert 'BAM_SORT_STATS_SAMTOOLS_BAM(ch_genome_bam, PREPARE_GENOME.out.fasta)' in wf
    assert 'save_mm_aux_bams' not in json.dumps(schema)
    assert 'save_mm_aux_bams' not in wf

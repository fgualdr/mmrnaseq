#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    multiqc_nf = (ROOT / 'modules/local/multiqc.nf').read_text()
    dump_nf = (ROOT / 'modules/nf-core/custom/dumpsoftwareversions/main.nf').read_text()
    resolver_nf = (ROOT / 'modules/local/mm_bayes_resolver.nf').read_text()
    workflow = (ROOT / 'workflows/rnaseq.nf').read_text()
    helper = (ROOT / 'lib/WorkflowRnaseq.groovy').read_text()
    sub = (ROOT / 'subworkflows/local/bam_filter_em.nf').read_text()
    merge_nf = (ROOT / 'modules/local/mm_qc_merge.nf').read_text()
    merge_py = (ROOT / 'bin/mm_merge_qc.py').read_text()
    template = (ROOT / 'lib/NfcoreTemplate.groovy').read_text()
    schema = json.loads((ROOT / 'nextflow_schema.json').read_text())

    # MultiQC itself and the software-version aggregator are upgraded to 1.35.
    for text in (multiqc_nf, dump_nf):
        assert 'multiqc=1.35' in text
        assert 'multiqc:1.35--pyhdfd78af_2' in text
        assert 'multiqc=1.14' not in text

    # Do NOT silently change the already validated resolver runtime just because it
    # happens to use an older MultiQC container as its Python/NumPy environment.
    assert 'multiqc=1.14' in resolver_nf
    assert 'multiqc=1.35' not in resolver_nf

    # Salmon QC must actually reach MultiQC and must not be reset to an empty channel.
    assert 'ch_salmon_multiqc = QUANTIFY_STAR_SALMON.out.results' in workflow
    assert len(re.findall(r'ch_salmon_multiqc\s*=\s*Channel\.empty\(\)', workflow)) == 1

    # Mapping acceptance moved after BAM_FILTER_EM and uses the final quantitative
    # fragment contract, not STAR unique mapping alone.
    assert 'getQuantitativeFragmentUsability' in helper
    assert 'getStarPercentMapped' not in helper
    assert 'quantitative_usable = trusted_unique + resolved' in helper
    assert "params.skip_em ? 'UNIQUE_ONLY' : 'BAYES_FDR'" in helper
    assert '--skip_em true must not resolve multimappers' in helper
    assert workflow.index('BAM_FILTER_EM(ch_bam_genome_bai') < workflow.index('getQuantitativeFragmentUsability')
    assert 'Quantitative usable fragments (%)' in workflow
    assert 'BAM_FILTER_EM.out.bam_quantitative' in workflow
    assert '.join(ch_mapping_pass, by: [0])' in workflow

    # Preseq now measures the same quantitative BAM used downstream.
    assert workflow.index('BAM_SORT_STATS_SAMTOOLS_BAM(ch_genome_bam') < workflow.index('PRESEQ_LCEXTRAP(ch_genome_bam)')

    # PATCH06 unique-only semantics survive: the subworkflow is still always called,
    # and skip_em changes the decision source rather than passing STAR multimappers.
    assert 'MM_SKIP_RESOLUTION' in sub
    assert 'ch_quantitative_bam = BAM_FILTER.out.bam' not in sub
    assert 'mm_realization_qc = MM_BAM_REALIZE.out.qc' in sub

    # PATCH05 k diagnostics are surfaced as a dedicated MultiQC table without
    # changing the current k search or resolver inference.
    assert "path('*.mmdensity_mqc.tsv'), emit: density_multiqc" in merge_nf
    assert 'mm_density_multiqc' in sub
    assert 'mmdensity_mqc.tsv' in merge_py
    assert '# id: mmdensity' in merge_py
    assert 'Gain vs previous k' in merge_py
    assert 'P95 CV radius (bp)' in merge_py

    # MultiQC process inputs and workflow call must stay positional and aligned.
    input_block = multiqc_nf.split('    input:\n', 1)[1].split('    output:\n', 1)[0]
    n_inputs = len(re.findall(r'^\s+path\b', input_block, flags=re.M))
    call_block = workflow.split('        MULTIQC (\n', 1)[1].split('        )\n', 1)[0]
    call_args = [ln for ln in call_block.splitlines() if ln.strip()]
    assert n_inputs == len(call_args), (n_inputs, len(call_args))
    for token in ['ch_fragment_usability_multiqc', 'ch_mmresolver_multiqc', 'ch_mmdensity_multiqc', 'ch_salmon_multiqc']:
        assert token in call_block, token

    # Cleaned workflow must not reintroduce dead aligners into the report.
    cfg_text = (ROOT / 'assets/multiqc_config.yml').read_text()
    assert '\n  - hisat2\n' not in cfg_text
    assert '\n  - rsem\n' not in cfg_text
    for token in ['mm_fragment_usability', 'mmresolver', 'mmdensity']:
        assert token in cfg_text
    if yaml is not None:
        cfg = yaml.safe_load(cfg_text)
        assert {'mm_fragment_usability', 'mmresolver', 'mmdensity'} <= set(cfg['custom_data'])

    # Schema / completion messaging now describe the actual threshold semantics.
    blob = json.dumps(schema)
    assert 'final quantitative BAM after multimapper handling' in blob
    assert 'trusted unique plus globally expected-FDR-resolved multimappers' in blob
    assert 'quantitative-usable-fragment threshold' in template
    assert 'pass_mapped_reads, pass_trimmed_reads, pass_strand_check' in workflow

    # Known-answer check that PATCH05 merge now produces the dedicated density table.
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        (d / 'contract.tsv').write_text(
            'metric\tvalue\n'
            'total_fragments\t100\nunique_fragments\t70\nmultimapper_fragments\t25\n'
            'excluded_fragments\t5\ncandidate_rows\t60\nanchor_rows\t140\n'
        )
        (d / 'resolver.tsv').write_text(
            'metric\tvalue\nmultimapper_fragments\t25\ncandidate_rows\t60\nresolved_fragments\t20\n'
            'unresolved_fragments\t5\nabsolute_score_failures\t1\ntarget_posterior_fdr\t0.05\n'
            'achieved_expected_fdr\t0.03\nequal_best_as_fraction\t0.4\nmedian_candidate_count\t2\n'
            'p95_candidate_count\t5\nmedian_pmax\t0.97\nmedian_entropy\t0.1\nmedian_as_delta\t0\n'
        )
        (d / 'realization.tsv').write_text(
            'metric\tvalue\nunique_fragments\t70\nresolved_fdr_fragments\t20\n'
            'resolved_highconf_fragments\t15\nresolved_lower_confidence_fragments\t5\n'
            'unresolved_fragments\t5\nexcluded_fragments\t5\nhighconf_min_pmax\t0.95\n'
            'full_bam_alignment_records\t220\nquantitative_bam_alignment_records\t180\n'
            'highconf_bam_alignment_records\t170\n'
        )
        (d / 'density.tsv').write_text(
            'orientation_group\tunique_anchor_weight\tbackground_density_per_bp\tselected_knn_mass\t'
            'knn_search_max\tknn_boundary_hit\tknn_cv_loglik\tunique_normalized_as_median\t'
            'unique_normalized_as_mad\tabsolute_normalized_as_cutoff\n'
            'fw\t140\t0.001\t256\t256\t1\t-2.0\t0.98\t0.01\t0.95\n'
        )
        (d / 'density_cv.tsv').write_text(
            'orientation_group\tcandidate_knn_mass\tknn_cv_loglik\tmedian_cv_radius_bp\t'
            'p95_cv_radius_bp\tsampled_anchor_positions\tselected\tsearch_boundary\n'
            'fw\t128\t-2.1\t100\t200\t1000\t0\t0\n'
            'fw\t256\t-2.0\t200\t400\t1000\t1\t1\n'
        )
        subprocess.run([
            sys.executable, str(ROOT / 'bin/mm_merge_qc.py'), '--sample', 'S',
            '--contract', str(d / 'contract.tsv'), '--resolver', str(d / 'resolver.tsv'),
            '--realization', str(d / 'realization.tsv'), '--density', str(d / 'density.tsv'),
            '--density-cv', str(d / 'density_cv.tsv'), '--output-prefix', str(d / 'S'),
        ], check=True)
        with (d / 'S.mmdensity_mqc.tsv').open() as handle:
            data_lines = [line for line in handle if not line.startswith('#')]
        rows = list(csv.DictReader(data_lines, delimiter='\t'))
        assert len(rows) == 1
        assert rows[0]['Sample-orientation'] == 'S|fw'
        assert rows[0]['Selected k'] == '256'
        assert rows[0]['Boundary hit'] == '1'
        assert rows[0]['P95 CV radius (bp)'] == '400'

    print('PASS PATCH07 MultiQC 1.35 + quantitative-fragment / resolver / k QC')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Synthetic skip-mode decision table: all multimappers are deliberately unresolved;
# excluded fragment remains excluded; trusted unique fragment has no decision row.
DECISIONS = '''fragment_id\tstatus\tselected_candidate_id\tselected_hi\tnh\tcandidate_count\torientation_group\tpmax\tposterior_error\tentropy\tbest_as\tsecond_as\tas_delta\ttied_best\tnormalized_as\tabsolute_score_pass
qM\tUNRESOLVED\t.\t0\t2\t2\t.\t0\t1\t0.693147\t100\t100\t0\t2\t1\t0
qX\tEXCLUDED\t.\t0\t1\t0\t.\t0\t1\t0\t0\t0\t0\t0\t0\t0
'''

SAM = '''@HD\tVN:1.6\tSO:queryname
qM\t99\tchr1\t101\t3\t50M\t=\t301\t250\tAAAA\t*\tNH:i:2\tHI:i:1\tAS:i:100
qM\t147\tchr1\t301\t3\t50M\t=\t101\t-250\tAAAA\t*\tNH:i:2\tHI:i:1\tAS:i:100
qM\t355\tchr2\t501\t3\t50M\t=\t701\t250\tAAAA\t*\tNH:i:2\tHI:i:2\tAS:i:100
qM\t403\tchr2\t701\t3\t50M\t=\t501\t-250\tAAAA\t*\tNH:i:2\tHI:i:2\tAS:i:100
qX\t2048\tchr6\t101\t20\t50M\t*\t0\t0\tAAAA\t*\tNH:i:1\tHI:i:1\tAS:i:100
qU\t99\tchr5\t101\t255\t50M\t=\t301\t250\tAAAA\t*\tNH:i:1\tHI:i:1\tAS:i:100
qU\t147\tchr5\t301\t255\t50M\t=\t101\t-250\tAAAA\t*\tNH:i:1\tHI:i:1\tAS:i:100
'''


def qrecords(text: str, qname: str) -> list[str]:
    return [x for x in text.splitlines() if x.startswith(qname + "\t")]


def qc_map(text: str) -> dict[str, str]:
    rows = [x.split("\t") for x in text.strip().splitlines()]
    assert rows[0] == ["metric", "value"]
    return dict(rows[1:])


def main() -> None:
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        (d / 'decisions.tsv').write_text(DECISIONS)
        (d / 'in.sam').write_text(SAM)

        subprocess.run([
            'awk',
            '-v', f'full_out={d}/full.sam',
            '-v', f'quantitative_out={d}/quant.sam',
            '-v', f'highconf_out={d}/high.sam',
            '-v', 'highconf_min=0.95',
            '-v', f'qc_out={d}/qc.tsv',
            '-f', str(ROOT / 'bin/mm_realize_bam.awk'),
            str(d / 'decisions.tsv'),
            str(d / 'in.sam'),
        ], check=True)

        full = (d / 'full.sam').read_text()
        quant = (d / 'quant.sam').read_text()
        high = (d / 'high.sam').read_text()
        qc = qc_map((d / 'qc.tsv').read_text())

        # --skip_em=true contract: no multimapper may enter the quantitative BAM,
        # even though STAR has a primary HI=1 alignment for qM.
        assert len(qrecords(quant, 'qM')) == 0
        assert len(qrecords(high, 'qM')) == 0
        assert len(qrecords(full, 'qM')) == 4
        assert all('ZR:Z:UNRESOLVED' in r for r in qrecords(full, 'qM'))

        # Fragment-level exclusions are not rescued merely because NH=1.
        assert len(qrecords(quant, 'qX')) == 0
        assert len(qrecords(high, 'qX')) == 0
        assert len(qrecords(full, 'qX')) == 1
        assert 'ZR:Z:EXCLUDED' in qrecords(full, 'qX')[0]

        # A trusted unique paired-end fragment remains whole in the quantitative set.
        assert len(qrecords(quant, 'qU')) == 2
        assert len(qrecords(high, 'qU')) == 2
        assert len(qrecords(full, 'qU')) == 2
        assert all('NH:i:1' in r for r in qrecords(quant, 'qU'))

        assert qc['unique_fragments'] == '1'
        assert qc['resolved_fdr_fragments'] == '0'
        assert qc['resolved_highconf_fragments'] == '0'
        assert qc['resolved_lower_confidence_fragments'] == '0'
        assert qc['unresolved_fragments'] == '1'
        assert qc['excluded_fragments'] == '1'
        assert qc['quantitative_bam_alignment_records'] == '2'
        assert qc['highconf_bam_alignment_records'] == '2'

    # Static integration: BAM_FILTER_EM is always called, while the subworkflow
    # switches only the decision source (Bayesian vs all-multimappers-unresolved).
    wf = (ROOT / 'workflows/rnaseq.nf').read_text()
    sub = (ROOT / 'subworkflows/local/bam_filter_em.nf').read_text()
    skipmod = (ROOT / 'modules/local/mm_skip_resolution.nf').read_text()

    assert 'BAM_FILTER_EM(ch_bam_genome_bai, PREPARE_GENOME.out.chrom_sizes)' in wf
    assert '--skip_em true : trusted unique fragments ONLY' in wf
    assert 'ch_quantitative_bam = BAM_FILTER.out.bam' not in sub
    assert 'MM_SKIP_RESOLUTION' in sub
    assert 'ch_resolution_for_realization = MM_SKIP_RESOLUTION.out.resolution' in sub
    assert 'ch_resolution_for_realization = MM_BAYES_RESOLVER.out.resolution' in sub
    assert 'MM_BAM_REALIZE(ch_realize_inputs)' in sub
    assert 'mm_realization_qc = MM_BAM_REALIZE.out.qc' in sub
    assert '"UNRESOLVED"' in skipmod

    # Picard must not run on the upstream STAR BAM in either mode.
    assert 'This applies to both --skip_em false and --skip_em true' in wf

    print('PASS PATCH06 skip_em=true trusted-unique-only quantitative BAM')


if __name__ == '__main__':
    main()

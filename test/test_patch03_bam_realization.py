#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

RES = '''fragment_id\tstatus\tselected_candidate_id\tselected_hi\tnh\tcandidate_count\torientation_group\tpmax\tposterior_error\tentropy\tbest_as\tsecond_as\tas_delta\ttied_best\tnormalized_as\tabsolute_score_pass
qRH\tRESOLVED\tqRH|HI=2\t2\t2\t2\t.\t0.99\t0.01\t0.05\t100\t90\t10\t1\t1\t1
qRL\tRESOLVED\tqRL|HI=1\t1\t2\t2\t.\t0.80\t0.20\t0.50\t100\t100\t0\t2\t1\t1
qU\tUNRESOLVED\tqU|HI=1\t1\t2\t2\t.\t0.55\t0.45\t0.69\t100\t100\t0\t2\t1\t1
qX\tEXCLUDED\t.\t0\t1\t0\t.\t0\t1\t0\t0\t0\t0\t0\t0\t0
'''

SAM = '''@HD\tVN:1.6\tSO:queryname
qRH\t99\tchr1\t101\t3\t50M\t=\t301\t250\tAAAA\t*\tNH:i:2\tHI:i:1\tAS:i:90
qRH\t147\tchr1\t301\t3\t50M\t=\t101\t-250\tAAAA\t*\tNH:i:2\tHI:i:1\tAS:i:90
qRH\t355\tchr2\t501\t3\t50M\t=\t701\t250\tAAAA\t*\tNH:i:2\tHI:i:2\tAS:i:100
qRH\t403\tchr2\t701\t3\t50M\t=\t501\t-250\tAAAA\t*\tNH:i:2\tHI:i:2\tAS:i:100
qRL\t0\tchr7\t101\t3\t50M\t*\t0\t0\tAAAA\t*\tNH:i:2\tHI:i:1\tAS:i:100
qRL\t256\tchr8\t101\t3\t50M\t*\t0\t0\tAAAA\t*\tNH:i:2\tHI:i:2\tAS:i:100
qU\t0\tchr3\t101\t3\t50M\t*\t0\t0\tAAAA\t*\tNH:i:2\tHI:i:1\tAS:i:100
qU\t256\tchr4\t101\t3\t50M\t*\t0\t0\tAAAA\t*\tNH:i:2\tHI:i:2\tAS:i:100
qX\t2048\tchr6\t101\t20\t50M\t*\t0\t0\tAAAA\t*\tNH:i:1\tHI:i:1\tAS:i:100
q1\t0\tchr5\t101\t255\t50M\t*\t0\t0\tAAAA\t*\tNH:i:1\tHI:i:1\tAS:i:100
'''


def parse_qc(text: str) -> dict[str, str]:
    rows = [line.split("\t") for line in text.strip().splitlines()]
    assert rows[0] == ["metric", "value"]
    return {k: v for k, v in rows[1:]}


def records(text: str, qname: str) -> list[str]:
    return [line for line in text.splitlines() if line.startswith(qname + "\t")]


def main() -> None:
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        (d / "r.tsv").write_text(RES)
        (d / "s.sam").write_text(SAM)

        subprocess.run(
            [
                "awk",
                "-v", f"full_out={d}/full.sam",
                "-v", f"quantitative_out={d}/quant.sam",
                "-v", f"highconf_out={d}/high.sam",
                "-v", "highconf_min=0.95",
                "-v", f"qc_out={d}/qc.tsv",
                "-f", str(ROOT / "bin/mm_realize_bam.awk"),
                str(d / "r.tsv"),
                str(d / "s.sam"),
            ],
            check=True,
        )

        full = (d / "full.sam").read_text()
        quant = (d / "quant.sam").read_text()
        high = (d / "high.sam").read_text()
        qc = parse_qc((d / "qc.tsv").read_text())

        # High-posterior resolved PE fragment: selected HI=2 only, both mates kept,
        # secondary bit cleared (355->99 and 403->147), present in all tiers.
        assert len(records(full, "qRH")) == 2
        assert len(records(quant, "qRH")) == 2
        assert len(records(high, "qRH")) == 2
        assert any(line.split("\t")[1] == "99" for line in records(quant, "qRH"))
        assert any(line.split("\t")[1] == "147" for line in records(quant, "qRH"))
        assert all("NH:i:1" in line and "HI:i:1" in line for line in records(quant, "qRH"))
        assert all("ZR:Z:RESOLVED_FDR" in line and "ZP:f:0.99" in line for line in records(quant, "qRH"))
        assert all("ZR:Z:RESOLVED_HIGHCONF" in line for line in records(high, "qRH"))

        # Globally FDR-resolved but individually lower-confidence fragment remains in
        # the quantitative BAM but is intentionally absent from the high-confidence BAM.
        assert len(records(full, "qRL")) == 1
        assert len(records(quant, "qRL")) == 1
        assert len(records(high, "qRL")) == 0
        assert "ZR:Z:RESOLVED_FDR" in records(quant, "qRL")[0]
        assert "ZP:f:0.8" in records(quant, "qRL")[0]

        # Unresolved alternatives are all retained only in full. They are labelled,
        # but no ZP is attached because fragment pmax is not candidate-specific.
        assert len(records(full, "qU")) == 2
        assert len(records(quant, "qU")) == 0
        assert len(records(high, "qU")) == 0
        assert all("ZR:Z:UNRESOLVED" in line for line in records(full, "qU"))
        assert all("ZP:f:" not in line for line in records(full, "qU"))

        # Excluded chimeric/supplementary-style fragment is archival only even if NH=1.
        assert len(records(full, "qX")) == 1
        assert len(records(quant, "qX")) == 0
        assert len(records(high, "qX")) == 0
        assert "ZR:Z:EXCLUDED" in records(full, "qX")[0]

        # Ordinary trusted unique fragment is unchanged and present in every tier.
        assert len(records(full, "q1")) == 1
        assert len(records(quant, "q1")) == 1
        assert len(records(high, "q1")) == 1

        assert qc["unique_fragments"] == "1"
        assert qc["resolved_fdr_fragments"] == "2"
        assert qc["resolved_highconf_fragments"] == "1"
        assert qc["resolved_lower_confidence_fragments"] == "1"
        assert qc["unresolved_fragments"] == "1"
        assert qc["excluded_fragments"] == "1"
        assert abs(float(qc["highconf_min_pmax"]) - 0.95) < 1e-12

        # Static integration checks: quantitative tier must feed the main workflow;
        # stricter high-confidence and full tiers are separate products.
        subwf = (ROOT / "subworkflows/local/bam_filter_em.nf").read_text()
        wf = (ROOT / "workflows/rnaseq.nf").read_text()
        assert "bam_quantitative = ch_quantitative_bam" in subwf
        assert "bam_highconf = ch_highconf_bam" in subwf
        assert "BAM_FILTER_EM.out.bam_quantitative" in wf
        assert "BAM_FILTER_EM.out.bam_highconf" in wf
        assert "MMRESOLVE_BUILD" not in subwf

    print("PASS PATCH03 Bayesian BAM tier realization")


if __name__ == "__main__":
    main()

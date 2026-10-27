#!/usr/bin/env python3
from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESOLVER = ROOT / "bin" / "mm_bayes_resolver.py"


def rows(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        (d / "chrom").write_text("chr1\t100000\n")
        (d / "uniq").write_text(
            "fragment_id\torientation_group\tas\tnM\tnm\tquery_bases\taligned_query_bases\tnormalized_as\tpaired\n"
            + "".join(f"u{i}\t.\t100\t0\t0\t100\t100\t1\t0\n" for i in range(8))
        )
        (d / "anchors").write_text(
            "fragment_id\tchrom\tpos0\torientation_group\tweight\tmate\tblock_index\tquery_bases\n"
            + "".join(f"u{i}\tchr1\t{1000+i*10}\t.\t1\tSE\t1\t100\n" for i in range(8))
        )
        hdr = (
            "fragment_id\tcandidate_id\thi\tnh\tpaired\torientation_group\tas\tnM\tnm\tquery_bases\t"
            "aligned_query_bases\tnormalized_as\tr1_chrom\tr1_pos1\tr1_flag\tr1_cigar\tr1_blocks\t"
            "r2_chrom\tr2_pos1\tr2_flag\tr2_cigar\tr2_blocks\tcandidate_signature\n"
        )
        body = (
            "f1\tf1|HI=1\t1\t2\t0\t.\t100\t0\t0\t100\t100\t1\tchr1\t1001\t0\t100M\tchr1:1000-1100:1\t.\t0\t0\t.\t.\tA\n"
            "f1\tf1|HI=2\t2\t2\t0\t.\t100\t0\t0\t100\t100\t1\tchr1\t9001\t0\t100M\tchr1:9000-9100:1\t.\t0\t0\t.\t.\tB\n"
            "f2\tf2|HI=1\t1\t2\t0\t.\t100\t0\t0\t100\t100\t1\tchr1\t1001\t0\t100M\tchr1:1000-1100:1\t.\t0\t0\t.\t.\tA\n"
            "f2\tf2|HI=2\t2\t2\t0\t.\t96\t1\t1\t100\t100\t0.96\tchr1\t9001\t0\t100M\tchr1:9000-9100:1\t.\t0\t0\t.\t.\tB\n"
        )
        (d / "cand").write_text(hdr + body)

        subprocess.run(
            [
                sys.executable,
                str(RESOLVER),
                "--candidates", str(d / "cand"),
                "--unique-fragments", str(d / "uniq"),
                "--anchors", str(d / "anchors"),
                "--chrom-sizes", str(d / "chrom"),
                "--resolution-fdr", "0.05",
                "--output-prefix", str(d / "x"),
                "--chunk-rows", "2",  # deliberately tiny: checks fragment-safe chunking
            ],
            check=True,
        )

        post = rows(d / "x.posterior_alignment.tsv")
        by = {}
        for row in post:
            by.setdefault(row["fragment_id"], []).append(row)
        assert abs(sum(float(x["posterior"]) for x in by["f1"]) - 1.0) < 1e-9
        assert float(by["f1"][0]["posterior"]) > float(by["f1"][1]["posterior"])
        assert float(by["f2"][0]["posterior"]) > float(by["f2"][1]["posterior"])

        resolution = rows(d / "x.resolution.tsv")
        assert len(resolution) == 2
        assert {r["fragment_id"] for r in resolution} == {"f1", "f2"}

        density = rows(d / "x.density_model.tsv")
        assert density and density[0]["orientation_group"] == "."
        assert "knn_search_max" in density[0]
        assert "knn_boundary_hit" in density[0]

        density_cv = rows(d / "x.density_cv.tsv")
        assert density_cv
        assert sum(int(r["selected"]) for r in density_cv) == 1
        assert sum(int(r["search_boundary"]) for r in density_cv) == 1
        assert all(float(r["median_cv_radius_bp"]) >= 0 for r in density_cv)
        assert all(float(r["p95_cv_radius_bp"]) >= float(r["median_cv_radius_bp"]) for r in density_cv)

        qc = {r["metric"]: r["value"] for r in rows(d / "x.resolver_qc.tsv")}
        assert qc["multimapper_fragments"] == "2"
        assert qc["candidate_rows"] == "4"

    print("PASS PATCH02 Python/NumPy STAR-candidate Bayesian resolver")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MERGE = ROOT / "bin" / "mm_merge_qc.py"


def read_long(path: Path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {r["metric"]: r["value"] for r in rows}


def main() -> None:
    resolver_py = (ROOT / "bin/mm_bayes_resolver.py").read_text()
    resolver_nf = (ROOT / "modules/local/mm_bayes_resolver.nf").read_text()
    subworkflow = (ROOT / "subworkflows/local/bam_filter_em.nf").read_text()
    config = (ROOT / "conf/modules.config").read_text()

    # PATCH05 is diagnostic only: retain the validated PATCH02P search grid while
    # exposing every tested point and explicit boundary/physical-scale diagnostics.
    assert "while k <= 256.0 and k < total" in resolver_py
    assert "write_density_cv" in resolver_py
    assert "median_cv_radius_bp" in resolver_py
    assert "p95_cv_radius_bp" in resolver_py
    assert "knn_search_max" in resolver_py
    assert "knn_boundary_hit" in resolver_py
    assert 'path("*.density_cv.tsv"), emit: density_cv' in resolver_nf

    # PATCH03 BAM semantics must remain untouched by the QC rebase.
    assert "bam = ch_quantitative_bam" in subworkflow
    assert "bam_quantitative = ch_quantitative_bam" in subworkflow
    assert "bam_highconf = ch_highconf_bam" in subworkflow
    assert "bam_full = ch_full_bam" in subworkflow
    assert "MM_QC_MERGE" in subworkflow
    assert ".join(MM_BAYES_RESOLVER.out.density_cv, by: [0])" in subworkflow
    assert "mm_density_cv" in subworkflow
    assert "mm_multiqc" in subworkflow
    assert "Final_BAM/mm_qc" in config

    with tempfile.TemporaryDirectory() as td:
        d = Path(td)
        (d / "contract.tsv").write_text(
            "metric\tvalue\n"
            "total_fragments\t100\n"
            "unique_fragments\t70\n"
            "multimapper_fragments\t25\n"
            "excluded_fragments\t5\n"
            "candidate_rows\t60\n"
            "anchor_rows\t140\n"
        )
        (d / "resolver.tsv").write_text(
            "metric\tvalue\n"
            "multimapper_fragments\t25\n"
            "candidate_rows\t60\n"
            "resolved_fragments\t20\n"
            "unresolved_fragments\t5\n"
            "absolute_score_failures\t1\n"
            "target_posterior_fdr\t0.05\n"
            "achieved_expected_fdr\t0.03\n"
            "equal_best_as_fraction\t0.4\n"
            "median_candidate_count\t2\n"
            "p95_candidate_count\t5\n"
            "median_pmax\t0.97\n"
            "median_entropy\t0.1\n"
            "median_as_delta\t0\n"
        )
        (d / "realization.tsv").write_text(
            "metric\tvalue\n"
            "unique_fragments\t70\n"
            "resolved_fdr_fragments\t20\n"
            "resolved_highconf_fragments\t15\n"
            "resolved_lower_confidence_fragments\t5\n"
            "unresolved_fragments\t5\n"
            "excluded_fragments\t5\n"
            "highconf_min_pmax\t0.95\n"
            "full_bam_alignment_records\t220\n"
            "quantitative_bam_alignment_records\t180\n"
            "highconf_bam_alignment_records\t170\n"
        )
        (d / "density.tsv").write_text(
            "orientation_group\tunique_anchor_weight\tbackground_density_per_bp\tselected_knn_mass\t"
            "knn_search_max\tknn_boundary_hit\tknn_cv_loglik\tunique_normalized_as_median\t"
            "unique_normalized_as_mad\tabsolute_normalized_as_cutoff\n"
            "fw\t140\t0.001\t256\t256\t1\t-2.0\t0.98\t0.01\t0.95\n"
        )
        (d / "density_cv.tsv").write_text(
            "orientation_group\tcandidate_knn_mass\tknn_cv_loglik\tmedian_cv_radius_bp\t"
            "p95_cv_radius_bp\tsampled_anchor_positions\tselected\tsearch_boundary\n"
            "fw\t128\t-2.1\t100\t200\t1000\t0\t0\n"
            "fw\t256\t-2.0\t200\t400\t1000\t1\t1\n"
        )

        cmd = [
            sys.executable, str(MERGE), "--sample", "S",
            "--contract", str(d / "contract.tsv"),
            "--resolver", str(d / "resolver.tsv"),
            "--realization", str(d / "realization.tsv"),
            "--density", str(d / "density.tsv"),
            "--density-cv", str(d / "density_cv.tsv"),
            "--output-prefix", str(d / "S"),
        ]
        subprocess.run(cmd, check=True)
        merged = read_long(d / "S.mm.qc.tsv")

        assert merged["qc_validation_status"] == "PASS"
        assert merged["density_k_status"] == "BOUNDARY_HIT"
        assert merged["knn_boundary_hit_any"] == "1"
        assert merged["knn_boundary_orientations"] == "fw"
        assert merged["selected_knn_mass_fw"] == "256"
        assert merged["selected_p95_cv_radius_bp_fw"] == "400"
        assert abs(float(merged["resolution_fraction_of_eligible_multimappers"]) - 0.8) < 1e-12
        assert abs(float(merged["highconf_fraction_of_resolved"]) - 0.75) < 1e-12
        assert abs(float(merged["knn_cv_gain_vs_previous_fw"]) - 0.1) < 1e-12

        mqc = (d / "S.mmresolver_mqc.tsv").read_text()
        assert "resolved_highconf_fragments" in mqc
        assert "knn_boundary_hit_any" in mqc
        assert "quantitative_bam_alignment_records" in mqc
        assert "confident_bam_alignment_records" not in mqc

        # Cross-stage contradictions must fail loudly rather than create a
        # superficially valid MultiQC table.
        bad = (d / "resolver_bad.tsv")
        bad.write_text((d / "resolver.tsv").read_text().replace("candidate_rows\t60", "candidate_rows\t61"))
        bad_cmd = cmd.copy()
        bad_cmd[bad_cmd.index(str(d / "resolver.tsv"))] = str(bad)
        proc = subprocess.run(bad_cmd, text=True, capture_output=True)
        assert proc.returncode != 0
        assert "candidate rows" in (proc.stdout + proc.stderr)

    print("PASS PATCH05 cross-stage resolver QC + density-k diagnostics")


if __name__ == "__main__":
    main()

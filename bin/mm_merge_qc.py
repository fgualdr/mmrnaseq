#!/usr/bin/env python3
"""Merge and validate multimapper-resolver QC for one sample.

PATCH 05 is intentionally diagnostic: it does not change posterior inference or
k selection.  It checks accounting across PATCH 01/02P/03 outputs, summarizes
the three BAM tiers, and exposes sample/orientation-specific density-model and
k-CV diagnostics for later MultiQC reporting.
"""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List


class QCError(RuntimeError):
    pass


def read_kv(path: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"metric", "value"}.issubset(reader.fieldnames):
            raise QCError(f"{path}: expected tab-separated metric/value columns")
        for row in reader:
            metric = (row.get("metric") or "").strip()
            if not metric:
                continue
            if metric in out:
                raise QCError(f"{path}: duplicate metric {metric!r}")
            out[metric] = (row.get("value") or "").strip()
    return out


def read_table(path: str) -> List[Dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise QCError(f"{path}: missing header")
        return list(reader)


def require(d: Dict[str, str], keys: Iterable[str], label: str) -> None:
    missing = [k for k in keys if k not in d]
    if missing:
        raise QCError(f"{label}: missing required metrics: {', '.join(missing)}")


def as_int(d: Dict[str, str], key: str, label: str) -> int:
    try:
        x = float(d[key])
    except (KeyError, ValueError) as exc:
        raise QCError(f"{label}: metric {key!r} is not numeric") from exc
    if not math.isfinite(x) or abs(x - round(x)) > 1e-8:
        raise QCError(f"{label}: metric {key!r} must be an integer, got {d[key]!r}")
    return int(round(x))


def as_float(d: Dict[str, str], key: str, label: str) -> float:
    try:
        x = float(d[key])
    except (KeyError, ValueError) as exc:
        raise QCError(f"{label}: metric {key!r} is not numeric") from exc
    if not math.isfinite(x):
        raise QCError(f"{label}: metric {key!r} is not finite, got {d[key]!r}")
    return x


def ratio(num: float, den: float) -> float:
    return num / den if den else 0.0


def fmt(x: float) -> str:
    return format(float(x), ".12g")


def safe_orientation(ori: str) -> str:
    ori = ori.strip()
    if ori in {"", "."}:
        return "unstranded"
    return "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in ori)


def assert_eq(a: int, b: int, what: str) -> None:
    if a != b:
        raise QCError(f"cross-stage accounting failure: {what}: {a} != {b}")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    p.add_argument("--contract", required=True)
    p.add_argument("--resolver", required=True)
    p.add_argument("--realization", required=True)
    p.add_argument("--density", required=True)
    p.add_argument("--density-cv", required=True)
    p.add_argument("--output-prefix", required=True)
    a = p.parse_args()

    c = read_kv(a.contract)
    r = read_kv(a.resolver)
    b = read_kv(a.realization)
    density = read_table(a.density)
    density_cv = read_table(a.density_cv)

    require(c, [
        "total_fragments", "unique_fragments", "multimapper_fragments",
        "excluded_fragments", "candidate_rows",
    ], "candidate contract QC")
    require(r, [
        "multimapper_fragments", "candidate_rows", "resolved_fragments",
        "unresolved_fragments", "absolute_score_failures", "target_posterior_fdr",
        "achieved_expected_fdr", "equal_best_as_fraction", "median_candidate_count",
        "p95_candidate_count", "median_pmax", "median_entropy", "median_as_delta",
    ], "resolver QC")
    require(b, [
        "unique_fragments", "resolved_fdr_fragments", "resolved_highconf_fragments",
        "resolved_lower_confidence_fragments", "unresolved_fragments", "excluded_fragments",
        "highconf_min_pmax", "full_bam_alignment_records",
        "quantitative_bam_alignment_records", "highconf_bam_alignment_records",
    ], "BAM realization QC")

    total = as_int(c, "total_fragments", "candidate contract QC")
    unique = as_int(c, "unique_fragments", "candidate contract QC")
    multi = as_int(c, "multimapper_fragments", "candidate contract QC")
    excluded = as_int(c, "excluded_fragments", "candidate contract QC")
    candidate_rows = as_int(c, "candidate_rows", "candidate contract QC")

    r_multi = as_int(r, "multimapper_fragments", "resolver QC")
    r_rows = as_int(r, "candidate_rows", "resolver QC")
    resolved = as_int(r, "resolved_fragments", "resolver QC")
    unresolved = as_int(r, "unresolved_fragments", "resolver QC")
    abs_fail = as_int(r, "absolute_score_failures", "resolver QC")
    target_fdr = as_float(r, "target_posterior_fdr", "resolver QC")
    achieved_fdr = as_float(r, "achieved_expected_fdr", "resolver QC")

    b_unique = as_int(b, "unique_fragments", "BAM realization QC")
    b_resolved = as_int(b, "resolved_fdr_fragments", "BAM realization QC")
    b_highconf = as_int(b, "resolved_highconf_fragments", "BAM realization QC")
    b_lower = as_int(b, "resolved_lower_confidence_fragments", "BAM realization QC")
    b_unresolved = as_int(b, "unresolved_fragments", "BAM realization QC")
    b_excluded = as_int(b, "excluded_fragments", "BAM realization QC")
    highconf_min = as_float(b, "highconf_min_pmax", "BAM realization QC")
    full_records = as_int(b, "full_bam_alignment_records", "BAM realization QC")
    quant_records = as_int(b, "quantitative_bam_alignment_records", "BAM realization QC")
    high_records = as_int(b, "highconf_bam_alignment_records", "BAM realization QC")

    # Exact cross-stage invariants.  These are scientific accounting checks, not
    # presentation-only QC: a contradiction means downstream quantification is unsafe.
    assert_eq(total, unique + multi + excluded, "total = unique + multimapper + excluded")
    assert_eq(multi, r_multi, "contract multimappers = resolver multimappers")
    assert_eq(candidate_rows, r_rows, "contract candidate rows = resolver candidate rows")
    assert_eq(multi, resolved + unresolved, "multimappers = resolved + unresolved")
    assert_eq(unique, b_unique, "contract unique = realization unique")
    assert_eq(resolved, b_resolved, "resolver resolved = realization global-FDR resolved")
    assert_eq(unresolved, b_unresolved, "resolver unresolved = realization unresolved")
    assert_eq(excluded, b_excluded, "contract excluded = realization excluded")
    assert_eq(b_resolved, b_highconf + b_lower, "resolved = high-confidence + lower-confidence")

    if abs_fail > unresolved:
        raise QCError("absolute-score failures cannot exceed unresolved fragments")
    if achieved_fdr > target_fdr + 1e-9:
        raise QCError(
            f"achieved expected FDR {achieved_fdr} exceeds target {target_fdr}"
        )
    if abs(highconf_min - (1.0 - target_fdr)) > 1e-9:
        raise QCError(
            "high-confidence pmax threshold is inconsistent with resolver FDR target: "
            f"{highconf_min} vs {1.0 - target_fdr}"
        )
    if not (full_records >= quant_records >= high_records):
        raise QCError(
            "BAM record nesting violated: expected full >= quantitative >= high-confidence, got "
            f"{full_records} >= {quant_records} >= {high_records}"
        )

    # Merge named inputs without allowing same-name metrics from different stages
    # to overwrite one another silently.
    merged: Dict[str, str] = {}
    for prefix, source in (("contract", c), ("resolver", r), ("realization", b)):
        for key, value in source.items():
            merged[f"{prefix}_{key}"] = value

    # Stable, stage-independent names used by MultiQC/downstream summaries.
    merged.update({
        "total_fragments": str(total),
        "unique_fragments": str(unique),
        "multimapper_fragments": str(multi),
        "excluded_fragments": str(excluded),
        "candidate_rows": str(candidate_rows),
        "resolved_fdr_fragments": str(resolved),
        "resolved_highconf_fragments": str(b_highconf),
        "resolved_lower_confidence_fragments": str(b_lower),
        "unresolved_fragments": str(unresolved),
        "absolute_score_failures": str(abs_fail),
        "target_posterior_fdr": fmt(target_fdr),
        "achieved_expected_fdr": fmt(achieved_fdr),
        "highconf_min_pmax": fmt(highconf_min),
        "full_bam_alignment_records": str(full_records),
        "quantitative_bam_alignment_records": str(quant_records),
        "highconf_bam_alignment_records": str(high_records),
        "resolution_fraction_of_eligible_multimappers": fmt(ratio(resolved, multi)),
        "highconf_fraction_of_resolved": fmt(ratio(b_highconf, resolved)),
        "lower_confidence_fraction_of_resolved": fmt(ratio(b_lower, resolved)),
        "unresolved_fraction_of_eligible_multimappers": fmt(ratio(unresolved, multi)),
        "absolute_score_failure_fraction_of_multimappers": fmt(ratio(abs_fail, multi)),
        "unique_fraction_of_input_fragments": fmt(ratio(unique, total)),
        "multimapper_fraction_of_input_fragments": fmt(ratio(multi, total)),
        "excluded_fraction_of_input_fragments": fmt(ratio(excluded, total)),
        "mean_candidate_rows_per_multimapper": fmt(ratio(candidate_rows, multi)),
        "expected_fdr_headroom": fmt(target_fdr - achieved_fdr),
        "qc_validation_status": "PASS",
    })

    # Keep established resolver summaries under stable names.
    for key in (
        "equal_best_as_fraction", "median_candidate_count", "p95_candidate_count",
        "median_pmax", "median_entropy", "median_as_delta",
    ):
        merged[key] = r[key]

    cv_by_ori: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in density_cv:
        ori = (row.get("orientation_group") or ".").strip() or "."
        cv_by_ori[ori].append(row)

    boundary_orientations: List[str] = []
    selected_ks: List[float] = []
    selected_p95_radii: List[float] = []
    density_mqc_rows: List[List[str]] = []

    required_density_cols = {
        "orientation_group", "unique_anchor_weight", "background_density_per_bp",
        "selected_knn_mass", "knn_search_max", "knn_boundary_hit", "knn_cv_loglik",
        "absolute_normalized_as_cutoff",
    }
    if density and not required_density_cols.issubset(density[0]):
        missing = sorted(required_density_cols - set(density[0]))
        raise QCError(f"density model missing PATCH05 columns: {', '.join(missing)}")

    for row in density:
        ori = (row.get("orientation_group") or ".").strip() or "."
        safe = safe_orientation(ori)
        rows = cv_by_ori.get(ori, [])
        if not rows:
            raise QCError(f"density CV table has no rows for orientation {ori!r}")

        selected_rows = [x for x in rows if int(float(x.get("selected", "0"))) == 1]
        if len(selected_rows) != 1:
            raise QCError(f"density CV orientation {ori!r} must have exactly one selected row")
        sel = selected_rows[0]

        selected_k = float(row["selected_knn_mass"])
        cv_selected_k = float(sel["candidate_knn_mass"])
        if abs(selected_k - cv_selected_k) > 1e-9:
            raise QCError(
                f"density model/CV selected-k mismatch for {ori}: {selected_k} != {cv_selected_k}"
            )

        search_max = float(row["knn_search_max"])
        boundary = int(float(row["knn_boundary_hit"]))
        expected_boundary = 1 if abs(selected_k - search_max) <= 1e-9 else 0
        if boundary != expected_boundary:
            raise QCError(f"density boundary flag inconsistent for orientation {ori!r}")

        scored = sorted(
            [(float(x["candidate_knn_mass"]), float(x["knn_cv_loglik"])) for x in rows],
            key=lambda z: z[0],
        )
        selected_index = next(i for i, (k, _s) in enumerate(scored) if abs(k - selected_k) <= 1e-9)
        gain_previous = 0.0
        if selected_index > 0:
            gain_previous = scored[selected_index][1] - scored[selected_index - 1][1]
        scores_desc = sorted((score for _k, score in scored), reverse=True)
        margin_second = scores_desc[0] - scores_desc[1] if len(scores_desc) > 1 else 0.0

        selected_ks.append(selected_k)
        selected_p95_radii.append(float(sel["p95_cv_radius_bp"]))
        if boundary:
            boundary_orientations.append(ori)

        merged[f"unique_anchor_weight_{safe}"] = row["unique_anchor_weight"]
        merged[f"background_density_per_bp_{safe}"] = row["background_density_per_bp"]
        merged[f"selected_knn_mass_{safe}"] = fmt(selected_k)
        merged[f"knn_search_max_{safe}"] = fmt(search_max)
        merged[f"knn_boundary_hit_{safe}"] = str(boundary)
        merged[f"knn_cv_loglik_{safe}"] = row["knn_cv_loglik"]
        merged[f"absolute_normalized_as_cutoff_{safe}"] = row["absolute_normalized_as_cutoff"]
        merged[f"selected_median_cv_radius_bp_{safe}"] = sel["median_cv_radius_bp"]
        merged[f"selected_p95_cv_radius_bp_{safe}"] = sel["p95_cv_radius_bp"]
        merged[f"knn_cv_gain_vs_previous_{safe}"] = fmt(gain_previous)
        merged[f"knn_cv_margin_to_second_best_{safe}"] = fmt(margin_second)
        merged[f"knn_cv_sampled_anchor_positions_{safe}"] = sel["sampled_anchor_positions"]

        density_mqc_rows.append([
            f"{a.sample}|{safe}",
            a.sample,
            ori,
            row["unique_anchor_weight"],
            row["background_density_per_bp"],
            fmt(selected_k),
            fmt(search_max),
            str(boundary),
            row["knn_cv_loglik"],
            fmt(gain_previous),
            sel["median_cv_radius_bp"],
            sel["p95_cv_radius_bp"],
            row["absolute_normalized_as_cutoff"],
        ])

    merged["density_orientation_count"] = str(len(density))
    merged["knn_boundary_hit_any"] = "1" if boundary_orientations else "0"
    merged["knn_boundary_orientation_count"] = str(len(boundary_orientations))
    merged["knn_boundary_orientations"] = ",".join(boundary_orientations) if boundary_orientations else "."
    merged["max_selected_knn_mass"] = fmt(max(selected_ks) if selected_ks else 0.0)
    merged["max_selected_p95_cv_radius_bp"] = fmt(max(selected_p95_radii) if selected_p95_radii else 0.0)
    merged["density_k_status"] = "BOUNDARY_HIT" if boundary_orientations else "INTERIOR_OPTIMUM"

    out = Path(a.output_prefix + ".mm.qc.tsv")
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "metric", "value"])
        for key in sorted(merged):
            writer.writerow([a.sample, key, merged[key]])

    # Compact numeric table intended for MultiQC PATCH 07.  All biological
    # counts are fragments; BAM alignment-record counts are explicitly named.
    cols = [
        "total_fragments", "unique_fragments", "multimapper_fragments", "excluded_fragments",
        "resolved_fdr_fragments", "resolved_highconf_fragments",
        "resolved_lower_confidence_fragments", "unresolved_fragments",
        "resolution_fraction_of_eligible_multimappers", "highconf_fraction_of_resolved",
        "absolute_score_failure_fraction_of_multimappers", "equal_best_as_fraction",
        "median_candidate_count", "p95_candidate_count", "median_as_delta", "median_pmax",
        "median_entropy", "target_posterior_fdr", "achieved_expected_fdr",
        "knn_boundary_hit_any", "knn_boundary_orientation_count", "max_selected_knn_mass",
        "max_selected_p95_cv_radius_bp", "full_bam_alignment_records",
        "quantitative_bam_alignment_records", "highconf_bam_alignment_records",
    ]
    mq = Path(a.output_prefix + ".mmresolver_mqc.tsv")
    with mq.open("w", newline="", encoding="utf-8") as handle:
        handle.write(
            "# id: mmresolver\n"
            "# section_name: Multimapper resolver\n"
            "# description: Fragment-level Bayesian STAR-candidate resolution and BAM-tier realization. Fragment counts and BAM alignment-record counts are kept distinct.\n"
            "# plot_type: table\n"
        )
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Sample"] + cols)
        writer.writerow([a.sample] + [merged.get(x, "0") for x in cols])

    dmq = Path(a.output_prefix + ".mmdensity_mqc.tsv")
    with dmq.open("w", newline="", encoding="utf-8") as handle:
        handle.write(
            "# id: mmdensity\n"
            "# section_name: Multimapper density bandwidth\n"
            "# description: Sample- and orientation-specific adaptive unique-anchor density model. A boundary hit means the selected k equals the largest k tested and should be reviewed before changing the resolver grid.\n"
            "# plot_type: table\n"
        )
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "Sample-orientation", "Sample", "Orientation", "Unique anchor weight",
            "Background density / bp", "Selected k", "Search max k", "Boundary hit",
            "CV log-likelihood", "Gain vs previous k", "Median CV radius (bp)",
            "P95 CV radius (bp)", "Absolute normalized AS cutoff",
        ])
        writer.writerows(density_mqc_rows)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except QCError as exc:
        raise SystemExit(f"MM_QC_MERGE ERROR: {exc}")

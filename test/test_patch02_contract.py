#!/usr/bin/env python3
"""Lightweight regression tests for the STAR candidate-table contract.

Run with: python tests/test_patch02_contract.py /path/to/sample.candidate_table.tsv
"""
from __future__ import annotations
import csv, sys
from collections import Counter, defaultdict

EXPECTED = ["alignment_id","fragment_id","target_id","chrom","start","end","hi","nh","as","nM"]

def main(path: str) -> int:
    with open(path, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        assert r.fieldnames == EXPECTED, f"Unexpected header: {r.fieldnames}"
        rows = list(r)
    assert rows, "Candidate table is empty"
    by_fragment = defaultdict(list)
    for row in rows:
        assert all(row[k] != "" for k in EXPECTED), f"Missing field: {row}"
        s, e, hi, nh = map(int, (row["start"], row["end"], row["hi"], row["nh"]))
        int(row["as"]); int(row["nM"])
        assert s >= 0 and e > s
        assert hi >= 1 and nh >= 1
        by_fragment[row["fragment_id"]].append(row)
    checked_multi = 0
    for frag, rr in by_fragment.items():
        nhs = {int(x["nh"]) for x in rr}
        assert len(nhs) == 1, f"NH differs within fragment {frag}: {nhs}"
        nh = next(iter(nhs))
        his = {int(x["hi"]) for x in rr}
        if nh > 1:
            checked_multi += 1
            assert min(his) >= 1 and max(his) <= nh, (frag, sorted(his), nh)
    print(f"PASS rows={len(rows)} fragments={len(by_fragment)} multimapper_fragments_checked={checked_multi}")
    return 0

if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_patch02_contract.py SAMPLE.candidate_table.tsv")
    raise SystemExit(main(sys.argv[1]))

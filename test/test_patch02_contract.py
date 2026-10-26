#!/usr/bin/env python3
"""Streaming regression test for the PATCH 02 STAR candidate-table contract.

Memory use is O(max candidates for one fragment), not O(total table size).
The PATCH 02 table is produced by a sort/join on ``alignment_id``; because
``alignment_id`` is ``fragment_id + _HI:i:<N>``, rows for one fragment are
contiguous.  We verify that ordering while streaming.

Usage:
    python3 test/test_patch02_contract.py SAMPLE.candidate_table.tsv
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

EXPECTED = [
    "alignment_id",
    "fragment_id",
    "target_id",
    "chrom",
    "start",
    "end",
    "hi",
    "nh",
    "as",
    "nM",
]


def fail(message: str, line_no: int | None = None) -> "NoReturn":
    prefix = f"line {line_no}: " if line_no is not None else ""
    raise AssertionError(prefix + message)


def finish_fragment(fragment_id: str | None, nh: int | None, his: set[int]) -> int:
    """Validate one completed fragment group; return 1 if it is multimapping."""
    if fragment_id is None:
        return 0
    assert nh is not None
    if not his:
        fail(f"fragment {fragment_id!r} has no candidate rows")
    if min(his) < 1 or max(his) > nh:
        fail(
            f"fragment {fragment_id!r}: HI outside 1..NH; "
            f"HI={sorted(his)[:20]} NH={nh}"
        )
    return int(nh > 1)


def main(path: str) -> int:
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(path)

    n_rows = 0
    n_fragments = 0
    n_multimapper_fragments = 0

    current_fragment: str | None = None
    current_nh: int | None = None
    current_his: set[int] = set()
    previous_fragment: str | None = None

    with p.open("r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            fail("candidate table is empty")

        if header != EXPECTED:
            fail(f"unexpected header: {header!r}; expected {EXPECTED!r}", 1)

        for line_no, row in enumerate(reader, start=2):
            if len(row) != len(EXPECTED):
                fail(f"expected 10 tab-separated fields, observed {len(row)}: {row!r}", line_no)
            if any(value == "" for value in row):
                fail(f"missing field: {row!r}", line_no)

            (
                alignment_id,
                fragment_id,
                target_id,
                chrom,
                start_s,
                end_s,
                hi_s,
                nh_s,
                as_s,
                nm_s,
            ) = row

            try:
                start = int(start_s)
                end = int(end_s)
                hi = int(hi_s)
                nh = int(nh_s)
                int(as_s)
                int(nm_s)
            except ValueError as exc:
                fail(f"non-integer numeric field ({exc}): {row!r}", line_no)

            if start < 0 or end <= start:
                fail(f"invalid interval {chrom}:{start}-{end}", line_no)
            if hi < 1 or nh < 1:
                fail(f"invalid HI/NH: HI={hi}, NH={nh}", line_no)
            if hi > nh:
                fail(f"HI={hi} exceeds NH={nh}", line_no)
            if not alignment_id.startswith(fragment_id + "_HI:i:"):
                fail(
                    f"alignment_id does not match fragment_id/HI convention: "
                    f"alignment_id={alignment_id!r}, fragment_id={fragment_id!r}",
                    line_no,
                )

            # PATCH 02 joins tables sorted by alignment_id. Since alignment_id starts
            # with fragment_id, each fragment must appear in one contiguous block.
            if current_fragment is None:
                current_fragment = fragment_id
                current_nh = nh
                current_his = set()
                n_fragments = 1
            elif fragment_id != current_fragment:
                n_multimapper_fragments += finish_fragment(
                    current_fragment, current_nh, current_his
                )
                if previous_fragment is not None and fragment_id < previous_fragment:
                    fail(
                        "fragment ordering is not monotonic; the streaming contract "
                        "requires PATCH 02 sort/join output",
                        line_no,
                    )
                previous_fragment = current_fragment
                current_fragment = fragment_id
                current_nh = nh
                current_his = set()
                n_fragments += 1

            if nh != current_nh:
                fail(
                    f"NH differs within fragment {fragment_id!r}: "
                    f"expected {current_nh}, observed {nh}",
                    line_no,
                )
            current_his.add(hi)
            n_rows += 1

    if n_rows == 0:
        fail("candidate table has a header but no data rows")

    n_multimapper_fragments += finish_fragment(
        current_fragment, current_nh, current_his
    )

    print(
        "PASS "
        f"rows={n_rows} "
        f"fragments={n_fragments} "
        f"multimapper_fragments_checked={n_multimapper_fragments}"
    )
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(
            "usage: python3 test/test_patch02_contract.py SAMPLE.candidate_table.tsv"
        )
    raise SystemExit(main(sys.argv[1]))

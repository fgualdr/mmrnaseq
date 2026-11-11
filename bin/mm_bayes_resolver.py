#!/usr/bin/env python3
"""STAR-candidate Bayesian multimapper resolver.

This is the readable Python/NumPy implementation of the PATCH 02 resolver.
It intentionally preserves the statistical contract of the former C++ program:

  posterior(candidate | fragment) ∝ exp(AS - best_AS) × local_unique_density

Local density is estimated separately by observed orientation from uniquely mapped
fragment anchors using a data-selected weighted k-nearest-neighbour window.
Assignments are hard-resolved only under an absolute alignment-quality gate and
a posterior expected-FDR criterion.

The implementation is designed for real RNA-seq scale:
  * input TSVs are streamed (plain text or .gz)
  * candidates are processed in fragment-complete chunks
  * expensive density calculations are vectorised with NumPy
  * only compact numeric arrays are retained for global FDR/QC
  * pandas is deliberately not used
"""
from __future__ import annotations

import argparse
import gzip
import math
import os
import sys
import tempfile
from array import array
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, MutableMapping, Optional, Sequence, Tuple

import numpy as np

EPS = 1e-15


def open_text(path: str):
    """Open plain or gzip-compressed text transparently."""
    return gzip.open(path, "rt", encoding="utf-8", newline="") if str(path).endswith(".gz") else open(path, "rt", encoding="utf-8", newline="")


def parse_header(line: str) -> Dict[str, int]:
    fields = line.rstrip("\n").split("\t")
    return {name: i for i, name in enumerate(fields)}


def require_columns(index: Dict[str, int], required: Sequence[str], label: str) -> None:
    missing = [x for x in required if x not in index]
    if missing:
        raise ValueError(f"{label} missing required columns: {', '.join(missing)}")


@dataclass
class AnchorIndex:
    positions: np.ndarray       # int64, unique sorted genomic positions
    prefix_mass: np.ndarray     # float64, length n+1; prefix_mass[0] == 0
    total_mass: float
    chrom_size: int


@dataclass
class OrientationModel:
    unique_weight: float = 0.0
    lambda0: float = 0.0
    k: float = 1.0
    median_norm: float = 0.0
    mad: float = 0.0
    absolute_cutoff: float = -1e9
    cv_score: float = float("-inf")


def read_chrom_sizes(path: str) -> Tuple[Dict[str, int], float]:
    sizes: Dict[str, int] = {}
    with open_text(path) as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Malformed chrom-sizes line {lineno}: {line.rstrip()}")
            sizes[fields[0]] = int(fields[1])
    genome_size = float(sum(sizes.values()))
    if genome_size <= 0:
        raise ValueError("Genome size is zero")
    return sizes, genome_size


def robust_models_from_unique_fragments(path: str, genome_size: float) -> Dict[str, OrientationModel]:
    """Estimate orientation-specific normalized-AS distributions compactly."""
    scores: MutableMapping[str, array] = defaultdict(lambda: array("f"))
    counts: MutableMapping[str, int] = defaultdict(int)

    with open_text(path) as handle:
        header = handle.readline()
        if not header:
            raise ValueError("Unique-fragment table is empty (missing header)")
        ix = parse_header(header)
        require_columns(ix, ["orientation_group", "normalized_as"], "unique-fragment table")
        for lineno, line in enumerate(handle, 2):
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            try:
                ori = f[ix["orientation_group"]]
                scores[ori].append(float(f[ix["normalized_as"]]))
                counts[ori] += 1
            except (IndexError, ValueError) as exc:
                raise ValueError(f"Malformed unique-fragment line {lineno}") from exc

    models: Dict[str, OrientationModel] = {}
    for ori in sorted(scores):
        vals = np.frombuffer(scores[ori], dtype=np.float32)
        if vals.size == 0:
            continue
        med = float(np.median(vals))
        robust_sd = float(np.median(np.abs(vals.astype(np.float64) - med)) * 1.4826)
        cutoff = med - 3.0 * robust_sd
        if robust_sd == 0.0:
            cutoff = float(np.min(vals))
        weight = float(counts[ori])
        models[ori] = OrientationModel(
            unique_weight=weight,
            lambda0=max(weight / genome_size, 1.0 / genome_size),
            median_norm=med,
            mad=robust_sd,
            absolute_cutoff=cutoff,
        )

    if not models:
        models["."] = OrientationModel(lambda0=1.0 / genome_size, absolute_cutoff=-1e9)
    return models


def load_anchor_indexes(path: str, chrom_sizes: Dict[str, int]) -> Dict[Tuple[str, str], AnchorIndex]:
    """Load weighted unique-fragment anchors in compact numeric arrays.

    Python tuples/objects are avoided for individual anchors. array.array keeps the
    raw representation close to C memory cost until each chromosome/orientation is
    converted to NumPy and indexed.
    """
    raw_pos: MutableMapping[Tuple[str, str], array] = defaultdict(lambda: array("q"))
    raw_w: MutableMapping[Tuple[str, str], array] = defaultdict(lambda: array("d"))

    with open_text(path) as handle:
        header = handle.readline()
        if not header:
            raise ValueError("Unique-anchor table is empty (missing header)")
        ix = parse_header(header)
        require_columns(ix, ["chrom", "pos0", "orientation_group", "weight"], "unique-anchor table")
        for lineno, line in enumerate(handle, 2):
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            try:
                key = (f[ix["chrom"]], f[ix["orientation_group"]])
                raw_pos[key].append(int(f[ix["pos0"]]))
                raw_w[key].append(float(f[ix["weight"]]))
            except (IndexError, ValueError) as exc:
                raise ValueError(f"Malformed unique-anchor line {lineno}") from exc

    anchors: Dict[Tuple[str, str], AnchorIndex] = {}
    for key in sorted(list(raw_pos.keys())):
        pos_buf = raw_pos.pop(key)
        w_buf = raw_w.pop(key)
        pos = np.asarray(pos_buf, dtype=np.int64)
        weight = np.asarray(w_buf, dtype=np.float64)
        if pos.size == 0:
            continue

        order = np.argsort(pos, kind="mergesort")
        pos = pos[order]
        weight = weight[order]
        uniq_pos, starts = np.unique(pos, return_index=True)
        uniq_weight = np.add.reduceat(weight, starts)
        prefix = np.empty(uniq_weight.size + 1, dtype=np.float64)
        prefix[0] = 0.0
        np.cumsum(uniq_weight, out=prefix[1:])
        chrom = key[0]
        anchors[key] = AnchorIndex(
            positions=uniq_pos.astype(np.int64, copy=False),
            prefix_mass=prefix,
            total_mass=float(prefix[-1]),
            chrom_size=int(chrom_sizes.get(chrom, 0)),
        )
    return anchors


def mass_window(anchor: AnchorIndex, x: np.ndarray, radius: np.ndarray) -> np.ndarray:
    lo = np.searchsorted(anchor.positions, x - radius, side="left")
    hi = np.searchsorted(anchor.positions, x + radius, side="right")
    return anchor.prefix_mass[hi] - anchor.prefix_mass[lo]


def radius_for_mass(anchor: AnchorIndex, x: np.ndarray, target_mass) -> np.ndarray:
    """Vectorised minimal symmetric radius containing target weighted anchor mass."""
    x = np.asarray(x, dtype=np.int64)
    if x.size == 0:
        return np.empty(0, dtype=np.int64)
    if anchor.positions.size == 0:
        return np.ones(x.size, dtype=np.int64)

    target = np.broadcast_to(np.asarray(target_mass, dtype=np.float64), x.shape)
    maxr = max(1, int(anchor.chrom_size))
    lo = np.zeros(x.size, dtype=np.int64)
    hi = np.full(x.size, maxr, dtype=np.int64)

    # If the full chromosome contains less than target mass, C++ returns maxr.
    full = mass_window(anchor, x, hi)
    active = full >= target
    if not np.any(active):
        return hi

    # Binary search only active queries. ~31 iterations for human chromosome sizes.
    while True:
        idx = np.flatnonzero(active & (lo < hi))
        if idx.size == 0:
            break
        mid = (lo[idx] + hi[idx]) // 2
        m = mass_window(anchor, x[idx], mid)
        enough = m >= target[idx]
        hi[idx[enough]] = mid[enough]
        lo[idx[~enough]] = mid[~enough] + 1

    out = hi
    out[out < 1] = 1
    return out


def density_at(anchor: Optional[AnchorIndex], x: np.ndarray, k: float, lambda0: float) -> np.ndarray:
    x = np.asarray(x, dtype=np.int64)
    if anchor is None or anchor.positions.size == 0:
        return np.full(x.size, max(lambda0, EPS), dtype=np.float64)
    radius = radius_for_mass(anchor, x, float(k))
    u = mass_window(anchor, x, radius)
    denom = 2.0 * radius.astype(np.float64) + 0.5 / max(lambda0, EPS)
    return np.maximum((u + 0.5) / denom, EPS)


def select_knn_mass(
    models: Dict[str, OrientationModel],
    anchors: Dict[Tuple[str, str], AnchorIndex],
    genome_size: float,
    max_samples: int = 20_000,
) -> List[Tuple[str, float, float, float, float, int, int, int]]:
    """Select orientation-specific weighted k by deterministic leave-one-out likelihood.

    The selected model is unchanged from PATCH 02P.  PATCH 05 additionally
    returns the complete evaluated CV curve so boundary optima and the physical
    genomic scale of each candidate k can be audited per sample/orientation.

    Returned rows are:
      orientation, k, cv_loglik, median_radius_bp, p95_radius_bp,
      sampled_anchor_positions, selected, search_boundary.
    """
    cv_rows: List[Tuple[str, float, float, float, float, int, int, int]] = []

    for ori in sorted(models):
        model = models[ori]
        keys = [key for key in sorted(anchors) if key[1] == ori]
        total = float(sum(anchors[key].total_mass for key in keys))
        point_count = int(sum(anchors[key].positions.size for key in keys))
        model.unique_weight = total
        model.lambda0 = max(total / genome_size, 1.0 / genome_size)

        ks: List[float] = []
        k = 4.0
        while k <= 256.0 and k < total:
            ks.append(k)
            k *= 2.0
        if not ks:
            ks = [max(1.0, total - 1.0)]

        if point_count == 0:
            model.k = ks[0]
            model.cv_score = float("-inf")
            cv_rows.append((ori, ks[0], float("-inf"), 0.0, 0.0, 0, 1, 1))
            continue

        stride = max(1, point_count // max_samples)
        samples: List[Tuple[AnchorIndex, np.ndarray]] = []
        seen = 0
        for key in keys:
            anchor = anchors[key]
            n = anchor.positions.size
            if n == 0:
                continue
            offsets = np.arange(n, dtype=np.int64)
            mask = ((seen + offsets) % stride) == 0
            idx = offsets[mask]
            if idx.size:
                samples.append((anchor, idx))
            seen += n

        best_score = float("-inf")
        best_k = ks[0]
        evaluated: List[Tuple[float, float, float, float, int]] = []
        for candidate_k in ks:
            score = 0.0
            wt = 0.0
            radius_parts: List[np.ndarray] = []
            n_sampled = 0
            for anchor, idx in samples:
                x = anchor.positions[idx]
                self_mass = anchor.prefix_mass[idx + 1] - anchor.prefix_mass[idx]
                radius = radius_for_mass(anchor, x, candidate_k + self_mass)
                u = np.maximum(0.0, mass_window(anchor, x, radius) - self_mass)
                d = np.maximum(
                    (u + 0.5) / (2.0 * radius.astype(np.float64) + 0.5 / model.lambda0),
                    EPS,
                )
                score += float(np.sum(self_mass * np.log(d)))
                wt += float(np.sum(self_mass))
                radius_parts.append(radius)
                n_sampled += int(radius.size)
            if wt > 0:
                score /= wt

            if radius_parts:
                radii = np.concatenate(radius_parts)
                median_radius = float(np.median(radii))
                p95_radius = float(np.quantile(radii, 0.95))
            else:
                median_radius = 0.0
                p95_radius = 0.0

            evaluated.append((candidate_k, score, median_radius, p95_radius, n_sampled))
            if score > best_score:
                best_score = score
                best_k = candidate_k

        model.k = best_k
        model.cv_score = best_score
        search_max = max(ks)
        for candidate_k, score, median_radius, p95_radius, n_sampled in evaluated:
            cv_rows.append((
                ori, candidate_k, score, median_radius, p95_radius, n_sampled,
                1 if candidate_k == best_k else 0,
                1 if candidate_k == search_max else 0,
            ))

    return cv_rows


def parse_blocks_into_queries(
    encoded: str,
    orientation: str,
    candidate_index: int,
    queries: MutableMapping[Tuple[str, str], Tuple[List[int], List[float], List[int]]],
) -> None:
    if not encoded or encoded == ".":
        return
    for item in encoded.split(","):
        try:
            chrom, interval, weight_text = item.rsplit(":", 2)
            start_text, end_text = interval.split("-", 1)
            start = int(start_text)
            end = int(end_text)
            weight = float(weight_text)
        except (ValueError, IndexError) as exc:
            raise ValueError(f"Malformed block encoding: {item}") from exc
        mid = (start + end - 1) // 2
        key = (chrom, orientation)
        if key not in queries:
            queries[key] = ([], [], [])
        queries[key][0].append(mid)
        queries[key][1].append(weight)
        queries[key][2].append(candidate_index)


class CandidateChunk:
    """A fragment-complete candidate chunk plus vectorisable block queries."""
    def __init__(self) -> None:
        self.fragment: List[str] = []
        self.candidate: List[str] = []
        self.hi: List[int] = []
        self.nh: List[int] = []
        self.orientation: List[str] = []
        self.as_score: List[int] = []
        self.normalized_as: List[float] = []
        self.signature: List[str] = []
        self.group_starts: List[int] = []
        self.block_queries: MutableMapping[Tuple[str, str], Tuple[List[int], List[float], List[int]]] = {}
        self.last_fragment: Optional[str] = None

    def __len__(self) -> int:
        return len(self.fragment)

    def add(self, f: List[str], ix: Dict[str, int]) -> None:
        frag = f[ix["fragment_id"]]
        if self.last_fragment != frag:
            self.group_starts.append(len(self.fragment))
            self.last_fragment = frag
        i = len(self.fragment)
        ori = f[ix["orientation_group"]]
        self.fragment.append(frag)
        self.candidate.append(f[ix["candidate_id"]])
        self.hi.append(int(f[ix["hi"]]))
        self.nh.append(int(f[ix["nh"]]))
        self.orientation.append(ori)
        self.as_score.append(int(f[ix["as"]]))
        self.normalized_as.append(float(f[ix["normalized_as"]]))
        self.signature.append(f[ix["candidate_signature"]])
        parse_blocks_into_queries(f[ix["r1_blocks"]], ori, i, self.block_queries)
        parse_blocks_into_queries(f[ix["r2_blocks"]], ori, i, self.block_queries)


def iter_candidate_chunks(path: str, chunk_rows: int) -> Iterator[CandidateChunk]:
    required = [
        "fragment_id", "candidate_id", "hi", "nh", "orientation_group", "as",
        "normalized_as", "r1_blocks", "r2_blocks", "candidate_signature",
    ]
    with open_text(path) as handle:
        header = handle.readline()
        if not header:
            raise ValueError("Candidate table is empty (missing header)")
        ix = parse_header(header)
        require_columns(ix, required, "candidate table")

        chunk = CandidateChunk()
        pending: Optional[List[str]] = None
        for lineno, line in enumerate(handle, 2):
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < len(ix):
                raise ValueError(f"Malformed candidate line {lineno}: expected at least {len(ix)} fields, got {len(f)}")

            frag = f[ix["fragment_id"]]
            # Flush only at a fragment boundary, never split one fragment across chunks.
            if len(chunk) >= chunk_rows and chunk.last_fragment is not None and frag != chunk.last_fragment:
                yield chunk
                chunk = CandidateChunk()
            chunk.add(f, ix)

        if len(chunk):
            yield chunk


def model_for_orientation(models: Dict[str, OrientationModel], ori: str, genome_size: float) -> OrientationModel:
    if ori in models:
        return models[ori]
    return OrientationModel(lambda0=1.0 / genome_size, k=1.0, absolute_cutoff=-1e9)


def process_chunk(
    chunk: CandidateChunk,
    anchors: Dict[Tuple[str, str], AnchorIndex],
    models: Dict[str, OrientationModel],
    genome_size: float,
):
    n = len(chunk)
    starts = np.asarray(chunk.group_starts, dtype=np.int64)
    counts = np.diff(np.append(starts, n)).astype(np.int64)

    log_density_sum = np.zeros(n, dtype=np.float64)
    weight_sum = np.zeros(n, dtype=np.float64)

    for key in sorted(chunk.block_queries):
        mids_l, weights_l, cand_l = chunk.block_queries[key]
        mids = np.asarray(mids_l, dtype=np.int64)
        weights = np.asarray(weights_l, dtype=np.float64)
        cand_idx = np.asarray(cand_l, dtype=np.int64)
        ori = key[1]
        model = model_for_orientation(models, ori, genome_size)
        dens = density_at(anchors.get(key), mids, model.k, model.lambda0)
        np.add.at(log_density_sum, cand_idx, weights * np.log(dens))
        np.add.at(weight_sum, cand_idx, weights)

    log_density = np.empty(n, dtype=np.float64)
    good = weight_sum > 0
    log_density[good] = log_density_sum[good] / weight_sum[good]
    if np.any(~good):
        fallback = np.fromiter(
            (math.log(model_for_orientation(models, chunk.orientation[i], genome_size).lambda0) for i in np.flatnonzero(~good)),
            dtype=np.float64,
            count=int(np.sum(~good)),
        )
        log_density[~good] = fallback

    as_arr = np.asarray(chunk.as_score, dtype=np.int32)
    norm_arr = np.asarray(chunk.normalized_as, dtype=np.float64)
    max_as = np.maximum.reduceat(as_arr, starts)
    max_as_row = np.repeat(max_as, counts)

    tied_best = np.add.reduceat((as_arr == max_as_row).astype(np.int32), starts)
    sentinel = np.iinfo(np.int32).min
    second_as = np.maximum.reduceat(np.where(as_arr < max_as_row, as_arr, sentinel), starts)
    second_as = np.where(second_as == sentinel, max_as, second_as)

    log_weight = (as_arr.astype(np.float64) - max_as_row.astype(np.float64)) + log_density
    max_log_weight = np.maximum.reduceat(log_weight, starts)
    centered = log_weight - np.repeat(max_log_weight, counts)
    exp_weight = np.exp(centered)
    sum_weight = np.add.reduceat(exp_weight, starts)
    posterior = exp_weight / np.repeat(sum_weight, counts)

    pmax = np.maximum.reduceat(posterior, starts)
    max_mask = posterior == np.repeat(pmax, counts)
    row_index = np.arange(n, dtype=np.int64)
    best_idx = np.minimum.reduceat(np.where(max_mask, row_index, n), starts)

    plogp = np.zeros(n, dtype=np.float64)
    positive = posterior > 0
    plogp[positive] = posterior[positive] * np.log(posterior[positive])
    entropy = -np.add.reduceat(plogp, starts)

    best_norm = norm_arr[best_idx]
    best_ori = [chunk.orientation[int(i)] for i in best_idx]
    cutoff = np.asarray([model_for_orientation(models, ori, genome_size).absolute_cutoff for ori in best_ori], dtype=np.float64)
    absolute_pass = best_norm >= cutoff

    return {
        "starts": starts,
        "counts": counts,
        "max_as": max_as,
        "second_as": second_as,
        "tied_best": tied_best,
        "log_density": log_density,
        "log_weight": log_weight,
        "posterior": posterior,
        "pmax": pmax,
        "best_idx": best_idx,
        "entropy": entropy,
        "best_norm": best_norm,
        "best_ori": best_ori,
        "absolute_pass": absolute_pass,
    }


def fmt(x: float) -> str:
    return format(float(x), ".12g")


def write_density_model(
    path: str,
    models: Dict[str, OrientationModel],
    cv_rows: Sequence[Tuple[str, float, float, float, float, int, int, int]],
) -> None:
    search_max = {}
    for ori, k, *_rest in cv_rows:
        search_max[ori] = max(search_max.get(ori, k), k)

    with open(path, "wt", encoding="utf-8") as out:
        out.write(
            "orientation_group\tunique_anchor_weight\tbackground_density_per_bp\tselected_knn_mass\t"
            "knn_search_max\tknn_boundary_hit\tknn_cv_loglik\tunique_normalized_as_median\t"
            "unique_normalized_as_mad\tabsolute_normalized_as_cutoff\n"
        )
        for ori in sorted(models):
            m = models[ori]
            kmax = float(search_max.get(ori, m.k))
            boundary = 1 if m.k == kmax else 0
            out.write(
                f"{ori}\t{fmt(m.unique_weight)}\t{fmt(m.lambda0)}\t{fmt(m.k)}\t{fmt(kmax)}\t{boundary}\t"
                f"{fmt(m.cv_score)}\t{fmt(m.median_norm)}\t{fmt(m.mad)}\t{fmt(m.absolute_cutoff)}\n"
            )


def write_density_cv(
    path: str,
    cv_rows: Sequence[Tuple[str, float, float, float, float, int, int, int]],
) -> None:
    with open(path, "wt", encoding="utf-8") as out:
        out.write(
            "orientation_group\tcandidate_knn_mass\tknn_cv_loglik\tmedian_cv_radius_bp\t"
            "p95_cv_radius_bp\tsampled_anchor_positions\tselected\tsearch_boundary\n"
        )
        for ori, k, score, median_radius, p95_radius, n_sampled, selected, boundary in cv_rows:
            out.write(
                f"{ori}\t{fmt(k)}\t{fmt(score)}\t{fmt(median_radius)}\t{fmt(p95_radius)}\t"
                f"{n_sampled}\t{selected}\t{boundary}\n"
            )


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--unique-fragments", required=True)
    parser.add_argument("--anchors", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--resolution-fdr", type=float, default=0.05)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--chunk-rows", type=int, default=100_000,
                        help="Maximum candidate rows per processing chunk; fragments are never split [default: 100000]")
    args = parser.parse_args(argv)

    if not (0.0 < args.resolution_fdr < 1.0):
        parser.error("--resolution-fdr must be in (0,1)")
    if args.chunk_rows < 1:
        parser.error("--chunk-rows must be >=1")

    chrom_sizes, genome_size = read_chrom_sizes(args.chrom_sizes)
    models = robust_models_from_unique_fragments(args.unique_fragments, genome_size)
    anchors = load_anchor_indexes(args.anchors, chrom_sizes)
    density_cv_rows = select_knn_mass(models, anchors, genome_size)

    prefix = args.output_prefix
    posterior_path = prefix + ".posterior_alignment.tsv"
    provisional_fd, provisional_path = tempfile.mkstemp(prefix="mmresolver_", suffix=".provisional.tsv", dir=".")
    os.close(provisional_fd)

    errors = array("d")
    abs_pass_arr = array("B")
    cand_counts = array("I")
    pmax_arr = array("d")
    entropy_arr = array("d")
    delta_arr = array("i")
    candidate_rows = 0
    equal_best_fragments = 0

    try:
        with open(posterior_path, "wt", encoding="utf-8", buffering=1024 * 1024) as post, \
             open(provisional_path, "wt", encoding="utf-8", buffering=1024 * 1024) as provisional:
            post.write(
                "fragment_id\tcandidate_id\thi\tnh\torientation_group\tas\tnormalized_as\tlocal_density\t"
                "log_sequence_weight\tposterior\tcandidate_signature\n"
            )
            provisional.write(
                "fragment_id\tselected_candidate_id\tselected_hi\tnh\tcandidate_count\torientation_group\t"
                "pmax\tposterior_error\tentropy\tbest_as\tsecond_as\tas_delta\ttied_best\tnormalized_as\tabsolute_score_pass\n"
            )

            for chunk in iter_candidate_chunks(args.candidates, args.chunk_rows):
                stat = process_chunk(chunk, anchors, models, genome_size)
                n = len(chunk)
                candidate_rows += n

                posterior_lines: List[str] = []
                max_as_row = np.repeat(stat["max_as"], stat["counts"])
                for i in range(n):
                    posterior_lines.append(
                        "\t".join([
                            chunk.fragment[i], chunk.candidate[i], str(chunk.hi[i]), str(chunk.nh[i]),
                            chunk.orientation[i], str(chunk.as_score[i]), fmt(chunk.normalized_as[i]),
                            fmt(math.exp(float(stat["log_density"][i]))),
                            str(int(chunk.as_score[i] - int(max_as_row[i]))),
                            fmt(stat["posterior"][i]), chunk.signature[i],
                        ])
                    )
                post.write("\n".join(posterior_lines) + "\n")

                for g, best_i_np in enumerate(stat["best_idx"]):
                    best_i = int(best_i_np)
                    count = int(stat["counts"][g])
                    pmax = float(stat["pmax"][g])
                    error = 1.0 - pmax
                    absolute_pass = bool(stat["absolute_pass"][g])
                    best_as = int(stat["max_as"][g])
                    second_as = int(stat["second_as"][g])
                    tied = int(stat["tied_best"][g])
                    if tied > 1:
                        equal_best_fragments += 1

                    provisional.write(
                        "\t".join([
                            chunk.fragment[best_i], chunk.candidate[best_i], str(chunk.hi[best_i]), str(chunk.nh[best_i]),
                            str(count), chunk.orientation[best_i], fmt(pmax), fmt(error), fmt(stat["entropy"][g]),
                            str(best_as), str(second_as), str(best_as - second_as), str(tied),
                            fmt(stat["best_norm"][g]), "1" if absolute_pass else "0",
                        ]) + "\n"
                    )
                    errors.append(error)
                    abs_pass_arr.append(1 if absolute_pass else 0)
                    cand_counts.append(count)
                    pmax_arr.append(pmax)
                    entropy_arr.append(float(stat["entropy"][g]))
                    delta_arr.append(best_as - second_as)

        # Global posterior expected-FDR selection.
        err = np.frombuffer(errors, dtype=np.float64)
        abs_ok = np.frombuffer(abs_pass_arr, dtype=np.uint8).astype(bool, copy=False)
        resolved_mask = np.zeros(err.size, dtype=bool)
        valid = np.flatnonzero(abs_ok)
        keep = 0
        achieved = 0.0
        if valid.size:
            order = valid[np.argsort(err[valid], kind="mergesort")]
            cumulative = np.cumsum(err[order])
            running_fdr = cumulative / np.arange(1, order.size + 1, dtype=np.float64)
            accepted = np.flatnonzero(running_fdr <= args.resolution_fdr)
            if accepted.size:
                keep = int(accepted[-1]) + 1
                selected = order[:keep]
                resolved_mask[selected] = True
                achieved = float(np.mean(err[selected]))

        resolution_path = prefix + ".resolution.tsv"
        abs_failures = int(np.sum(~abs_ok)) if err.size else 0
        with open(provisional_path, "rt", encoding="utf-8") as inp, \
             open(resolution_path, "wt", encoding="utf-8", buffering=1024 * 1024) as out:
            inp.readline()  # provisional header
            out.write(
                "fragment_id\tstatus\tselected_candidate_id\tselected_hi\tnh\tcandidate_count\torientation_group\t"
                "pmax\tposterior_error\tentropy\tbest_as\tsecond_as\tas_delta\ttied_best\tnormalized_as\tabsolute_score_pass\n"
            )
            for i, line in enumerate(inp):
                fields = line.rstrip("\n").split("\t")
                status = "RESOLVED" if resolved_mask[i] else "UNRESOLVED"
                out.write(fields[0] + "\t" + status + "\t" + "\t".join(fields[1:]) + "\n")

        write_density_model(prefix + ".density_model.tsv", models, density_cv_rows)
        write_density_cv(prefix + ".density_cv.tsv", density_cv_rows)

        nfrag = err.size
        ccount = np.frombuffer(cand_counts, dtype=np.uint32)
        pmax_np = np.frombuffer(pmax_arr, dtype=np.float64)
        entropy_np = np.frombuffer(entropy_arr, dtype=np.float64)
        delta_np = np.frombuffer(delta_arr, dtype=np.int32)

        def q(arr_: np.ndarray, prob: float) -> float:
            return float(np.quantile(arr_, prob)) if arr_.size else 0.0

        with open(prefix + ".resolver_qc.tsv", "wt", encoding="utf-8") as qc:
            qc.write("metric\tvalue\n")
            qc.write(f"multimapper_fragments\t{nfrag}\n")
            qc.write(f"candidate_rows\t{candidate_rows}\n")
            qc.write(f"resolved_fragments\t{int(np.sum(resolved_mask))}\n")
            qc.write(f"unresolved_fragments\t{int(nfrag - np.sum(resolved_mask))}\n")
            qc.write(f"absolute_score_failures\t{abs_failures}\n")
            qc.write(f"target_posterior_fdr\t{fmt(args.resolution_fdr)}\n")
            qc.write(f"achieved_expected_fdr\t{fmt(achieved)}\n")
            qc.write(f"equal_best_as_fraction\t{fmt(equal_best_fragments / nfrag if nfrag else 0.0)}\n")
            qc.write(f"median_candidate_count\t{fmt(q(ccount, 0.5))}\n")
            qc.write(f"p95_candidate_count\t{fmt(q(ccount, 0.95))}\n")
            qc.write(f"median_pmax\t{fmt(q(pmax_np, 0.5))}\n")
            qc.write(f"median_entropy\t{fmt(q(entropy_np, 0.5))}\n")
            qc.write(f"median_as_delta\t{fmt(q(delta_np, 0.5))}\n")

    finally:
        try:
            os.unlink(provisional_path)
        except FileNotFoundError:
            pass

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        raise SystemExit(141)

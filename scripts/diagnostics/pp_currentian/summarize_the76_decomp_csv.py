#!/usr/bin/env python3
"""Summarize matched-row score and feature parity from THE-76 decomposition CSV."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path


RECO_BINS = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]
FEATURE_FIELDS = [
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e11_over_e33",
    "e32_over_e35",
    "cluster_prob",
    "npb_score",
]


def as_float(row: dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, "nan"))
    except Exception:
        return float("nan")


def as_int(row: dict[str, str], key: str) -> int:
    try:
        return int(float(row.get(key, "-999")))
    except Exception:
        return -999


def finite(x: float) -> bool:
    return math.isfinite(x)


def find_bin(x: float) -> int:
    for i in range(len(RECO_BINS) - 1):
        if RECO_BINS[i] < x < RECO_BINS[i + 1]:
            return i
    return -1


def label(ib: int) -> str:
    return f"{RECO_BINS[ib]}-{RECO_BINS[ib + 1]}" if 0 <= ib < len(RECO_BINS) - 1 else "out"


def mean(xs: list[float]) -> float:
    vals = [x for x in xs if finite(x)]
    return statistics.mean(vals) if vals else float("nan")


def median(xs: list[float]) -> float:
    vals = [x for x in xs if finite(x)]
    return statistics.median(vals) if vals else float("nan")


def fmt(x: float) -> str:
    return "nan" if not finite(x) else f"{x:.6g}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv_path")
    args = ap.parse_args()

    path = Path(args.csv_path)
    rows = list(csv.DictReader(path.open()))
    exact = [r for r in rows if as_int(r, "matched_ppg12_pass") == 1]
    near = [r for r in rows if r.get("category") == "same_track_pass_cluster_dr002_to_dr010"]

    print(f"# {path}")
    print(f"rows={len(rows)} exact_matches={len(exact)} near_same_track={len(near)}")

    cats = Counter(r.get("category", "") for r in rows)
    print("categories=" + ",".join(f"{k}:{v}" for k, v in cats.most_common()))

    tag_matrix: Counter[tuple[int, int]] = Counter()
    by_bin: dict[int, Counter] = defaultdict(Counter)
    score_diffs: list[float] = []
    for r in exact:
        rj_tag = as_int(r, "tag")
        pp_tag = as_int(r, "ppg12_tag")
        tag_matrix[(rj_tag, pp_tag)] += 1
        ib = find_bin(as_float(r, "cluster_Et"))
        by_bin[ib]["total"] += 1
        by_bin[ib][f"rj_tag:{rj_tag}"] += 1
        by_bin[ib][f"ppg12_tag:{pp_tag}"] += 1
        if rj_tag == pp_tag:
            by_bin[ib]["tag_agree"] += 1
        else:
            by_bin[ib]["tag_disagree"] += 1
        diff = as_float(r, "bdt_score") - as_float(r, "ppg12_selected_bdt_score")
        if finite(diff):
            score_diffs.append(diff)

    print("tag_matrix_rj_to_ppg12=" + ",".join(f"{k[0]}->{k[1]}:{v}" for k, v in sorted(tag_matrix.items())))
    print(
        "score_diff_rj_minus_ppg12="
        f"n={len(score_diffs)} mean={fmt(mean(score_diffs))} median={fmt(median(score_diffs))} "
        f"mean_abs={fmt(mean([abs(x) for x in score_diffs]))}"
    )

    print("## matched_by_bin")
    print("bin,total,tag_agree,tag_disagree,rj_tight,rj_nontight,ppg12_tight,ppg12_nontight,score_diff_mean")
    for ib in range(len(RECO_BINS) - 1):
        subset = [r for r in exact if find_bin(as_float(r, "cluster_Et")) == ib]
        diffs = [
            as_float(r, "bdt_score") - as_float(r, "ppg12_selected_bdt_score")
            for r in subset
        ]
        c = by_bin.get(ib, Counter())
        print(
            f"{label(ib)},{c.get('total', 0)},{c.get('tag_agree', 0)},{c.get('tag_disagree', 0)},"
            f"{c.get('rj_tag:1', 0)},{c.get('rj_tag:2', 0)},"
            f"{c.get('ppg12_tag:1', 0)},{c.get('ppg12_tag:2', 0)},"
            f"{fmt(mean(diffs))}"
        )

    print("## feature_diff_exact_matches")
    print("feature,n,mean_rj,mean_ppg12,mean_diff,mean_abs_diff")
    for feature in FEATURE_FIELDS:
        pairs = []
        for r in exact:
            rv = as_float(r, feature)
            pv = as_float(r, f"ppg12_{feature}")
            if finite(rv) and finite(pv):
                pairs.append((rv, pv))
        diffs = [a - b for a, b in pairs]
        print(
            f"{feature},{len(pairs)},{fmt(mean([a for a, _ in pairs]))},"
            f"{fmt(mean([b for _, b in pairs]))},{fmt(mean(diffs))},"
            f"{fmt(mean([abs(x) for x in diffs]))}"
        )

    disagreements = [r for r in exact if as_int(r, "tag") != as_int(r, "ppg12_tag")]
    if disagreements:
        print("## tag_disagreement_examples")
        print("row,event,track,cluster_Et,rj_tag,ppg12_tag,rj_score,ppg12_score,dr")
        for r in disagreements[:20]:
            print(
                f"{r.get('row')},{r.get('eventnumber')},{r.get('truth_track_id')},"
                f"{r.get('cluster_Et')},{r.get('tag')},{r.get('ppg12_tag')},"
                f"{r.get('bdt_score')},{r.get('ppg12_selected_bdt_score')},"
                f"{r.get('nearest_pass_dr')}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Compare PPG12 data slimtree clusters with RecoilJets pp-data diagnostic rows.

THE-76 uses this after a tiny pp-data diagnostic canary with
RJ_PP_PHOTONID_TRAINING_TREE=1.  It is intentionally read-only: the script
matches clusters by run/event and nearest eta/phi, then reports whether
PhotonClusterBuilder/RecoilJets shower-shape inputs and raw topo isolation
match PPG12 CaloAna24 for the same pp data clusters. ABCD flips deliberately
hold the RecoilJets PPG12-compatible tight/non-tight tag fixed so this helper
isolates the effect of the two raw-Eiso values.
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import ROOT  # type: ignore


ROOT.gROOT.SetBatch(True)

FEATURES = [
    ("cluster_Et", "cluster_Et"),
    ("cluster_Et_score_input", "cluster_Et"),
    ("cluster_Eta", "cluster_Eta"),
    ("cluster_Phi", "cluster_Phi"),
    ("vertexz", "vertexz"),
    ("cluster_weta_cogx", "cluster_weta_cogx"),
    ("cluster_wphi_cogx", "cluster_wphi_cogx"),
    ("e11_over_e33", "e11_over_e33"),
    ("e32_over_e35", "e32_over_e35"),
    ("cluster_et1", "cluster_et1"),
    ("cluster_et2", "cluster_et2"),
    ("cluster_et3", "cluster_et3"),
    ("cluster_et4", "cluster_et4"),
    ("cluster_prob", "cluster_prob"),
    ("npb_score", "npb_score"),
    ("tight_bdt_score", "bdt_base_v3E"),
    ("ppg12_raw_eiso", "cluster_iso_topo_04"),
]

ABCD_CATEGORIES = ("A", "B", "C", "D")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rj-root", required=True, help="RecoilJets ROOT with AuAuPhotonIDTrainingTree")
    parser.add_argument(
        "--ppg12-root-glob",
        required=True,
        help="PPG12 data slimtree ROOT path or glob, e.g. part_*_with_bdt_split.root",
    )
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-md", required=True)
    parser.add_argument(
        "--out-summary-csv",
        help="Compact one-row period/tag summary (default: <out-csv stem>_summary.csv)",
    )
    parser.add_argument("--period", default="unknown", help="Canary period label, e.g. 0mrad or 1p5mrad")
    parser.add_argument("--tag", default="unknown", help="Canary configuration/CDB tag label")
    parser.add_argument("--max-rj-rows", type=int, default=0, help="0 means scan all RecoilJets rows")
    parser.add_argument("--max-ppg12-events", type=int, default=0, help="0 means scan all PPG12 entries")
    parser.add_argument(
        "--include-noncommon",
        action="store_true",
        help=(
            "Also compare rows failing ppg12_common_pass or outside tight/non-tight tags 1/2; "
            "the default restricts the ABCD audit to canonical common-pass candidates"
        ),
    )
    parser.add_argument("--et-min", type=float, default=22.0)
    parser.add_argument("--et-max", type=float, default=28.0)
    parser.add_argument("--eta-max", type=float, default=0.7)
    parser.add_argument("--max-dr", type=float, default=0.02)
    parser.add_argument(
        "--max-det",
        type=float,
        default=0.5,
        help="Maximum |RecoilJets ET - PPG12 ET| for run/kinematics fallback matching",
    )
    parser.add_argument("--node", default="CLUSTERINFO_CEMC")
    return parser.parse_args()


def finite(x: Any) -> bool:
    try:
        return math.isfinite(float(x))
    except Exception:
        return False


def dphi(a: float, b: float) -> float:
    x = a - b
    while x > math.pi:
        x -= 2.0 * math.pi
    while x <= -math.pi:
        x += 2.0 * math.pi
    return x


def isolation_region(raw_eiso: float, et: float) -> str:
    """Classify the PPG12 data raw-Eiso window using strict inequalities."""
    if not (finite(raw_eiso) and finite(et)):
        return "INVALID"
    iso_upper = 0.490 + 0.037 * et
    noniso_lower = iso_upper + 0.8
    if -20.0 < raw_eiso < iso_upper:
        return "ISO"
    if noniso_lower < raw_eiso < 20.0:
        return "NONISO"
    if iso_upper <= raw_eiso <= noniso_lower:
        return "GAP"
    return "OUTSIDE"


def abcd_category(tight_tag: int, raw_eiso: float, et: float) -> str:
    """Return A/B/C/D while holding the supplied tight/non-tight tag fixed."""
    region = isolation_region(raw_eiso, et)
    if tight_tag == 1:
        if region == "ISO":
            return "A"
        if region == "NONISO":
            return "B"
    elif tight_tag == 2:
        if region == "ISO":
            return "C"
        if region == "NONISO":
            return "D"
    else:
        return "OTHER_ID"
    return region


def branch_names(tree: Any) -> set[str]:
    return {b.GetName() for b in tree.GetListOfBranches()}


def get_scalar(tree: Any, name: str, default: float = float("nan")) -> float:
    if not hasattr(tree, name):
        return default
    try:
        return float(getattr(tree, name))
    except Exception:
        return default


def get_int(tree: Any, name: str, default: int = -1) -> int:
    if not hasattr(tree, name):
        return default
    try:
        return int(getattr(tree, name))
    except Exception:
        return default


def open_rj_tree(path: str) -> Any:
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open RecoilJets ROOT: {path}")
    tree = f.Get("AuAuPhotonIDTrainingTree")
    if not tree:
        raise RuntimeError(f"Missing AuAuPhotonIDTrainingTree in {path}")
    tree._owner_file = f
    return tree


def open_ppg12_chain(pattern: str) -> Any:
    paths = sorted(glob.glob(pattern))
    if not paths and Path(pattern).exists():
        paths = [pattern]
    if not paths:
        raise RuntimeError(f"No PPG12 ROOT files matched: {pattern}")
    chain = ROOT.TChain("slimtree")
    for path in paths:
        chain.Add(path)
    if chain.GetEntries() <= 0:
        raise RuntimeError(f"PPG12 chain has no entries for {pattern}")
    return chain


def read_rj_rows(tree: Any, args: argparse.Namespace) -> tuple[list[dict[str, float]], set[tuple[int, int]]]:
    names = branch_names(tree)
    required = {
        "run",
        "eventnumber",
        "cluster_Et",
        "cluster_Eta",
        "cluster_Phi",
        "e11_over_e33",
        "ppg12_raw_eiso",
        "ppg12_tight_tag",
    }
    missing = sorted(required - names)
    if missing:
        raise RuntimeError(f"RecoilJets tree missing required branches: {missing}")

    rows: list[dict[str, float]] = []
    events: set[tuple[int, int]] = set()
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        et = get_scalar(tree, "cluster_Et")
        eta = get_scalar(tree, "cluster_Eta")
        if not (args.et_min < et < args.et_max and abs(eta) < args.eta_max):
            continue
        row = {
            "source_entry": i,
            "run": get_int(tree, "run"),
            "eventnumber": get_int(tree, "eventnumber"),
            "cluster_index": get_int(tree, "cluster_index"),
            "cluster_Et": et,
            "cluster_Et_score_input": get_scalar(tree, "cluster_Et_score_input", et),
            "cluster_Eta": eta,
            "cluster_Phi": get_scalar(tree, "cluster_Phi"),
            "vertexz": get_scalar(tree, "vertexz"),
            "cluster_weta_cogx": get_scalar(tree, "cluster_weta_cogx"),
            "cluster_wphi_cogx": get_scalar(tree, "cluster_wphi_cogx"),
            "e11_over_e33": get_scalar(tree, "e11_over_e33"),
            "e32_over_e35": get_scalar(tree, "e32_over_e35"),
            "cluster_et1": get_scalar(tree, "cluster_et1"),
            "cluster_et2": get_scalar(tree, "cluster_et2"),
            "cluster_et3": get_scalar(tree, "cluster_et3"),
            "cluster_et4": get_scalar(tree, "cluster_et4"),
            "cluster_prob": get_scalar(tree, "cluster_prob"),
            "npb_score": get_scalar(tree, "npb_score"),
            "tight_bdt_score": get_scalar(tree, "tight_bdt_score"),
            "ppg12_common_pass": get_int(tree, "ppg12_common_pass"),
            "ppg12_tight_tag": get_int(tree, "ppg12_tight_tag"),
            "ppg12_raw_eiso": get_scalar(tree, "ppg12_raw_eiso"),
        }
        if row["run"] <= 0 or row["eventnumber"] <= 0:
            continue
        if not args.include_noncommon and (
            int(row["ppg12_common_pass"]) != 1
            or int(row["ppg12_tight_tag"]) not in (1, 2)
        ):
            continue
        rows.append(row)
        events.add((int(row["run"]), int(row["eventnumber"])))
        if args.max_rj_rows > 0 and len(rows) >= args.max_rj_rows:
            break
    return rows, events


def array_value(tree: Any, branch: str, idx: int, default: float = float("nan")) -> float:
    if not hasattr(tree, branch):
        return default
    arr = getattr(tree, branch)
    try:
        return float(arr[idx])
    except Exception:
        return default


def read_ppg12_candidates(
    chain: Any,
    wanted_events: set[tuple[int, int]],
    args: argparse.Namespace,
) -> tuple[
    dict[tuple[int, int], list[dict[str, float]]],
    dict[int, list[dict[str, float]]],
]:
    names = branch_names(chain)
    node = args.node
    ncluster_name = f"ncluster_{node}"
    required = {
        "runnumber",
        "eventnumber",
        "vertexz",
        ncluster_name,
        f"cluster_Et_{node}",
        f"cluster_Eta_{node}",
        f"cluster_Phi_{node}",
        f"cluster_e11_{node}",
        f"cluster_e33_{node}",
        f"cluster_iso_topo_04_{node}",
    }
    missing = sorted(required - names)
    if missing:
        raise RuntimeError(f"PPG12 slimtree missing required branches: {missing}")

    # The PPG12 reference tree has hundreds of branches.  This audit needs a
    # small, fixed family; pruning here makes repeated exact canary comparisons
    # fast without changing the selected rows or values.
    active = {
        "runnumber",
        "eventnumber",
        "vertexz",
        ncluster_name,
    }
    for stem in (
        "cluster_Et",
        "cluster_Eta",
        "cluster_Phi",
        "cluster_e11",
        "cluster_e33",
        "cluster_e32",
        "cluster_e35",
        "cluster_weta_cogx",
        "cluster_wphi_cogx",
        "cluster_et1",
        "cluster_et2",
        "cluster_et3",
        "cluster_et4",
        "cluster_prob",
        "cluster_npb_score",
        "cluster_bdt",
        "cluster_iso_topo_04",
    ):
        active.add(f"{stem}_{node}")
    # The BDT branch has the model suffix after the node name.
    active.discard(f"cluster_bdt_{node}")
    active.add(f"cluster_bdt_{node}_base_v3E")
    chain.SetBranchStatus("*", 0)
    for name in sorted(active & names):
        chain.SetBranchStatus(name, 1)

    by_event: dict[tuple[int, int], list[dict[str, float]]] = defaultdict(list)
    by_run: dict[int, list[dict[str, float]]] = defaultdict(list)
    wanted_runs = {run for run, _ in wanted_events}
    stop = chain.GetEntries()
    if args.max_ppg12_events > 0:
        stop = min(stop, args.max_ppg12_events)

    for entry in range(stop):
        chain.GetEntry(entry)
        run = get_int(chain, "runnumber")
        event = get_int(chain, "eventnumber")
        event_key = (run, event)
        if run not in wanted_runs:
            continue
        vertexz = get_scalar(chain, "vertexz")
        ncluster = min(get_int(chain, ncluster_name, 0), 20000)
        for ic in range(ncluster):
            et = array_value(chain, f"cluster_Et_{node}", ic)
            eta = array_value(chain, f"cluster_Eta_{node}", ic)
            if not (args.et_min < et < args.et_max and abs(eta) < args.eta_max):
                continue
            phi = array_value(chain, f"cluster_Phi_{node}", ic)
            e11 = array_value(chain, f"cluster_e11_{node}", ic)
            e33 = array_value(chain, f"cluster_e33_{node}", ic)
            e32 = array_value(chain, f"cluster_e32_{node}", ic)
            e35 = array_value(chain, f"cluster_e35_{node}", ic)
            row = {
                "source_entry": entry,
                "cluster_index": ic,
                "run": run,
                "eventnumber": event,
                "cluster_Et": et,
                "cluster_Eta": eta,
                "cluster_Phi": phi,
                "vertexz": vertexz,
                "cluster_weta_cogx": array_value(chain, f"cluster_weta_cogx_{node}", ic),
                "cluster_wphi_cogx": array_value(chain, f"cluster_wphi_cogx_{node}", ic),
                "e11_over_e33": e11 / e33 if e33 > 0 else float("nan"),
                "e32_over_e35": e32 / e35 if e35 > 0 else float("nan"),
                "cluster_et1": array_value(chain, f"cluster_et1_{node}", ic),
                "cluster_et2": array_value(chain, f"cluster_et2_{node}", ic),
                "cluster_et3": array_value(chain, f"cluster_et3_{node}", ic),
                "cluster_et4": array_value(chain, f"cluster_et4_{node}", ic),
                "cluster_prob": array_value(chain, f"cluster_prob_{node}", ic),
                "npb_score": array_value(chain, f"cluster_npb_score_{node}", ic),
                "bdt_base_v3E": array_value(chain, f"cluster_bdt_{node}_base_v3E", ic),
                "cluster_iso_topo_04": array_value(chain, f"cluster_iso_topo_04_{node}", ic),
            }
            by_event[event_key].append(row)
            by_run[run].append(row)
    return by_event, by_run


def match_rows(
    rj_rows: list[dict[str, float]],
    ppg12_by_event: dict[tuple[int, int], list[dict[str, float]]],
    ppg12_by_run: dict[int, list[dict[str, float]]],
    max_dr: float,
    max_det: float,
) -> tuple[
    list[tuple[dict[str, float], dict[str, float], float, str]],
    list[dict[str, float]],
]:
    matches: list[tuple[dict[str, float], dict[str, float], float, str]] = []
    misses: list[dict[str, float]] = []
    used: set[tuple[int, int]] = set()
    for rj in rj_rows:
        event_key = (int(rj["run"]), int(rj["eventnumber"]))
        candidates = ppg12_by_event.get(event_key, [])
        match_basis = "event"
        if not candidates:
            # CaloAna24's `eventnumber` and RecoilJets' EventHeader sequence are
            # not the same identifier in the paired pp DSTs.  Fall back to a
            # strict same-run ET/eta/phi match and record that basis explicitly.
            candidates = ppg12_by_run.get(event_key[0], [])
            match_basis = "run_kinematics"
        best_dist = float("inf")
        best_score = float("inf")
        best_row: dict[str, float] | None = None
        for ppg in candidates:
            key = (int(ppg["source_entry"]), int(ppg["cluster_index"]))
            if key in used:
                continue
            deta = float(rj["cluster_Eta"]) - float(ppg["cluster_Eta"])
            dph = dphi(float(rj["cluster_Phi"]), float(ppg["cluster_Phi"]))
            dist = math.hypot(deta, dph)
            det = abs(float(rj["cluster_Et"]) - float(ppg["cluster_Et"]))
            if dist > max_dr or det > max_det:
                continue
            score = (dist / max_dr) ** 2 + (det / max_det) ** 2
            if score < best_score:
                best_score = score
                best_dist = dist
                best_row = ppg
        if best_row is None:
            misses.append(rj)
            continue
        used.add((int(best_row["source_entry"]), int(best_row["cluster_index"])))
        matches.append((rj, best_row, best_dist, match_basis))
    return matches, misses


def mean(values: list[float]) -> float:
    clean = [x for x in values if finite(x)]
    return sum(clean) / len(clean) if clean else float("nan")


def write_outputs(
    matches: list[tuple[dict[str, float], dict[str, float], float, str]],
    misses: list[dict[str, float]],
    out_csv: Path,
    out_md: Path,
    out_summary_csv: Path,
    period: str,
    tag: str,
) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_summary_csv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "period",
        "tag",
        "run",
        "eventnumber",
        "rj_entry",
        "ppg12_entry",
        "rj_cluster_index",
        "ppg12_cluster_index",
        "match_basis",
        "dr",
        "tight_tag_fixed",
        "rj_iso_upper",
        "ppg12_iso_upper",
        "rj_noniso_lower",
        "ppg12_noniso_lower",
        "rj_isolation_region",
        "ppg12_isolation_region",
        "isolation_region_flip",
        "rj_abcd_category",
        "ppg12_abcd_category",
        "abcd_category_flip",
        "abcd_letter_flip",
    ]
    for rj_name, ppg_name in FEATURES:
        fieldnames += [f"rj_{rj_name}", f"ppg12_{ppg_name}", f"delta_{rj_name}"]

    deltas_by_feature: dict[str, list[float]] = defaultdict(list)
    category_counts: dict[str, Counter[str]] = {"rj": Counter(), "ppg12": Counter()}
    category_transitions: Counter[tuple[str, str]] = Counter()
    region_transitions: Counter[tuple[str, str]] = Counter()
    isolation_region_flips = 0
    abcd_category_flips = 0
    abcd_letter_flips = 0
    match_basis_counts: Counter[str] = Counter()
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for rj, ppg, dist, match_basis in matches:
            row: dict[str, Any] = {
                "period": period,
                "tag": tag,
                "run": int(rj["run"]),
                "eventnumber": int(rj["eventnumber"]),
                "rj_entry": int(rj["source_entry"]),
                "ppg12_entry": int(ppg["source_entry"]),
                "rj_cluster_index": int(rj["cluster_index"]),
                "ppg12_cluster_index": int(ppg["cluster_index"]),
                "match_basis": match_basis,
                "dr": dist,
            }
            match_basis_counts[match_basis] += 1
            tight_tag = int(rj["ppg12_tight_tag"])
            rj_et = float(rj["cluster_Et"])
            ppg12_et = float(ppg["cluster_Et"])
            rj_raw_eiso = float(rj["ppg12_raw_eiso"])
            ppg12_raw_eiso = float(ppg["cluster_iso_topo_04"])
            rj_region = isolation_region(rj_raw_eiso, rj_et)
            ppg12_region = isolation_region(ppg12_raw_eiso, ppg12_et)
            rj_category = abcd_category(tight_tag, rj_raw_eiso, rj_et)
            ppg12_category = abcd_category(tight_tag, ppg12_raw_eiso, ppg12_et)
            region_flip = rj_region != ppg12_region
            category_flip = rj_category != ppg12_category
            letter_flip = (
                rj_category in ABCD_CATEGORIES
                and ppg12_category in ABCD_CATEGORIES
                and category_flip
            )
            row.update(
                {
                    "tight_tag_fixed": tight_tag,
                    "rj_iso_upper": 0.490 + 0.037 * rj_et,
                    "ppg12_iso_upper": 0.490 + 0.037 * ppg12_et,
                    "rj_noniso_lower": 1.290 + 0.037 * rj_et,
                    "ppg12_noniso_lower": 1.290 + 0.037 * ppg12_et,
                    "rj_isolation_region": rj_region,
                    "ppg12_isolation_region": ppg12_region,
                    "isolation_region_flip": int(region_flip),
                    "rj_abcd_category": rj_category,
                    "ppg12_abcd_category": ppg12_category,
                    "abcd_category_flip": int(category_flip),
                    "abcd_letter_flip": int(letter_flip),
                }
            )
            category_counts["rj"][rj_category] += 1
            category_counts["ppg12"][ppg12_category] += 1
            category_transitions[(rj_category, ppg12_category)] += 1
            region_transitions[(rj_region, ppg12_region)] += 1
            isolation_region_flips += int(region_flip)
            abcd_category_flips += int(category_flip)
            abcd_letter_flips += int(letter_flip)
            for rj_name, ppg_name in FEATURES:
                rj_val = float(rj.get(rj_name, float("nan")))
                ppg_val = float(ppg.get(ppg_name, float("nan")))
                delta = rj_val - ppg_val if finite(rj_val) and finite(ppg_val) else float("nan")
                row[f"rj_{rj_name}"] = rj_val
                row[f"ppg12_{ppg_name}"] = ppg_val
                row[f"delta_{rj_name}"] = delta
                if finite(delta):
                    deltas_by_feature[rj_name].append(delta)
            writer.writerow(row)

    eiso_deltas = deltas_by_feature.get("ppg12_raw_eiso", [])
    eiso_abs_deltas = [abs(value) for value in eiso_deltas]
    summary_row: dict[str, Any] = {
        "period": period,
        "tag": tag,
        "matched": len(matches),
        "unmatched_rj": len(misses),
        "mean_rj_minus_ppg12_raw_eiso": mean(eiso_deltas),
        "median_rj_minus_ppg12_raw_eiso": (
            statistics.median(eiso_deltas) if eiso_deltas else float("nan")
        ),
        "mean_abs_raw_eiso_delta": mean(eiso_abs_deltas),
        "max_abs_raw_eiso_delta": max(eiso_abs_deltas) if eiso_abs_deltas else float("nan"),
        "isolation_region_flips": isolation_region_flips,
        "abcd_category_flips": abcd_category_flips,
        "abcd_letter_flips": abcd_letter_flips,
        "match_basis_counts": ";".join(
            f"{basis}:{count}" for basis, count in sorted(match_basis_counts.items())
        ),
        "category_transitions": ";".join(
            f"{source}->{target}:{count}"
            for (source, target), count in sorted(category_transitions.items())
        ),
        "region_transitions": ";".join(
            f"{source}->{target}:{count}"
            for (source, target), count in sorted(region_transitions.items())
        ),
    }
    for side in ("rj", "ppg12"):
        for category in (*ABCD_CATEGORIES, "GAP", "OUTSIDE", "OTHER_ID", "INVALID"):
            summary_row[f"{side}_{category}"] = category_counts[side][category]
    with out_summary_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_row))
        writer.writeheader()
        writer.writerow(summary_row)

    lines = [
        "# THE-97 pp-data same-cluster raw-Eiso parity",
        "",
        f"- Period: `{period}`",
        f"- Canary tag: `{tag}`",
        f"- Matched clusters: {len(matches)}",
        f"- Unmatched RecoilJets rows: {len(misses)}",
        f"- Match bases: {dict(sorted(match_basis_counts.items()))}",
        f"- CSV: `{out_csv}`",
        f"- Summary CSV: `{out_summary_csv}`",
        "- ABCD flips hold the RecoilJets `ppg12_tight_tag` fixed; they isolate raw-Eiso/ET-window effects rather than photon-ID retagging.",
        "",
        "## Raw Eiso And ABCD Summary",
        "",
        f"- Mean RJ - PPG12 raw Eiso: {mean(eiso_deltas):+.6g} GeV",
        f"- Mean absolute raw-Eiso difference: {mean(eiso_abs_deltas):.6g} GeV",
        f"- Isolation-region flips: {isolation_region_flips}/{len(matches)}",
        f"- ABCD-category flips (including gap/outside): {abcd_category_flips}/{len(matches)}",
        f"- A/B/C/D letter-to-letter flips: {abcd_letter_flips}/{len(matches)}",
        "",
        "| source | A | B | C | D | gap | outside | other ID | invalid |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        f"| RecoilJets | {category_counts['rj']['A']} | {category_counts['rj']['B']} | {category_counts['rj']['C']} | {category_counts['rj']['D']} | {category_counts['rj']['GAP']} | {category_counts['rj']['OUTSIDE']} | {category_counts['rj']['OTHER_ID']} | {category_counts['rj']['INVALID']} |",
        f"| PPG12 Eiso | {category_counts['ppg12']['A']} | {category_counts['ppg12']['B']} | {category_counts['ppg12']['C']} | {category_counts['ppg12']['D']} | {category_counts['ppg12']['GAP']} | {category_counts['ppg12']['OUTSIDE']} | {category_counts['ppg12']['OTHER_ID']} | {category_counts['ppg12']['INVALID']} |",
        "",
        "### Category Transition Matrix",
        "",
        "| RJ category | PPG12 category | count |",
        "| --- | --- | ---: |",
    ]
    for (source, target), count in sorted(category_transitions.items()):
        lines.append(f"| {source} | {target} | {count} |")

    lines += [
        "",
        "## Mean Differences",
        "",
        "| feature | mean RJ - PPG12 | mean abs | max abs |",
        "| --- | ---: | ---: | ---: |",
    ]
    for rj_name, _ in FEATURES:
        vals = deltas_by_feature.get(rj_name, [])
        abs_vals = [abs(v) for v in vals if finite(v)]
        lines.append(
            f"| {rj_name} | {mean(vals):+.6g} | {mean(abs_vals):.6g} | "
            f"{(max(abs_vals) if abs_vals else float('nan')):.6g} |"
        )

    lines += [
        "",
        "## First Matched Rows",
        "",
        "| run | event | dR | RJ idx | PPG12 idx | RJ E11/E33 | PPG12 E11/E33 | delta | RJ weta | PPG12 weta | RJ et2/et3/et4 | PPG12 et2/et3/et4 |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    for rj, ppg, dist, _match_basis in matches[:25]:
        lines.append(
            f"| {int(rj['run'])} | {int(rj['eventnumber'])} | {dist:.5g} | "
            f"{int(rj['cluster_index'])} | {int(ppg['cluster_index'])} | "
            f"{float(rj['e11_over_e33']):.6g} | {float(ppg['e11_over_e33']):.6g} | "
            f"{float(rj['e11_over_e33']) - float(ppg['e11_over_e33']):+.3g} | "
            f"{float(rj['cluster_weta_cogx']):.6g} | {float(ppg['cluster_weta_cogx']):.6g} | "
            f"{float(rj['cluster_et2']):.4g}/{float(rj['cluster_et3']):.4g}/{float(rj['cluster_et4']):.4g} | "
            f"{float(ppg['cluster_et2']):.4g}/{float(ppg['cluster_et3']):.4g}/{float(ppg['cluster_et4']):.4g} |"
        )
    out_md.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    rj_tree = open_rj_tree(args.rj_root)
    ppg12_chain = open_ppg12_chain(args.ppg12_root_glob)
    rj_rows, wanted_events = read_rj_rows(rj_tree, args)
    if not rj_rows:
        raise RuntimeError("No RecoilJets rows survived the selected ET/eta cuts")
    ppg12_by_event, ppg12_by_run = read_ppg12_candidates(
        ppg12_chain, wanted_events, args
    )
    matches, misses = match_rows(
        rj_rows,
        ppg12_by_event,
        ppg12_by_run,
        args.max_dr,
        args.max_det,
    )
    out_csv = Path(args.out_csv)
    out_summary_csv = (
        Path(args.out_summary_csv)
        if args.out_summary_csv
        else out_csv.with_name(f"{out_csv.stem}_summary.csv")
    )
    write_outputs(
        matches,
        misses,
        out_csv,
        Path(args.out_md),
        out_summary_csv,
        args.period,
        args.tag,
    )
    print(
        f"rj_rows={len(rj_rows)} wanted_events={len(wanted_events)} "
        f"period={args.period} tag={args.tag} matched={len(matches)} unmatched={len(misses)} "
        f"out_md={args.out_md} out_summary_csv={out_summary_csv}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

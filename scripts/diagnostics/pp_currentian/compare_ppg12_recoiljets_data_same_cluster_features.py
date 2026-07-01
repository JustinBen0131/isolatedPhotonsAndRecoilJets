#!/usr/bin/env python3
"""Compare PPG12 data slimtree clusters with RecoilJets pp-data diagnostic rows.

THE-76 uses this after a tiny pp-data diagnostic canary with
RJ_PP_PHOTONID_TRAINING_TREE=1.  It is intentionally read-only: the script
matches clusters by run/event and nearest eta/phi, then reports whether
PhotonClusterBuilder/RecoilJets shower-shape inputs match PPG12 CaloAna24 for
the same pp data clusters.
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
from collections import defaultdict
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
]


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
    parser.add_argument("--max-rj-rows", type=int, default=1000)
    parser.add_argument("--max-ppg12-events", type=int, default=0, help="0 means scan all PPG12 entries")
    parser.add_argument("--et-min", type=float, default=22.0)
    parser.add_argument("--et-max", type=float, default=28.0)
    parser.add_argument("--eta-max", type=float, default=0.7)
    parser.add_argument("--max-dr", type=float, default=0.02)
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
    required = {"run", "eventnumber", "cluster_Et", "cluster_Eta", "cluster_Phi", "e11_over_e33"}
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
        }
        if row["run"] <= 0 or row["eventnumber"] <= 0:
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


def read_ppg12_by_event(
    chain: Any,
    wanted_events: set[tuple[int, int]],
    args: argparse.Namespace,
) -> dict[tuple[int, int], list[dict[str, float]]]:
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
    }
    missing = sorted(required - names)
    if missing:
        raise RuntimeError(f"PPG12 slimtree missing required branches: {missing}")

    out: dict[tuple[int, int], list[dict[str, float]]] = defaultdict(list)
    stop = chain.GetEntries()
    if args.max_ppg12_events > 0:
        stop = min(stop, args.max_ppg12_events)

    for entry in range(stop):
        chain.GetEntry(entry)
        run = get_int(chain, "runnumber")
        event = get_int(chain, "eventnumber")
        event_key = (run, event)
        if event_key not in wanted_events:
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
            }
            out[event_key].append(row)
    return out


def match_rows(
    rj_rows: list[dict[str, float]],
    ppg12_by_event: dict[tuple[int, int], list[dict[str, float]]],
    max_dr: float,
) -> tuple[list[tuple[dict[str, float], dict[str, float], float]], list[dict[str, float]]]:
    matches: list[tuple[dict[str, float], dict[str, float], float]] = []
    misses: list[dict[str, float]] = []
    used: set[tuple[int, int, int]] = set()
    for rj in rj_rows:
        event_key = (int(rj["run"]), int(rj["eventnumber"]))
        best_idx = -1
        best_dist = float("inf")
        best_row: dict[str, float] | None = None
        for idx, ppg in enumerate(ppg12_by_event.get(event_key, [])):
            key = (event_key[0], event_key[1], idx)
            if key in used:
                continue
            deta = float(rj["cluster_Eta"]) - float(ppg["cluster_Eta"])
            dph = dphi(float(rj["cluster_Phi"]), float(ppg["cluster_Phi"]))
            dist = math.hypot(deta, dph)
            if dist < best_dist:
                best_idx = idx
                best_dist = dist
                best_row = ppg
        if best_row is None or best_dist > max_dr:
            misses.append(rj)
            continue
        used.add((event_key[0], event_key[1], best_idx))
        matches.append((rj, best_row, best_dist))
    return matches, misses


def mean(values: list[float]) -> float:
    clean = [x for x in values if finite(x)]
    return sum(clean) / len(clean) if clean else float("nan")


def write_outputs(
    matches: list[tuple[dict[str, float], dict[str, float], float]],
    misses: list[dict[str, float]],
    out_csv: Path,
    out_md: Path,
) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_md.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "run",
        "eventnumber",
        "rj_entry",
        "ppg12_entry",
        "rj_cluster_index",
        "ppg12_cluster_index",
        "dr",
    ]
    for rj_name, ppg_name in FEATURES:
        fieldnames += [f"rj_{rj_name}", f"ppg12_{ppg_name}", f"delta_{rj_name}"]

    deltas_by_feature: dict[str, list[float]] = defaultdict(list)
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for rj, ppg, dist in matches:
            row: dict[str, Any] = {
                "run": int(rj["run"]),
                "eventnumber": int(rj["eventnumber"]),
                "rj_entry": int(rj["source_entry"]),
                "ppg12_entry": int(ppg["source_entry"]),
                "rj_cluster_index": int(rj["cluster_index"]),
                "ppg12_cluster_index": int(ppg["cluster_index"]),
                "dr": dist,
            }
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

    lines = [
        "# THE-76 pp-data same-cluster feature parity",
        "",
        f"- Matched clusters: {len(matches)}",
        f"- Unmatched RecoilJets rows: {len(misses)}",
        f"- CSV: `{out_csv}`",
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
    for rj, ppg, dist in matches[:25]:
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
    ppg12_by_event = read_ppg12_by_event(ppg12_chain, wanted_events, args)
    matches, misses = match_rows(rj_rows, ppg12_by_event, args.max_dr)
    write_outputs(matches, misses, Path(args.out_csv), Path(args.out_md))
    print(
        f"rj_rows={len(rj_rows)} wanted_events={len(wanted_events)} "
        f"matched={len(matches)} unmatched={len(misses)} out_md={args.out_md}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

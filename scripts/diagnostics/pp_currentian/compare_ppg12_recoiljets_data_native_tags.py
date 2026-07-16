#!/usr/bin/env python3
"""Read-only candidate/native-ABCD parity check for pp-data canaries.

The existing same-cluster data helper intentionally fixes the RecoilJets ID
tag when comparing isolation.  That is the wrong diagnostic for a shower-
shape canary, whose expected effect is to change the base-v3E inputs, scores,
and therefore native tight/non-tight membership.  This companion recomputes
the PPG12 common/tight/non-tight decision from each side's own features and
reports native A/B/C/D membership, including the 10--24 GeV bin counts.

The default output mode is safe for a streamed SDCC invocation: per-candidate
CSV goes to stdout and the JSON summary goes to stderr.  ROOT inputs are never
modified.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, TextIO


RECO_BINS = (10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0)
FEATURE_PAIRS = (
    ("cluster_Et", "cluster_Et"),
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
)
SHAPE_SCORE_FEATURES = (
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "e11_over_e33",
    "e32_over_e35",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "cluster_prob",
    "npb_score",
    "tight_bdt_score",
)


def finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except Exception:
        return False


def open_interval(value: float, low: float, high: float) -> bool:
    return finite(value) and low < value < high


def wrapped_dphi(lhs: float, rhs: float) -> float:
    value = lhs - rhs
    while value > math.pi:
        value -= 2.0 * math.pi
    while value <= -math.pi:
        value += 2.0 * math.pi
    return value


def thresholds(et: float) -> tuple[float, float, float]:
    return (
        0.815625 - 0.0015625 * et,
        0.7333333333333333 - 0.01333333333333333 * et,
        0.684375 + 0.0015625 * et,
    )


def classify(row: dict[str, float], score: float, et: float) -> tuple[int, bool]:
    """Return exact PPG12 tag (0 fail, 1 tight, 2 non-tight, 3 neither)."""
    e11e33 = row["e11_over_e33"]
    e32e35 = row["e32_over_e35"]
    weta = row["cluster_weta_cogx"]
    wphi = row["cluster_wphi_cogx"]
    wr = wphi / weta if weta != 0.0 else float("nan")
    common = (
        open_interval(row["cluster_prob"], 0.0, 1.0)
        and open_interval(e11e33, 0.0, 0.98)
        and finite(wr)
        and wr > 0.0
        and finite(weta)
        and weta < 2.0
        and row["npb_score"] > 0.5
    )
    if not common:
        return 0, False

    tight_min, nontight_min, nontight_max = thresholds(et)
    tight_prob = open_interval(row["cluster_prob"], 0.0, 1.0)
    tight_weta = open_interval(weta, 0.0, 1.0)
    tight_wphi = open_interval(wphi, 0.0, 1.0)
    tight_bdt = finite(score) and tight_min < score < 1.0
    tight = (
        tight_prob
        and tight_weta
        and tight_wphi
        and open_interval(row["cluster_et1"], 0.5, 1.0)
        and open_interval(row["cluster_et2"], 0.0, 1.0)
        and open_interval(row["cluster_et3"], 0.0, 1.0)
        and open_interval(row["cluster_et4"], 0.0, 1.0)
        and open_interval(e11e33, 0.0, 1.0)
        and open_interval(e32e35, 0.8, 1.0)
        and tight_bdt
    )
    if tight:
        return 1, True

    nontight_shape = (
        tight_prob
        and tight_weta
        and tight_wphi
        and open_interval(row["cluster_et1"], 0.6, 1.0)
        and open_interval(row["cluster_et4"], 0.0, 1.0)
        and open_interval(e11e33, 0.0, 1.0)
        and open_interval(e32e35, 0.8, 1.0)
    )
    nfail = int(not tight_weta) + int(not tight_prob) + int(not tight_bdt)
    if (
        nontight_shape
        and finite(score)
        and nontight_min < score < nontight_max
        and nfail > 0
    ):
        return 2, True
    return 3, True


def isolation_region(raw_eiso: float, et: float) -> str:
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


def native_category(tag: int, raw_eiso: float, et: float) -> str:
    region = isolation_region(raw_eiso, et)
    if tag == 1:
        return {"ISO": "A", "NONISO": "B"}.get(region, region)
    if tag == 2:
        return {"ISO": "C", "NONISO": "D"}.get(region, region)
    return "OTHER_ID" if region in ("ISO", "NONISO") else region


def find_bin(et: float) -> int:
    for index, (low, high) in enumerate(zip(RECO_BINS[:-1], RECO_BINS[1:])):
        if low < et < high:
            return index
    return -1


def branch_names(tree: Any) -> set[str]:
    return {branch.GetName() for branch in tree.GetListOfBranches()}


def scalar(tree: Any, name: str, default: float = float("nan")) -> float:
    try:
        return float(getattr(tree, name))
    except Exception:
        return default


def integer(tree: Any, name: str, default: int = -1) -> int:
    try:
        return int(getattr(tree, name))
    except Exception:
        return default


def array_value(tree: Any, name: str, index: int) -> float:
    try:
        return float(getattr(tree, name)[index])
    except Exception:
        return float("nan")


def open_rj_tree(ROOT: Any, pattern: str, tree_name: str) -> Any:
    paths = sorted(glob.glob(pattern))
    if not paths and Path(pattern).exists():
        paths = [pattern]
    if not paths:
        raise RuntimeError(f"No RecoilJets ROOT files match: {pattern}")
    if len(paths) == 1:
        root_file = ROOT.TFile.Open(paths[0])
        if not root_file or root_file.IsZombie():
            raise RuntimeError(f"Could not open RecoilJets ROOT: {paths[0]}")
        tree = root_file.Get(tree_name)
        if not tree:
            raise RuntimeError(f"Missing {tree_name}: {paths[0]}")
        tree._owner_file = root_file
        return tree
    chain = ROOT.TChain(tree_name)
    for path in paths:
        chain.Add(path)
    if chain.GetEntries() <= 0:
        raise RuntimeError(f"RecoilJets chain has no entries: {pattern}")
    return chain


def read_rj(tree: Any, args: argparse.Namespace) -> list[dict[str, float]]:
    required = {
        "run", "eventnumber", "cluster_Et", "cluster_Eta", "cluster_Phi",
        "vertexz", "cluster_weta_cogx", "cluster_wphi_cogx",
        "e11_over_e33", "e32_over_e35", "cluster_et1", "cluster_et2",
        "cluster_et3", "cluster_et4", "cluster_prob", "npb_score",
        "tight_bdt_score", "ppg12_common_pass", "ppg12_tight_tag",
        "ppg12_raw_eiso",
    }
    missing = sorted(required - branch_names(tree))
    if missing:
        raise RuntimeError(f"RecoilJets tree missing branches: {missing}")
    rows: list[dict[str, float]] = []
    for entry in range(tree.GetEntries()):
        tree.GetEntry(entry)
        et = scalar(tree, "cluster_Et")
        eta = scalar(tree, "cluster_Eta")
        if not (args.et_min < et < args.et_max and abs(eta) < args.eta_max):
            continue
        row = {
            "source_entry": entry,
            "run": integer(tree, "run"),
            "eventnumber": integer(tree, "eventnumber"),
            "cluster_index": integer(tree, "cluster_index"),
            "cluster_Et": et,
            "cluster_Eta": eta,
            "cluster_Phi": scalar(tree, "cluster_Phi"),
            "vertexz": scalar(tree, "vertexz"),
            "cluster_weta_cogx": scalar(tree, "cluster_weta_cogx"),
            "cluster_wphi_cogx": scalar(tree, "cluster_wphi_cogx"),
            "e11_over_e33": scalar(tree, "e11_over_e33"),
            "e32_over_e35": scalar(tree, "e32_over_e35"),
            "cluster_et1": scalar(tree, "cluster_et1"),
            "cluster_et2": scalar(tree, "cluster_et2"),
            "cluster_et3": scalar(tree, "cluster_et3"),
            "cluster_et4": scalar(tree, "cluster_et4"),
            "cluster_prob": scalar(tree, "cluster_prob"),
            "npb_score": scalar(tree, "npb_score"),
            "tight_bdt_score": scalar(tree, "tight_bdt_score"),
            "ppg12_common_pass": integer(tree, "ppg12_common_pass"),
            "ppg12_tight_tag": integer(tree, "ppg12_tight_tag"),
            "ppg12_raw_eiso": scalar(tree, "ppg12_raw_eiso"),
        }
        if row["run"] > 0:
            rows.append(row)
    return rows


def open_ppg_chain(ROOT: Any, pattern: str, wanted_runs: set[int]) -> Any:
    paths = sorted(glob.glob(pattern))
    if not paths and Path(pattern).exists():
        paths = [pattern]
    if not paths:
        raise RuntimeError(f"No PPG12 files match: {pattern}")
    # A full PPG12 data glob contains many large array trees.  First scan only
    # the scalar run branch so unrelated files never load cluster arrays.
    selected_paths: list[str] = []
    for path in paths:
        root_file = ROOT.TFile.Open(path)
        if not root_file or root_file.IsZombie():
            raise RuntimeError(f"Could not open PPG12 ROOT: {path}")
        tree = root_file.Get("slimtree")
        if not tree:
            root_file.Close()
            raise RuntimeError(f"Missing slimtree: {path}")
        if "runnumber" not in branch_names(tree):
            root_file.Close()
            raise RuntimeError(f"Missing runnumber: {path}")
        tree.SetBranchStatus("*", 0)
        tree.SetBranchStatus("runnumber", 1)
        run_selection = " || ".join(
            f"runnumber=={run}" for run in sorted(wanted_runs)
        )
        contains_wanted_run = tree.GetEntries(run_selection) > 0
        root_file.Close()
        if contains_wanted_run:
            selected_paths.append(path)
    if not selected_paths:
        raise RuntimeError(
            f"No PPG12 files contain requested runs {sorted(wanted_runs)}: {pattern}"
        )
    chain = ROOT.TChain("slimtree")
    for path in selected_paths:
        chain.Add(path)
    return chain


def read_ppg(
    chain: Any,
    wanted_runs: set[int],
    args: argparse.Namespace,
) -> tuple[dict[tuple[int, int], list[dict[str, float]]], dict[int, list[dict[str, float]]]]:
    node = args.node
    scalar_required = {"runnumber", "eventnumber", "vertexz", f"ncluster_{node}"}
    array_stems = (
        "cluster_Et", "cluster_Eta", "cluster_Phi", "cluster_e11",
        "cluster_e33", "cluster_e32", "cluster_e35", "cluster_weta_cogx",
        "cluster_wphi_cogx", "cluster_et1", "cluster_et2", "cluster_et3",
        "cluster_et4", "cluster_prob", "cluster_npb_score",
        "cluster_iso_topo_04",
    )
    arrays = {f"{stem}_{node}" for stem in array_stems}
    arrays.add(f"cluster_bdt_{node}_base_v3E")
    names = branch_names(chain)
    missing = sorted((scalar_required | arrays) - names)
    if missing:
        raise RuntimeError(f"PPG12 slimtree missing branches: {missing}")
    chain.SetBranchStatus("*", 0)
    for name in sorted(scalar_required | arrays):
        chain.SetBranchStatus(name, 1)

    by_event: dict[tuple[int, int], list[dict[str, float]]] = defaultdict(list)
    by_run: dict[int, list[dict[str, float]]] = defaultdict(list)
    for source_entry in range(chain.GetEntries()):
        chain.GetEntry(source_entry)
        run = integer(chain, "runnumber")
        if run not in wanted_runs:
            continue
        event = integer(chain, "eventnumber")
        vertexz = scalar(chain, "vertexz")
        ncluster = min(integer(chain, f"ncluster_{node}", 0), 20000)
        for cluster_index in range(ncluster):
            et = array_value(chain, f"cluster_Et_{node}", cluster_index)
            eta = array_value(chain, f"cluster_Eta_{node}", cluster_index)
            if not (args.et_min < et < args.et_max and abs(eta) < args.eta_max):
                continue
            e11 = array_value(chain, f"cluster_e11_{node}", cluster_index)
            e33 = array_value(chain, f"cluster_e33_{node}", cluster_index)
            e32 = array_value(chain, f"cluster_e32_{node}", cluster_index)
            e35 = array_value(chain, f"cluster_e35_{node}", cluster_index)
            row = {
                "source_entry": source_entry,
                "cluster_index": cluster_index,
                "run": run,
                "eventnumber": event,
                "cluster_Et": et,
                "cluster_Eta": eta,
                "cluster_Phi": array_value(chain, f"cluster_Phi_{node}", cluster_index),
                "vertexz": vertexz,
                "cluster_weta_cogx": array_value(chain, f"cluster_weta_cogx_{node}", cluster_index),
                "cluster_wphi_cogx": array_value(chain, f"cluster_wphi_cogx_{node}", cluster_index),
                "e11_over_e33": e11 / e33 if e33 > 0.0 else float("nan"),
                "e32_over_e35": e32 / e35 if e35 > 0.0 else float("nan"),
                "cluster_et1": array_value(chain, f"cluster_et1_{node}", cluster_index),
                "cluster_et2": array_value(chain, f"cluster_et2_{node}", cluster_index),
                "cluster_et3": array_value(chain, f"cluster_et3_{node}", cluster_index),
                "cluster_et4": array_value(chain, f"cluster_et4_{node}", cluster_index),
                "cluster_prob": array_value(chain, f"cluster_prob_{node}", cluster_index),
                "npb_score": array_value(chain, f"cluster_npb_score_{node}", cluster_index),
                "bdt_base_v3E": array_value(chain, f"cluster_bdt_{node}_base_v3E", cluster_index),
                "cluster_iso_topo_04": array_value(chain, f"cluster_iso_topo_04_{node}", cluster_index),
            }
            by_event[(run, event)].append(row)
            by_run[run].append(row)
    return by_event, by_run


def match_rows(
    rj_rows: list[dict[str, float]],
    by_event: dict[tuple[int, int], list[dict[str, float]]],
    by_run: dict[int, list[dict[str, float]]],
    max_dr: float,
    max_det: float,
) -> tuple[list[tuple[dict[str, float], dict[str, float], float, str]], list[dict[str, float]]]:
    matches: list[tuple[dict[str, float], dict[str, float], float, str]] = []
    misses: list[dict[str, float]] = []
    used: set[tuple[int, int]] = set()
    for rj in rj_rows:
        key = (int(rj["run"]), int(rj["eventnumber"]))
        candidates = by_event.get(key, [])
        basis = "event"
        if not candidates:
            candidates = by_run.get(key[0], [])
            basis = "run_kinematics"
        best: tuple[float, float, dict[str, float]] | None = None
        for ppg in candidates:
            ppg_key = (int(ppg["source_entry"]), int(ppg["cluster_index"]))
            if ppg_key in used:
                continue
            deta = rj["cluster_Eta"] - ppg["cluster_Eta"]
            dphi = wrapped_dphi(rj["cluster_Phi"], ppg["cluster_Phi"])
            dr = math.hypot(deta, dphi)
            det = abs(rj["cluster_Et"] - ppg["cluster_Et"])
            if dr > max_dr or det > max_det:
                continue
            score = (dr / max_dr) ** 2 + (det / max_det) ** 2
            if best is None or score < best[0]:
                best = (score, dr, ppg)
        if best is None:
            misses.append(rj)
            continue
        ppg = best[2]
        used.add((int(ppg["source_entry"]), int(ppg["cluster_index"])))
        matches.append((rj, ppg, best[1], basis))
    return matches, misses


def abs_delta(rj: dict[str, float], ppg: dict[str, float], rj_name: str, ppg_name: str) -> float:
    if rj_name == "cluster_Phi":
        return abs(wrapped_dphi(rj[rj_name], ppg[ppg_name]))
    return abs(rj[rj_name] - ppg[ppg_name])


def emit(
    matches: list[tuple[dict[str, float], dict[str, float], float, str]],
    misses: list[dict[str, float]],
    args: argparse.Namespace,
    csv_stream: TextIO,
) -> dict[str, Any]:
    fields = [
        "period", "stage", "run", "rj_entry", "ppg_entry", "match_basis", "dr",
        "rj_stored_common", "rj_recomputed_common", "ppg_recomputed_common",
        "rj_stored_tag", "rj_recomputed_tag", "ppg_recomputed_tag",
        "rj_region", "ppg_region", "rj_native_abcd", "ppg_native_abcd",
        "et_bin",
    ]
    for rj_name, ppg_name in FEATURE_PAIRS:
        fields.extend((f"rj_{rj_name}", f"ppg_{ppg_name}", f"abs_delta_{rj_name}"))
    writer = csv.DictWriter(csv_stream, fieldnames=fields)
    writer.writeheader()

    tag_mismatch = 0
    stored_recompute_mismatch = 0
    stored_common_recompute_mismatch = 0
    common_mismatch = 0
    category_mismatch = 0
    basis_counts: Counter[str] = Counter()
    rj_counts: Counter[tuple[int, str]] = Counter()
    ppg_counts: Counter[tuple[int, str]] = Counter()
    max_abs: dict[str, float] = defaultdict(float)
    for rj, ppg, dr, basis in matches:
        rj_tag, rj_common = classify(rj, rj["tight_bdt_score"], rj["cluster_Et"])
        ppg_tag, ppg_common = classify(ppg, ppg["bdt_base_v3E"], ppg["cluster_Et"])
        rj_region = isolation_region(rj["ppg12_raw_eiso"], rj["cluster_Et"])
        ppg_region = isolation_region(ppg["cluster_iso_topo_04"], ppg["cluster_Et"])
        rj_category = native_category(rj_tag, rj["ppg12_raw_eiso"], rj["cluster_Et"])
        ppg_category = native_category(ppg_tag, ppg["cluster_iso_topo_04"], ppg["cluster_Et"])
        bin_index = find_bin(ppg["cluster_Et"])
        if bin_index >= 0:
            rj_counts[(bin_index, rj_category)] += 1
            ppg_counts[(bin_index, ppg_category)] += 1
        tag_mismatch += int(rj_tag != ppg_tag)
        stored_recompute_mismatch += int(int(rj["ppg12_tight_tag"]) != rj_tag)
        stored_common_recompute_mismatch += int(
            int(rj["ppg12_common_pass"]) != int(rj_common)
        )
        common_mismatch += int(rj_common != ppg_common)
        category_mismatch += int(rj_category != ppg_category)
        basis_counts[basis] += 1
        output: dict[str, Any] = {
            "period": args.period,
            "stage": args.stage,
            "run": int(rj["run"]),
            "rj_entry": int(rj["source_entry"]),
            "ppg_entry": int(ppg["source_entry"]),
            "match_basis": basis,
            "dr": dr,
            "rj_stored_common": int(rj["ppg12_common_pass"]),
            "rj_recomputed_common": int(rj_common),
            "ppg_recomputed_common": int(ppg_common),
            "rj_stored_tag": int(rj["ppg12_tight_tag"]),
            "rj_recomputed_tag": rj_tag,
            "ppg_recomputed_tag": ppg_tag,
            "rj_region": rj_region,
            "ppg_region": ppg_region,
            "rj_native_abcd": rj_category,
            "ppg_native_abcd": ppg_category,
            "et_bin": (
                f"{RECO_BINS[bin_index]:g}-{RECO_BINS[bin_index + 1]:g}"
                if bin_index >= 0 else "OUTSIDE"
            ),
        }
        for rj_name, ppg_name in FEATURE_PAIRS:
            delta = abs_delta(rj, ppg, rj_name, ppg_name)
            output[f"rj_{rj_name}"] = rj[rj_name]
            output[f"ppg_{ppg_name}"] = ppg[ppg_name]
            output[f"abs_delta_{rj_name}"] = delta
            if finite(delta):
                max_abs[rj_name] = max(max_abs[rj_name], delta)
        writer.writerow(output)

    bins: list[dict[str, Any]] = []
    for index in range(len(RECO_BINS) - 1):
        row: dict[str, Any] = {
            "bin": f"{RECO_BINS[index]:g}-{RECO_BINS[index + 1]:g}"
        }
        for side, counts in (("rj", rj_counts), ("ppg", ppg_counts)):
            for category in "ABCD":
                row[f"{side}_{category}"] = counts[(index, category)]
        bins.append(row)

    kinematic_max = max(
        max_abs.get("cluster_Et", 0.0),
        max_abs.get("cluster_Eta", 0.0),
        max_abs.get("cluster_Phi", 0.0),
    )
    shape_score_max = max((max_abs.get(name, 0.0) for name in SHAPE_SCORE_FEATURES), default=0.0)
    failures: list[str] = []
    if misses:
        failures.append(f"unmatched_rj={len(misses)}")
    if kinematic_max > args.kinematic_tolerance:
        failures.append(f"kinematic_max={kinematic_max:.9g}")
    if max_abs.get("ppg12_raw_eiso", 0.0) > args.eiso_tolerance:
        failures.append(f"eiso_max={max_abs['ppg12_raw_eiso']:.9g}")
    if shape_score_max > args.feature_tolerance:
        failures.append(f"shape_score_max={shape_score_max:.9g}")
    if stored_recompute_mismatch:
        failures.append(f"stored_recompute_tag_mismatch={stored_recompute_mismatch}")
    if stored_common_recompute_mismatch:
        failures.append(
            "stored_recompute_common_mismatch="
            f"{stored_common_recompute_mismatch}"
        )
    if common_mismatch:
        failures.append(f"common_mismatch={common_mismatch}")
    if tag_mismatch:
        failures.append(f"native_tag_mismatch={tag_mismatch}")
    if category_mismatch:
        failures.append(f"native_abcd_mismatch={category_mismatch}")

    return {
        "period": args.period,
        "stage": args.stage,
        "rj_root": args.rj_root,
        "ppg12_root_glob": args.ppg12_root_glob,
        "selected_rj": len(matches) + len(misses),
        "matched": len(matches),
        "unmatched_rj": len(misses),
        "match_basis": dict(sorted(basis_counts.items())),
        "rj_stored_vs_recomputed_tag_mismatches": stored_recompute_mismatch,
        "rj_stored_vs_recomputed_common_mismatches": stored_common_recompute_mismatch,
        "rj_vs_ppg_common_mismatches": common_mismatch,
        "rj_vs_ppg_native_tag_mismatches": tag_mismatch,
        "rj_vs_ppg_native_abcd_mismatches": category_mismatch,
        "ppg_b_candidates": sum(ppg_counts[(i, "B")] for i in range(len(RECO_BINS) - 1)),
        "ppg_d_candidates": sum(ppg_counts[(i, "D")] for i in range(len(RECO_BINS) - 1)),
        "max_abs_delta": dict(sorted(max_abs.items())),
        "binwise_abcd": bins,
        "closure_pass": not failures,
        "closure_failures": failures,
    }


def self_test() -> None:
    base = {
        "cluster_weta_cogx": 0.2,
        "cluster_wphi_cogx": 0.2,
        "e11_over_e33": 0.9,
        "e32_over_e35": 0.9,
        "cluster_et1": 0.8,
        "cluster_et2": 0.5,
        "cluster_et3": 0.5,
        "cluster_et4": 0.5,
        "cluster_prob": 0.8,
        "npb_score": 0.8,
    }
    assert classify(base, 0.9, 15.0) == (1, True)
    assert classify(base, 0.6, 15.0) == (2, True)
    assert classify(base, 0.75, 15.0) == (3, True)
    failed = dict(base, npb_score=0.4)
    assert classify(failed, 0.9, 15.0) == (0, False)
    upper = 0.490 + 0.037 * 15.0
    assert isolation_region(upper - 1e-6, 15.0) == "ISO"
    assert isolation_region(upper, 15.0) == "GAP"
    assert isolation_region(upper + 0.8, 15.0) == "GAP"
    assert isolation_region(upper + 0.800001, 15.0) == "NONISO"
    assert native_category(1, 0.0, 15.0) == "A"
    assert native_category(2, 3.0, 15.0) == "D"
    assert find_bin(11.0) == 0 and find_bin(23.0) == 6 and find_bin(12.0) == -1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rj-root")
    parser.add_argument("--ppg12-root-glob")
    parser.add_argument("--tree", default="AuAuPhotonIDTrainingTree")
    parser.add_argument("--node", default="CLUSTERINFO_CEMC")
    parser.add_argument("--period", default="unknown")
    parser.add_argument("--stage", default="unknown")
    parser.add_argument("--et-min", type=float, default=10.0)
    parser.add_argument("--et-max", type=float, default=24.0)
    parser.add_argument("--eta-max", type=float, default=0.7)
    parser.add_argument("--max-dr", type=float, default=1e-5)
    parser.add_argument("--max-det", type=float, default=1e-4)
    parser.add_argument("--kinematic-tolerance", type=float, default=2e-6)
    parser.add_argument("--eiso-tolerance", type=float, default=2e-6)
    parser.add_argument("--feature-tolerance", type=float, default=5e-6)
    parser.add_argument("--out-csv", default="-", help="- writes CSV to stdout")
    parser.add_argument("--out-summary-json", default="-", help="- writes JSON to stderr")
    parser.add_argument("--assert-closure", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def output_stream(path: str, default: TextIO) -> tuple[TextIO, bool]:
    if path == "-":
        return default, False
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output.open("w"), True


def main() -> int:
    args = parse_args()
    if args.self_test:
        self_test()
        print("native_tag_self_test=PASS")
        return 0
    if not args.rj_root or not args.ppg12_root_glob:
        raise SystemExit("--rj-root and --ppg12-root-glob are required")
    import ROOT  # type: ignore

    ROOT.gROOT.SetBatch(True)
    rj_tree = open_rj_tree(ROOT, args.rj_root, args.tree)
    rj_rows = read_rj(rj_tree, args)
    wanted_runs = {int(row["run"]) for row in rj_rows}
    ppg_chain = open_ppg_chain(ROOT, args.ppg12_root_glob, wanted_runs)
    by_event, by_run = read_ppg(ppg_chain, wanted_runs, args)
    matches, misses = match_rows(rj_rows, by_event, by_run, args.max_dr, args.max_det)

    csv_stream, close_csv = output_stream(args.out_csv, sys.stdout)
    try:
        summary = emit(matches, misses, args, csv_stream)
    finally:
        if close_csv:
            csv_stream.close()
    summary_stream, close_summary = output_stream(args.out_summary_json, sys.stderr)
    try:
        json.dump(summary, summary_stream, indent=2, sort_keys=True)
        summary_stream.write("\n")
    finally:
        if close_summary:
            summary_stream.close()
    if args.assert_closure and not summary["closure_pass"]:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

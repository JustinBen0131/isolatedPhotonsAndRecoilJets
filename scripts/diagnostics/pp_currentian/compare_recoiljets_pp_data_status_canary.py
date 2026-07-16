#!/usr/bin/env python3
"""Compare two RecoilJets pp-data reconstruction canaries candidate by candidate.

This read-only diagnostic is intended for one-variable reconstruction A/B tests,
especially the THE-97 tower-status OFF/ON canary.  It matches
``AuAuPhotonIDTrainingTree`` candidates one-to-one within each run/event by
nearest eta/phi, then reports status-ON minus status-OFF changes in candidate
kinematics, raw PPG12 topo isolation, isolation-region membership, and ABCD
membership while holding the status-OFF tight/non-tight tag fixed.

The PPG12 raw-isolation windows are evaluated with their strict inequalities:

* ISO: ``-20 < Eiso < 0.490 + 0.037*ET``
* GAP: the closed 0.8 GeV interval between ISO and NONISO
* NONISO: ``0.490 + 0.037*ET + 0.8 < Eiso < 20``

No ROOT file is modified.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

import ROOT  # type: ignore


ROOT.gROOT.SetBatch(True)

ABCD = ("A", "B", "C", "D")
REGIONS = ("ISO", "GAP", "NONISO", "OUTSIDE", "INVALID")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--status-off-root", required=True)
    parser.add_argument("--status-on-root", required=True)
    parser.add_argument("--out-csv", required=True, help="Per-candidate CSV")
    parser.add_argument("--out-summary-csv", required=True, help="One-row summary CSV")
    parser.add_argument("--out-md", required=True, help="Human-readable summary")
    parser.add_argument("--tree", default="AuAuPhotonIDTrainingTree")
    parser.add_argument("--off-label", default="status-OFF")
    parser.add_argument("--on-label", default="status-ON")
    parser.add_argument("--period", default="unknown")
    parser.add_argument("--tag", default="unknown")
    parser.add_argument("--et-min", type=float, default=10.0)
    parser.add_argument("--et-max", type=float, default=24.0)
    parser.add_argument("--eta-max", type=float, default=0.7)
    parser.add_argument("--max-dr", type=float, default=0.02)
    parser.add_argument(
        "--max-rows-per-input",
        type=int,
        default=0,
        help="0 scans every tree row; intended only as a diagnostic throttle",
    )
    return parser.parse_args()


def finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except Exception:
        return False


def wrapped_dphi(lhs: float, rhs: float) -> float:
    value = lhs - rhs
    while value > math.pi:
        value -= 2.0 * math.pi
    while value <= -math.pi:
        value += 2.0 * math.pi
    return value


def isolation_region(raw_eiso: float, et: float) -> str:
    """Classify the raw-Eiso window with exact PPG12 strict boundaries."""
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
    """Classify A/B/C/D using a supplied, possibly fixed, ID tag."""
    region = isolation_region(raw_eiso, et)
    if tight_tag == 1:
        return {"ISO": "A", "NONISO": "B"}.get(region, region)
    if tight_tag == 2:
        return {"ISO": "C", "NONISO": "D"}.get(region, region)
    return "OTHER_ID" if region in ("ISO", "NONISO") else region


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


def branch_names(tree: Any) -> set[str]:
    return {branch.GetName() for branch in tree.GetListOfBranches()}


def open_tree(path: str, tree_name: str) -> Any:
    root_file = ROOT.TFile.Open(path)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open ROOT file: {path}")
    tree = root_file.Get(tree_name)
    if not tree:
        raise RuntimeError(f"Missing {tree_name} in {path}")
    tree._owner_file = root_file
    return tree


def read_rows(
    tree: Any,
    *,
    et_min: float,
    et_max: float,
    eta_max: float,
    max_rows: int,
) -> tuple[list[dict[str, Any]], int]:
    names = branch_names(tree)
    required = {
        "run",
        "eventnumber",
        "cluster_Et",
        "cluster_Eta",
        "cluster_Phi",
        "ppg12_raw_eiso",
        "ppg12_tight_tag",
    }
    missing = sorted(required - names)
    if missing:
        raise RuntimeError(f"Tree is missing required branches: {missing}")

    rows: list[dict[str, Any]] = []
    total_entries = int(tree.GetEntries())
    stop = total_entries if max_rows <= 0 else min(total_entries, max_rows)
    for entry in range(stop):
        tree.GetEntry(entry)
        et = get_scalar(tree, "cluster_Et")
        eta = get_scalar(tree, "cluster_Eta")
        if not (finite(et) and finite(eta)):
            continue
        if not (et_min < et < et_max and abs(eta) < eta_max):
            continue
        run = get_int(tree, "run")
        event = get_int(tree, "eventnumber")
        if run <= 0 or event <= 0:
            continue
        raw_eiso = get_scalar(tree, "ppg12_raw_eiso")
        tight_tag = get_int(tree, "ppg12_tight_tag")
        rows.append(
            {
                "entry": entry,
                "run": run,
                "eventnumber": event,
                "cluster_index": get_int(tree, "cluster_index"),
                "et": et,
                "eta": eta,
                "phi": get_scalar(tree, "cluster_Phi"),
                "raw_eiso": raw_eiso,
                "tight_tag": tight_tag,
                "common_pass": get_int(tree, "ppg12_common_pass"),
                "region": isolation_region(raw_eiso, et),
                "native_abcd": abcd_category(tight_tag, raw_eiso, et),
            }
        )
    return rows, total_entries


def group_by_event(rows: Iterable[dict[str, Any]]) -> dict[tuple[int, int], list[dict[str, Any]]]:
    grouped: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[(int(row["run"]), int(row["eventnumber"]))].append(row)
    return grouped


def match_rows(
    off_rows: list[dict[str, Any]],
    on_rows: list[dict[str, Any]],
    max_dr: float,
) -> tuple[
    list[tuple[dict[str, Any], dict[str, Any], float]],
    list[dict[str, Any]],
    list[dict[str, Any]],
]:
    """One-to-one nearest eta/phi matching within identical run/event keys."""
    off_by_event = group_by_event(off_rows)
    on_by_event = group_by_event(on_rows)
    matches: list[tuple[dict[str, Any], dict[str, Any], float]] = []
    unmatched_off: list[dict[str, Any]] = []
    unmatched_on: list[dict[str, Any]] = []

    for event_key in sorted(set(off_by_event) | set(on_by_event)):
        offs = off_by_event.get(event_key, [])
        ons = on_by_event.get(event_key, [])
        # Sort all possible edges, then take the globally closest unused pair.
        # This is deterministic and avoids source-order-dependent nearest-neighbor
        # reuse when an event has multiple nearby clusters.
        edges: list[tuple[float, int, int]] = []
        for off_index, off in enumerate(offs):
            for on_index, on in enumerate(ons):
                deta = float(on["eta"]) - float(off["eta"])
                dph = wrapped_dphi(float(on["phi"]), float(off["phi"]))
                edges.append((math.hypot(deta, dph), off_index, on_index))
        edges.sort()
        used_off: set[int] = set()
        used_on: set[int] = set()
        for dr, off_index, on_index in edges:
            if dr > max_dr:
                break
            if off_index in used_off or on_index in used_on:
                continue
            used_off.add(off_index)
            used_on.add(on_index)
            matches.append((offs[off_index], ons[on_index], dr))
        unmatched_off.extend(row for index, row in enumerate(offs) if index not in used_off)
        unmatched_on.extend(row for index, row in enumerate(ons) if index not in used_on)
    return matches, unmatched_off, unmatched_on


def clean(values: Iterable[float]) -> list[float]:
    return [float(value) for value in values if finite(value)]


def mean(values: Iterable[float]) -> float:
    valid = clean(values)
    return sum(valid) / len(valid) if valid else float("nan")


def median(values: Iterable[float]) -> float:
    valid = clean(values)
    return statistics.median(valid) if valid else float("nan")


def maximum_abs(values: Iterable[float]) -> float:
    valid = [abs(value) for value in clean(values)]
    return max(valid) if valid else float("nan")


def transition_text(counter: Counter[tuple[str, str]]) -> str:
    return ";".join(
        f"{source}->{target}:{count}"
        for (source, target), count in sorted(counter.items())
    )


def fmt(value: float) -> str:
    return f"{value:+.8g}" if finite(value) else "nan"


def write_outputs(
    *,
    args: argparse.Namespace,
    off_total_entries: int,
    on_total_entries: int,
    off_rows: list[dict[str, Any]],
    on_rows: list[dict[str, Any]],
    matches: list[tuple[dict[str, Any], dict[str, Any], float]],
    unmatched_off: list[dict[str, Any]],
    unmatched_on: list[dict[str, Any]],
) -> None:
    out_csv = Path(args.out_csv)
    out_summary_csv = Path(args.out_summary_csv)
    out_md = Path(args.out_md)
    for path in (out_csv, out_summary_csv, out_md):
        path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "match_status",
        "period",
        "tag",
        "run",
        "eventnumber",
        "off_entry",
        "on_entry",
        "off_cluster_index",
        "on_cluster_index",
        "dr",
        "off_et",
        "on_et",
        "delta_et_on_minus_off",
        "off_eta",
        "on_eta",
        "delta_eta_on_minus_off",
        "off_phi",
        "on_phi",
        "delta_phi_on_minus_off",
        "off_raw_eiso",
        "on_raw_eiso",
        "delta_raw_eiso_on_minus_off",
        "off_common_pass",
        "on_common_pass",
        "off_tight_tag",
        "on_tight_tag",
        "tight_tag_changed",
        "off_region",
        "on_region",
        "region_transition",
        "fixed_off_tag",
        "off_fixed_abcd",
        "on_fixed_abcd",
        "fixed_abcd_transition",
        "fixed_abcd_changed",
        "off_native_abcd",
        "on_native_abcd",
        "native_abcd_transition",
    ]

    delta_et: list[float] = []
    delta_eta: list[float] = []
    delta_phi: list[float] = []
    delta_eiso: list[float] = []
    region_transitions: Counter[tuple[str, str]] = Counter()
    fixed_abcd_transitions: Counter[tuple[str, str]] = Counter()
    native_abcd_transitions: Counter[tuple[str, str]] = Counter()
    off_region_counts: Counter[str] = Counter()
    on_region_counts: Counter[str] = Counter()
    off_abcd_counts: Counter[str] = Counter()
    on_fixed_abcd_counts: Counter[str] = Counter()
    tight_tag_changes = 0
    common_pass_changes = 0
    region_changes = 0
    fixed_abcd_changes = 0

    def unmatched_row(side: str, row: dict[str, Any]) -> dict[str, Any]:
        result: dict[str, Any] = {name: "" for name in fieldnames}
        result.update(
            {
                "match_status": f"unmatched_{side}",
                "period": args.period,
                "tag": args.tag,
                "run": row["run"],
                "eventnumber": row["eventnumber"],
                f"{side}_entry": row["entry"],
                f"{side}_cluster_index": row["cluster_index"],
                f"{side}_et": row["et"],
                f"{side}_eta": row["eta"],
                f"{side}_phi": row["phi"],
                f"{side}_raw_eiso": row["raw_eiso"],
                f"{side}_common_pass": row["common_pass"],
                f"{side}_tight_tag": row["tight_tag"],
                f"{side}_region": row["region"],
                f"{side}_native_abcd": row["native_abcd"],
            }
        )
        return result

    with out_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for off, on, dr in matches:
            d_et = float(on["et"]) - float(off["et"])
            d_eta = float(on["eta"]) - float(off["eta"])
            d_phi = wrapped_dphi(float(on["phi"]), float(off["phi"]))
            d_eiso = float(on["raw_eiso"]) - float(off["raw_eiso"])
            fixed_tag = int(off["tight_tag"])
            off_fixed_abcd = abcd_category(fixed_tag, float(off["raw_eiso"]), float(off["et"]))
            on_fixed_abcd = abcd_category(fixed_tag, float(on["raw_eiso"]), float(on["et"]))

            off_region = str(off["region"])
            on_region = str(on["region"])
            off_native = str(off["native_abcd"])
            on_native = str(on["native_abcd"])
            region_transitions[(off_region, on_region)] += 1
            fixed_abcd_transitions[(off_fixed_abcd, on_fixed_abcd)] += 1
            native_abcd_transitions[(off_native, on_native)] += 1
            off_region_counts[off_region] += 1
            on_region_counts[on_region] += 1
            off_abcd_counts[off_fixed_abcd] += 1
            on_fixed_abcd_counts[on_fixed_abcd] += 1
            tight_tag_changes += int(int(off["tight_tag"]) != int(on["tight_tag"]))
            common_pass_changes += int(int(off["common_pass"]) != int(on["common_pass"]))
            region_changes += int(off_region != on_region)
            fixed_abcd_changes += int(off_fixed_abcd != on_fixed_abcd)
            delta_et.append(d_et)
            delta_eta.append(d_eta)
            delta_phi.append(d_phi)
            delta_eiso.append(d_eiso)

            writer.writerow(
                {
                    "match_status": "matched",
                    "period": args.period,
                    "tag": args.tag,
                    "run": off["run"],
                    "eventnumber": off["eventnumber"],
                    "off_entry": off["entry"],
                    "on_entry": on["entry"],
                    "off_cluster_index": off["cluster_index"],
                    "on_cluster_index": on["cluster_index"],
                    "dr": dr,
                    "off_et": off["et"],
                    "on_et": on["et"],
                    "delta_et_on_minus_off": d_et,
                    "off_eta": off["eta"],
                    "on_eta": on["eta"],
                    "delta_eta_on_minus_off": d_eta,
                    "off_phi": off["phi"],
                    "on_phi": on["phi"],
                    "delta_phi_on_minus_off": d_phi,
                    "off_raw_eiso": off["raw_eiso"],
                    "on_raw_eiso": on["raw_eiso"],
                    "delta_raw_eiso_on_minus_off": d_eiso,
                    "off_common_pass": off["common_pass"],
                    "on_common_pass": on["common_pass"],
                    "off_tight_tag": off["tight_tag"],
                    "on_tight_tag": on["tight_tag"],
                    "tight_tag_changed": int(int(off["tight_tag"]) != int(on["tight_tag"])),
                    "off_region": off_region,
                    "on_region": on_region,
                    "region_transition": f"{off_region}->{on_region}",
                    "fixed_off_tag": fixed_tag,
                    "off_fixed_abcd": off_fixed_abcd,
                    "on_fixed_abcd": on_fixed_abcd,
                    "fixed_abcd_transition": f"{off_fixed_abcd}->{on_fixed_abcd}",
                    "fixed_abcd_changed": int(off_fixed_abcd != on_fixed_abcd),
                    "off_native_abcd": off_native,
                    "on_native_abcd": on_native,
                    "native_abcd_transition": f"{off_native}->{on_native}",
                }
            )
        for row in unmatched_off:
            writer.writerow(unmatched_row("off", row))
        for row in unmatched_on:
            writer.writerow(unmatched_row("on", row))

    summary: dict[str, Any] = {
        "period": args.period,
        "tag": args.tag,
        "off_label": args.off_label,
        "on_label": args.on_label,
        "status_off_root": args.status_off_root,
        "status_on_root": args.status_on_root,
        "tree": args.tree,
        "off_total_tree_entries": off_total_entries,
        "on_total_tree_entries": on_total_entries,
        "off_selected_candidates": len(off_rows),
        "on_selected_candidates": len(on_rows),
        "matched_candidates": len(matches),
        "unmatched_off_candidates": len(unmatched_off),
        "unmatched_on_candidates": len(unmatched_on),
        "tight_tag_changes": tight_tag_changes,
        "common_pass_changes": common_pass_changes,
        "isolation_region_changes": region_changes,
        "fixed_abcd_changes": fixed_abcd_changes,
        "mean_delta_et_on_minus_off": mean(delta_et),
        "median_delta_et_on_minus_off": median(delta_et),
        "max_abs_delta_et": maximum_abs(delta_et),
        "mean_delta_eta_on_minus_off": mean(delta_eta),
        "max_abs_delta_eta": maximum_abs(delta_eta),
        "mean_delta_phi_on_minus_off": mean(delta_phi),
        "max_abs_delta_phi": maximum_abs(delta_phi),
        "mean_delta_raw_eiso_on_minus_off": mean(delta_eiso),
        "median_delta_raw_eiso_on_minus_off": median(delta_eiso),
        "mean_abs_delta_raw_eiso": mean(abs(value) for value in delta_eiso),
        "max_abs_delta_raw_eiso": maximum_abs(delta_eiso),
        "region_transitions": transition_text(region_transitions),
        "fixed_abcd_transitions": transition_text(fixed_abcd_transitions),
        "native_abcd_transitions": transition_text(native_abcd_transitions),
    }
    for region in REGIONS:
        summary[f"off_region_{region}"] = off_region_counts[region]
        summary[f"on_region_{region}"] = on_region_counts[region]
    for category in (*ABCD, "GAP", "OUTSIDE", "OTHER_ID", "INVALID"):
        summary[f"off_fixed_{category}"] = off_abcd_counts[category]
        summary[f"on_fixed_{category}"] = on_fixed_abcd_counts[category]

    with out_summary_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary))
        writer.writeheader()
        writer.writerow(summary)

    lines = [
        "# RecoilJets pp-data tower-status A/B candidate comparison",
        "",
        f"- Period: `{args.period}`",
        f"- Tag: `{args.tag}`",
        f"- {args.off_label} input: `{args.status_off_root}`",
        f"- {args.on_label} input: `{args.status_on_root}`",
        f"- Tree: `{args.tree}`",
        f"- Selection: `{args.et_min} < ET < {args.et_max} GeV`, `|eta| < {args.eta_max}`",
        f"- Match: identical run/event, globally nearest unused eta/phi pair with `dR <= {args.max_dr}`",
        "- Direction for every delta: status-ON minus status-OFF.",
        "- Fixed ABCD transitions hold the status-OFF `ppg12_tight_tag` fixed, so they isolate ET/raw-Eiso effects.",
        "",
        "## Candidate accounting",
        "",
        "| quantity | count |",
        "| --- | ---: |",
        f"| OFF total tree entries | {off_total_entries} |",
        f"| ON total tree entries | {on_total_entries} |",
        f"| OFF selected candidates | {len(off_rows)} |",
        f"| ON selected candidates | {len(on_rows)} |",
        f"| matched candidates | {len(matches)} |",
        f"| unmatched OFF | {len(unmatched_off)} |",
        f"| unmatched ON | {len(unmatched_on)} |",
        "",
        "## Matched-candidate deltas",
        "",
        "| quantity | mean | median | maximum absolute |",
        "| --- | ---: | ---: | ---: |",
        f"| ET [GeV] | {fmt(mean(delta_et))} | {fmt(median(delta_et))} | {fmt(maximum_abs(delta_et))} |",
        f"| eta | {fmt(mean(delta_eta))} | {fmt(median(delta_eta))} | {fmt(maximum_abs(delta_eta))} |",
        f"| phi | {fmt(mean(delta_phi))} | {fmt(median(delta_phi))} | {fmt(maximum_abs(delta_phi))} |",
        f"| raw Eiso [GeV] | {fmt(mean(delta_eiso))} | {fmt(median(delta_eiso))} | {fmt(maximum_abs(delta_eiso))} |",
        "",
        f"- Tight/non-tight tag changes: {tight_tag_changes}/{len(matches)}",
        f"- PPG12 common-pass changes: {common_pass_changes}/{len(matches)}",
        f"- Isolation-region changes: {region_changes}/{len(matches)}",
        f"- Fixed-tag ABCD/category changes: {fixed_abcd_changes}/{len(matches)}",
        "",
        "## Isolation-region transitions",
        "",
        "| OFF | ON | count |",
        "| --- | --- | ---: |",
    ]
    lines.extend(
        f"| {source} | {target} | {count} |"
        for (source, target), count in sorted(region_transitions.items())
    )
    lines += [
        "",
        "## Fixed-tag ABCD/category transitions",
        "",
        "| OFF | ON with OFF tag fixed | count |",
        "| --- | --- | ---: |",
    ]
    lines.extend(
        f"| {source} | {target} | {count} |"
        for (source, target), count in sorted(fixed_abcd_transitions.items())
    )
    lines += [
        "",
        "## Output files",
        "",
        f"- Per-candidate CSV: `{out_csv}`",
        f"- Summary CSV: `{out_summary_csv}`",
    ]
    out_md.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    if not (args.et_min < args.et_max):
        raise ValueError("--et-min must be less than --et-max")
    if args.eta_max <= 0 or args.max_dr <= 0:
        raise ValueError("--eta-max and --max-dr must be positive")

    off_tree = open_tree(args.status_off_root, args.tree)
    on_tree = open_tree(args.status_on_root, args.tree)
    off_rows, off_total_entries = read_rows(
        off_tree,
        et_min=args.et_min,
        et_max=args.et_max,
        eta_max=args.eta_max,
        max_rows=args.max_rows_per_input,
    )
    on_rows, on_total_entries = read_rows(
        on_tree,
        et_min=args.et_min,
        et_max=args.et_max,
        eta_max=args.eta_max,
        max_rows=args.max_rows_per_input,
    )
    matches, unmatched_off, unmatched_on = match_rows(off_rows, on_rows, args.max_dr)
    write_outputs(
        args=args,
        off_total_entries=off_total_entries,
        on_total_entries=on_total_entries,
        off_rows=off_rows,
        on_rows=on_rows,
        matches=matches,
        unmatched_off=unmatched_off,
        unmatched_on=unmatched_on,
    )
    print(
        f"off_selected={len(off_rows)} on_selected={len(on_rows)} "
        f"matched={len(matches)} unmatched_off={len(unmatched_off)} "
        f"unmatched_on={len(unmatched_on)} out_md={args.out_md}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

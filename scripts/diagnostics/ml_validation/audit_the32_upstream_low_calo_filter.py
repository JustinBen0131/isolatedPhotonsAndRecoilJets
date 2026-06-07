#!/usr/bin/env python3
"""Dry-run audit for the THE-32 upstream low-calo event-quality filter.

This script intentionally reads only event-quality columns from
AuAuPhotonIDTrainingTree ROOT files.  It does not train, score, inspect truth
labels, inspect candidate shower variables, or use weights/splits/WP behavior.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


TREE_NAME = "AuAuPhotonIDTrainingTree"
GLOBAL_EVENT_KEY_COLUMNS = ["source_sample", "input_file", "run", "evt"]
DIRECT_CUT_COLUMNS = ["centrality", "event_calo_log10_total_energy_plus1"]
COMPONENT_COLUMNS = [
    "event_calo_cemc_energy",
    "event_calo_ihcal_energy",
    "event_calo_ohcal_energy",
    "event_calo_total_energy",
]
CENTRALITY_BINS = [
    ("0-20", 0.0, 20.0),
    ("20-50", 20.0, 50.0),
    ("50-80", 50.0, 80.0),
]
SAMPLE_ALIASES = (
    ("run28_embeddedPhoton12", ("run28_embeddedPhoton12", "embeddedPhoton12", "Photon12")),
    ("run28_embeddedPhoton20", ("run28_embeddedPhoton20", "embeddedPhoton20", "Photon20")),
    ("run28_embeddedJet12", ("run28_embeddedJet12", "embeddedJet12", "Jet12")),
    ("run28_embeddedJet20", ("run28_embeddedJet20", "embeddedJet20", "Jet20")),
    ("run28_embeddedJet30", ("run28_embeddedJet30", "embeddedJet30", "Jet30")),
    ("run28_embeddedJet40", ("run28_embeddedJet40", "embeddedJet40", "Jet40")),
)
SAMPLE_TO_CODE = {sample: idx + 1 for idx, (sample, _) in enumerate(SAMPLE_ALIASES)}
CODE_TO_SAMPLE = {idx + 1: sample for idx, (sample, _) in enumerate(SAMPLE_ALIASES)}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_ready(obj):
    if isinstance(obj, dict):
        return {key: json_ready(value) for key, value in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_ready(value) for value in obj]
    if isinstance(obj, np.ndarray):
        return json_ready(obj.tolist())
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    return obj


def git_status() -> str:
    try:
        return subprocess.check_output(["git", "status", "--short"], text=True, stderr=subprocess.STDOUT)
    except Exception as exc:
        return f"git status unavailable: {exc}"


def expand_manifest(path: Path) -> list[Path]:
    paths = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line and not line.startswith("#"):
            paths.append(Path(line))
    if not paths:
        raise SystemExit(f"No ROOT files listed in {path}")
    return paths


def infer_source_sample(path: Path) -> str:
    text = str(path)
    for sample, aliases in SAMPLE_ALIASES:
        if any(alias in text for alias in aliases):
            return sample
    return "unknown"


def load_cut(path: Path) -> dict:
    if not path.is_file():
        raise SystemExit(f"Cut JSON does not exist: {path}")
    payload = json.loads(path.read_text())
    envelope = payload.get("envelope")
    if not isinstance(envelope, list) or not envelope:
        raise SystemExit(f"Cut JSON has no envelope table: {path}")
    rows = []
    for item in envelope:
        rows.append(
            {
                "cent_lo": float(item["cent_lo"]),
                "cent_hi": float(item["cent_hi"]),
                "threshold": float(item["threshold"]),
                "median": float(item.get("median", math.nan)),
                "mad_sigma": float(item.get("mad_sigma", math.nan)),
                "quantile_floor": float(item.get("quantile_floor", math.nan)),
                "n_events": int(item.get("n_events", 0)),
                "status": str(item.get("status", "unknown")),
            }
        )
    rows.sort(key=lambda row: (row["cent_lo"], row["cent_hi"]))
    return {
        "path": str(path),
        "sha256": file_sha256(path),
        "schema": payload.get("schema"),
        "source_blind": bool(payload.get("source_blind", False)),
        "truth_blind": bool(payload.get("truth_blind", False)),
        "bdt_score_blind": bool(payload.get("bdt_score_blind", False)),
        "centrality_source": payload.get("centrality_source"),
        "cut_variable": payload.get("cut_variable"),
        "mad_scale": payload.get("mad_scale"),
        "quantile_floor": payload.get("quantile_floor"),
        "slice_width_percent": payload.get("slice_width_percent"),
        "envelope": rows,
    }


def thresholds_for_cent(cent: np.ndarray, envelope: list[dict]) -> np.ndarray:
    out = np.full(len(cent), np.nan, dtype="float64")
    for row in envelope:
        mask = (cent >= float(row["cent_lo"])) & (cent < float(row["cent_hi"]))
        out[mask] = float(row["threshold"])
    return out


def update_fraction_counts(store: dict, key: str, total_candidates: int, rejected_candidates: int, total_events: int, rejected_events: int) -> None:
    item = store[key]
    item["total_candidates"] += int(total_candidates)
    item["rejected_candidates"] += int(rejected_candidates)
    item["total_events"] += int(total_events)
    item["rejected_events"] += int(rejected_events)


def fraction_rows(store: dict, label_name: str) -> list[dict]:
    rows = []
    for key in sorted(store):
        item = store[key]
        total_events = item["total_events"]
        total_candidates = item["total_candidates"]
        rejected_events = item["rejected_events"]
        rejected_candidates = item["rejected_candidates"]
        event_frac = rejected_events / total_events if total_events else math.nan
        cand_frac = rejected_candidates / total_candidates if total_candidates else math.nan
        rows.append(
            {
                label_name: key,
                "total_events": total_events,
                "rejected_events": rejected_events,
                "event_rejection_fraction": event_frac,
                "total_candidates": total_candidates,
                "rejected_candidates": rejected_candidates,
                "candidate_rejection_fraction": cand_frac,
                "candidate_minus_event_fraction": (
                    cand_frac - event_frac if math.isfinite(event_frac) and math.isfinite(cand_frac) else math.nan
                ),
            }
        )
    return rows


def centrality_label(cent: float) -> str | None:
    for label, lo, hi in CENTRALITY_BINS:
        if cent >= lo and cent < hi:
            return label
    return None


def safe_median(values: list[float]) -> float:
    if not values:
        return math.nan
    arr = np.asarray(values, dtype="float64")
    arr = arr[np.isfinite(arr)]
    return float(np.median(arr)) if len(arr) else math.nan


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--cut-json", type=Path, required=True)
    ap.add_argument("--output-json", type=Path, required=True)
    ap.add_argument("--event-table-npz", type=Path, required=True)
    ap.add_argument("--rejection-csv", type=Path, required=True)
    ap.add_argument("--component-csv", type=Path, required=True)
    ap.add_argument("--tree-name", default=TREE_NAME)
    ap.add_argument("--progress-every", type=int, default=250)
    args = ap.parse_args()

    try:
        import pandas as pd
        import uproot
    except ImportError as exc:
        raise SystemExit("This audit requires uproot and pandas in the ML environment.") from exc

    root_paths = expand_manifest(args.manifest)
    cut = load_cut(args.cut_json)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.event_table_npz.parent.mkdir(parents=True, exist_ok=True)

    cent_by_event: list[np.ndarray] = []
    log_by_event: list[np.ndarray] = []
    threshold_by_event: list[np.ndarray] = []
    cemc_by_event: list[np.ndarray] = []
    ihcal_by_event: list[np.ndarray] = []
    ohcal_by_event: list[np.ndarray] = []
    total_by_event: list[np.ndarray] = []
    source_code_by_event: list[np.ndarray] = []
    rejected_by_event: list[np.ndarray] = []
    candidate_count_by_event: list[np.ndarray] = []
    file_index_by_event: list[np.ndarray] = []
    run_by_event: list[np.ndarray] = []
    evt_by_event: list[np.ndarray] = []

    source_counts = defaultdict(lambda: {"total_events": 0, "rejected_events": 0, "total_candidates": 0, "rejected_candidates": 0})
    cent_counts = defaultdict(lambda: {"total_events": 0, "rejected_events": 0, "total_candidates": 0, "rejected_candidates": 0})
    comp_values: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    full_event_keys: set[str] = set()
    simple_event_keys: set[str] = set()
    duplicate_full_event_keys = 0
    candidate_rows = 0
    rejected_candidate_rows = 0
    skipped_missing_tree: list[str] = []
    skipped_missing_branches: list[dict] = []
    branches_read_seen: set[str] = set()

    requested = [
        "run",
        "evt",
        "centrality",
        "event_calo_log10_total_energy_plus1",
        "event_calo_cemc_energy",
        "event_calo_ihcal_energy",
        "event_calo_ohcal_energy",
        "event_calo_total_energy",
        "source_sample",
    ]

    for file_index, path in enumerate(root_paths):
        if args.progress_every > 0 and (file_index == 0 or (file_index + 1) % args.progress_every == 0 or file_index + 1 == len(root_paths)):
            print(f"[THE32 upstream audit] reading {file_index + 1}/{len(root_paths)} {path}", flush=True)
        source_fallback = infer_source_sample(path)
        try:
            with uproot.open(path) as root_file:
                if args.tree_name not in root_file:
                    skipped_missing_tree.append(str(path))
                    continue
                tree = root_file[args.tree_name]
                keys = set(tree.keys())
                required = {"run", "evt", "centrality", "event_calo_log10_total_energy_plus1"}
                missing = sorted(required.difference(keys))
                if missing:
                    skipped_missing_branches.append({"path": str(path), "missing": missing})
                    continue
                present = [name for name in requested if name in keys]
                branches_read_seen.update(present)
                frame = tree.arrays(present, library="pd")
        except Exception as exc:
            skipped_missing_branches.append({"path": str(path), "error": str(exc)})
            continue

        if len(frame) == 0:
            continue
        if "source_sample" in frame.columns:
            source_arr = frame["source_sample"].astype(str).to_numpy(dtype=object)
        else:
            source_arr = np.full(len(frame), source_fallback, dtype=object)
        run = frame["run"].to_numpy(dtype="int64", copy=False)
        evt = frame["evt"].to_numpy(dtype="int64", copy=False)
        cent = frame["centrality"].to_numpy(dtype="float64", copy=False)
        log_calo = frame["event_calo_log10_total_energy_plus1"].to_numpy(dtype="float64", copy=False)
        thresholds = thresholds_for_cent(cent, cut["envelope"])
        in_range = np.isfinite(cent) & np.isfinite(log_calo) & np.isfinite(thresholds)
        row_below = in_range & (log_calo < thresholds)

        key_frame = pd.DataFrame(
            {
                "source_sample": source_arr,
                "run": run,
                "evt": evt,
                "centrality": cent,
                "log_calo": log_calo,
                "threshold": thresholds,
                "row_below": row_below,
                "candidate_count": np.ones(len(frame), dtype="int64"),
            }
        )
        for col in ["event_calo_cemc_energy", "event_calo_ihcal_energy", "event_calo_ohcal_energy"]:
            key_frame[col] = frame[col].to_numpy(dtype="float64", copy=False) if col in frame.columns else np.full(len(frame), np.nan)
        if "event_calo_total_energy" in frame.columns:
            key_frame["event_calo_total_energy"] = frame["event_calo_total_energy"].to_numpy(dtype="float64", copy=False)
        else:
            key_frame["event_calo_total_energy"] = (
                key_frame["event_calo_cemc_energy"].to_numpy(dtype="float64")
                + key_frame["event_calo_ihcal_energy"].to_numpy(dtype="float64")
                + key_frame["event_calo_ohcal_energy"].to_numpy(dtype="float64")
            )

        grouped = key_frame.groupby(["source_sample", "run", "evt"], sort=False, dropna=False)
        events = grouped.agg(
            centrality=("centrality", "first"),
            log_calo=("log_calo", "first"),
            threshold=("threshold", "first"),
            reject=("row_below", "any"),
            candidate_count=("candidate_count", "sum"),
            event_calo_cemc_energy=("event_calo_cemc_energy", "first"),
            event_calo_ihcal_energy=("event_calo_ihcal_energy", "first"),
            event_calo_ohcal_energy=("event_calo_ohcal_energy", "first"),
            event_calo_total_energy=("event_calo_total_energy", "first"),
        ).reset_index()
        events["source_sample"] = events["source_sample"].astype(str)
        source_codes = np.asarray([SAMPLE_TO_CODE.get(src, 0) for src in events["source_sample"].to_numpy(dtype=object)], dtype="int16")
        event_reject = events["reject"].to_numpy(dtype=bool)
        event_candidates = events["candidate_count"].to_numpy(dtype="int64")
        candidate_rows += int(event_candidates.sum())
        rejected_candidate_rows += int(event_candidates[event_reject].sum())

        for row in events.itertuples(index=False):
            full_key = f"{row.source_sample}|{path}|{int(row.run)}|{int(row.evt)}"
            simple_key = f"{row.source_sample}|{int(row.run)}|{int(row.evt)}"
            if full_key in full_event_keys:
                duplicate_full_event_keys += 1
            full_event_keys.add(full_key)
            simple_event_keys.add(simple_key)
            cent_label = centrality_label(float(row.centrality))
            rejected = bool(row.reject)
            sample_label = str(row.source_sample)
            update_fraction_counts(source_counts, sample_label, int(row.candidate_count), int(row.candidate_count) if rejected else 0, 1, 1 if rejected else 0)
            if cent_label is not None:
                update_fraction_counts(cent_counts, cent_label, int(row.candidate_count), int(row.candidate_count) if rejected else 0, 1, 1 if rejected else 0)
                for component, value in [
                    ("CEMC", row.event_calo_cemc_energy),
                    ("IHCal", row.event_calo_ihcal_energy),
                    ("OHCal", row.event_calo_ohcal_energy),
                    ("total_calo", row.event_calo_total_energy),
                    ("log10_total_plus1", row.log_calo),
                ]:
                    comp_values[(cent_label, component, "rejected" if rejected else "retained")].append(float(value))

        cent_by_event.append(events["centrality"].to_numpy(dtype="float32"))
        log_by_event.append(events["log_calo"].to_numpy(dtype="float32"))
        threshold_by_event.append(events["threshold"].to_numpy(dtype="float32"))
        cemc_by_event.append(events["event_calo_cemc_energy"].to_numpy(dtype="float32"))
        ihcal_by_event.append(events["event_calo_ihcal_energy"].to_numpy(dtype="float32"))
        ohcal_by_event.append(events["event_calo_ohcal_energy"].to_numpy(dtype="float32"))
        total_by_event.append(events["event_calo_total_energy"].to_numpy(dtype="float32"))
        source_code_by_event.append(source_codes)
        rejected_by_event.append(event_reject)
        candidate_count_by_event.append(event_candidates.astype("int64"))
        file_index_by_event.append(np.full(len(events), file_index, dtype="int32"))
        run_by_event.append(events["run"].to_numpy(dtype="int64"))
        evt_by_event.append(events["evt"].to_numpy(dtype="int64"))

    if not cent_by_event:
        raise SystemExit("No usable event rows found for low-calo audit.")

    event_cent = np.concatenate(cent_by_event)
    event_log = np.concatenate(log_by_event)
    event_threshold = np.concatenate(threshold_by_event)
    event_reject = np.concatenate(rejected_by_event).astype(bool)
    event_candidate_count = np.concatenate(candidate_count_by_event).astype("int64")
    retained_below = (~event_reject) & np.isfinite(event_log) & np.isfinite(event_threshold) & (event_log < event_threshold)
    unique_events = int(len(event_cent))
    rejected_events = int(event_reject.sum())

    component_rows = []
    for cent_label, _, _ in CENTRALITY_BINS:
        for component in ["CEMC", "IHCal", "OHCal", "total_calo", "log10_total_plus1"]:
            retained_median = safe_median(comp_values[(cent_label, component, "retained")])
            rejected_median = safe_median(comp_values[(cent_label, component, "rejected")])
            component_rows.append(
                {
                    "centrality_bin": cent_label,
                    "component": component,
                    "retained_median": retained_median,
                    "rejected_median": rejected_median,
                    "rejected_over_retained": (
                        rejected_median / retained_median
                        if math.isfinite(retained_median) and retained_median != 0.0 and math.isfinite(rejected_median)
                        else math.nan
                    ),
                }
            )

    rejection_rows = []
    for row in fraction_rows(cent_counts, "centrality_bin"):
        row["axis"] = "centrality"
        rejection_rows.append(row)
    for row in fraction_rows(source_counts, "source_sample"):
        row["axis"] = "source_sample"
        rejection_rows.append(row)

    write_csv(args.rejection_csv, rejection_rows)
    write_csv(args.component_csv, component_rows)
    np.savez(
        args.event_table_npz,
        centrality=event_cent,
        log_calo=event_log,
        threshold=event_threshold,
        cemc=np.concatenate(cemc_by_event).astype("float32"),
        ihcal=np.concatenate(ihcal_by_event).astype("float32"),
        ohcal=np.concatenate(ohcal_by_event).astype("float32"),
        total_calo=np.concatenate(total_by_event).astype("float32"),
        source_code=np.concatenate(source_code_by_event).astype("int16"),
        rejected=event_reject,
        candidate_count=event_candidate_count,
        file_index=np.concatenate(file_index_by_event).astype("int32"),
        run=np.concatenate(run_by_event).astype("int64"),
        evt=np.concatenate(evt_by_event).astype("int64"),
        source_code_samples=np.asarray([CODE_TO_SAMPLE[i] for i in sorted(CODE_TO_SAMPLE)], dtype=object),
    )

    audit = {
        "schema": "THE32_UPSTREAM_LOW_CALO_FAST_DRYRUN_AUDIT_V1",
        "manifest": str(args.manifest),
        "input_file_count": len(root_paths),
        "tree_name": args.tree_name,
        "cut_json_path": str(args.cut_json),
        "cut_json_sha256": cut["sha256"],
        "cut_json_schema": cut["schema"],
        "source_blind": cut["source_blind"],
        "truth_blind": cut["truth_blind"],
        "bdt_score_blind": cut["bdt_score_blind"],
        "variables_used_for_threshold_derivation_and_application": list(DIRECT_CUT_COLUMNS),
        "variables_used_only_for_event_identity": list(GLOBAL_EVENT_KEY_COLUMNS),
        "variables_used_only_for_outcome_audit": ["source_sample", *COMPONENT_COLUMNS],
        "variables_explicitly_not_read": [
            "is_signal",
            "truth photon label",
            "truth isolation label",
            "BDT score",
            "cluster_Et",
            "cluster_Eta",
            "candidate shower-shape variables",
            "train/test split",
            "sample weight",
            "WP80 behavior",
        ],
        "branches_read_seen": sorted(branches_read_seen),
        "boundary_convention": "centrality slice uses cent_lo <= centrality < cent_hi; event is rejected only when log10(total calo + 1) < threshold; equality is retained",
        "threshold_table": cut["envelope"],
        "candidate_rows_before": int(candidate_rows),
        "candidate_rows_after": int(candidate_rows - rejected_candidate_rows),
        "candidate_rows_rejected": int(rejected_candidate_rows),
        "globally_unique_events_before": unique_events,
        "globally_unique_events_after": int(unique_events - rejected_events),
        "globally_unique_events_rejected": rejected_events,
        "retained_below_envelope_events": int(retained_below.sum()),
        "retained_below_envelope_candidates": int(event_candidate_count[retained_below].sum()),
        "global_event_key_audit": {
            "present": True,
            "event_key_columns": list(GLOBAL_EVENT_KEY_COLUMNS),
            "candidate_key_columns": [*GLOBAL_EVENT_KEY_COLUMNS, "input_tree_entry"],
            "event_dedup_relies_only_on_source_run_evt": False,
            "globally_unique_events": unique_events,
            "duplicate_global_event_keys": duplicate_full_event_keys,
            "source_run_evt_unique_keys": int(len(simple_event_keys)),
            "source_run_evt_repeated_after_global_event_dedup": int(unique_events - len(simple_event_keys)),
            "candidate_rows": int(candidate_rows),
            "candidate_multiplicity_mean": float(candidate_rows / unique_events) if unique_events else math.nan,
        },
        "rejection_by_centrality": [row for row in rejection_rows if row["axis"] == "centrality"],
        "rejection_by_source": [row for row in rejection_rows if row["axis"] == "source_sample"],
        "component_medians_by_centrality": component_rows,
        "outputs": {
            "audit_json": str(args.output_json),
            "event_table_npz": str(args.event_table_npz),
            "rejection_csv": str(args.rejection_csv),
            "component_csv": str(args.component_csv),
        },
        "skipped_missing_tree": skipped_missing_tree,
        "skipped_missing_branches_or_errors": skipped_missing_branches,
        "command": " ".join(sys.argv),
        "git_status_short": git_status(),
    }
    if int(retained_below.sum()) != 0:
        bad_idx = np.flatnonzero(retained_below)[:10]
        audit["retained_below_examples"] = [
            {
                "centrality": float(event_cent[i]),
                "log_calo": float(event_log[i]),
                "threshold": float(event_threshold[i]),
            }
            for i in bad_idx
        ]

    args.output_json.write_text(json.dumps(json_ready(audit), indent=2, sort_keys=True) + "\n")
    print(f"[OK] wrote dry-run audit JSON: {args.output_json}", flush=True)
    print(f"[OK] wrote event table NPZ: {args.event_table_npz}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Canonical numeric gate for THE-76 PPG12 photon+jet Fig.5 stitching.

This checker intentionally validates the strict source-stage Fig.5 contract:
one photon5/photon10/photon20 source population, pre-firstEventCuts/vz fill,
PPG12 slimTree denominators, PPG12 sample ownership windows, and no global
scale. It does not validate downstream final-yield period/SI-DI products.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import statistics
from typing import Any


REPO = Path(__file__).resolve().parents[3]
DEFAULT_POINTS = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig5_prevz_strict_20260705_224243/"
    "strict_stitched_photonjet_prevz/"
    "photon_data_over_fit_sdcc_vs_current_overlay_prevz_strict_points.csv"
)
DEFAULT_MANIFEST = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig5_prevz_strict_20260705_224243/"
    "strict_stitched_photonjet_prevz/"
    "photon_data_over_fit_sdcc_vs_current_overlay_prevz_strict_manifest.json"
)
DEFAULT_REPORT = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig5_prevz_strict_20260705_224243/"
    "strict_stitched_photonjet_prevz/photon_fig5_canonical_check.json"
)

EXPECTED_DENOMINATORS = {
    "photon5": 9986593,
    "photon10": 9998561,
    "photon20": 9999988,
}
EXPECTED_BIN_COUNTS = {
    "photon5": 8,
    "photon10": 16,
    "photon20": 36,
}
EXPECTED_WINDOWS = {
    "photon5": (10.0, 14.0),
    "photon10": (14.0, 22.0),
    "photon20": (22.0, 40.0),
}


def as_float(row: dict[str, str], key: str) -> float:
    try:
        return float(row[key])
    except Exception as exc:  # pragma: no cover - defensive error context
        raise ValueError(f"cannot parse {key}={row.get(key)!r}") from exc


def as_int(row: dict[str, str], key: str) -> int:
    try:
        return int(float(row[key]))
    except Exception as exc:  # pragma: no cover - defensive error context
        raise ValueError(f"cannot parse {key}={row.get(key)!r}") from exc


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=Path, default=DEFAULT_POINTS)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--json-output", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--max-abs-current-over-sdcc-minus-one", type=float, default=1.0e-4)
    parser.add_argument("--max-rms-current-over-sdcc-minus-one", type=float, default=2.0e-5)
    parser.add_argument("--expected-bins", type=int, default=60)
    return parser


def check_manifest(path: Path) -> tuple[list[str], dict[str, Any]]:
    failures: list[str] = []
    if not path.exists():
        return [f"missing manifest: {path}"], {}
    payload = json.loads(path.read_text())
    if payload.get("status") != "ok":
        failures.append(f"manifest status is not ok: {payload.get('status')!r}")
    if payload.get("campaign_tag") != "the76_ppg12_fig5_prevz_strict_20260705_224243":
        failures.append(f"unexpected campaign_tag: {payload.get('campaign_tag')!r}")
    if payload.get("histogram") != "SIM/h_ppPhotonStitch_ppg12Fig5PreVz_maxPhotonPt_kept":
        failures.append(f"unexpected histogram: {payload.get('histogram')!r}")

    samples = payload.get("samples") or {}
    for sample in EXPECTED_DENOMINATORS:
        sample_payload = samples.get(sample) or {}
        if sample_payload.get("files") != 1000:
            failures.append(f"{sample}: expected 1000 source files, found {sample_payload.get('files')!r}")
        if sample_payload.get("missing_histogram") != 0:
            failures.append(f"{sample}: missing_histogram={sample_payload.get('missing_histogram')!r}")
        if sample_payload.get("zombie") != 0:
            failures.append(f"{sample}: zombie={sample_payload.get('zombie')!r}")
    return failures, payload


def check_rows(rows: list[dict[str, str]], args: argparse.Namespace) -> tuple[list[str], dict[str, Any]]:
    failures: list[str] = []
    if len(rows) != args.expected_bins:
        failures.append(f"expected {args.expected_bins} rows, found {len(rows)}")

    ratios = [as_float(row, "current_over_sdcc") for row in rows]
    max_abs = max((abs(ratio - 1.0) for ratio in ratios), default=math.inf)
    rms = math.sqrt(sum((ratio - 1.0) ** 2 for ratio in ratios) / len(ratios)) if ratios else math.inf
    if max_abs > args.max_abs_current_over_sdcc_minus_one:
        failures.append(
            "max |current_over_sdcc - 1| "
            f"{max_abs:.12g} exceeds {args.max_abs_current_over_sdcc_minus_one:.12g}"
        )
    if rms > args.max_rms_current_over_sdcc_minus_one:
        failures.append(
            "RMS(current_over_sdcc - 1) "
            f"{rms:.12g} exceeds {args.max_rms_current_over_sdcc_minus_one:.12g}"
        )

    rows_by_sample: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        rows_by_sample.setdefault(row["sample"], []).append(row)

    for sample, expected_count in EXPECTED_BIN_COUNTS.items():
        sample_rows = rows_by_sample.get(sample, [])
        if len(sample_rows) != expected_count:
            failures.append(f"{sample}: expected {expected_count} stitched bins, found {len(sample_rows)}")
            continue
        low_expected, high_expected = EXPECTED_WINDOWS[sample]
        lows = [as_float(row, "bin_low") for row in sample_rows]
        highs = [as_float(row, "bin_high") for row in sample_rows]
        if min(lows) < low_expected or max(highs) > high_expected:
            failures.append(f"{sample}: bin window {min(lows)}-{max(highs)} outside {low_expected}-{high_expected}")
        for row in sample_rows:
            if row.get("current_value_mode") != "prevz_raw_count_times_xsec_over_ppg12_slimtree_entries":
                failures.append(f"{sample}: unexpected current_value_mode={row.get('current_value_mode')!r}")
            if "strict pre-vz" not in row.get("comparison_scope", ""):
                failures.append(f"{sample}: comparison_scope does not state strict pre-vz")
            if "no period/SI-DI aggregate" not in row.get("comparison_scope", ""):
                failures.append(f"{sample}: comparison_scope does not exclude period/SI-DI aggregate")
            if "no global scale" not in row.get("comparison_scope", ""):
                failures.append(f"{sample}: comparison_scope does not exclude global scale")
            if as_int(row, "current_source_files") != 1000:
                failures.append(f"{sample}: current_source_files={row.get('current_source_files')!r}")
            if as_int(row, "current_missing_histogram_files") != 0:
                failures.append(f"{sample}: current_missing_histogram_files={row.get('current_missing_histogram_files')!r}")
            if as_int(row, "current_zombie_files") != 0:
                failures.append(f"{sample}: current_zombie_files={row.get('current_zombie_files')!r}")
            if as_int(row, "ppg12_slimtree_entries") != EXPECTED_DENOMINATORS[sample]:
                failures.append(
                    f"{sample}: ppg12_slimtree_entries={row.get('ppg12_slimtree_entries')!r}, "
                    f"expected {EXPECTED_DENOMINATORS[sample]}"
                )

    per_sample = {}
    for sample, sample_rows in sorted(rows_by_sample.items()):
        sample_ratios = [as_float(row, "current_over_sdcc") for row in sample_rows]
        per_sample[sample] = {
            "bins": len(sample_ratios),
            "min_current_over_sdcc": min(sample_ratios),
            "max_current_over_sdcc": max(sample_ratios),
            "max_abs_current_over_sdcc_minus_one": max(abs(ratio - 1.0) for ratio in sample_ratios),
            "rms_current_over_sdcc_minus_one": math.sqrt(
                sum((ratio - 1.0) ** 2 for ratio in sample_ratios) / len(sample_ratios)
            ),
        }

    summary = {
        "rows": len(rows),
        "min_current_over_sdcc": min(ratios) if ratios else None,
        "max_current_over_sdcc": max(ratios) if ratios else None,
        "mean_current_over_sdcc": sum(ratios) / len(ratios) if ratios else None,
        "median_current_over_sdcc": statistics.median(ratios) if ratios else None,
        "max_abs_current_over_sdcc_minus_one": max_abs,
        "rms_current_over_sdcc_minus_one": rms,
        "per_sample": per_sample,
    }
    return failures, summary


def main() -> int:
    args = build_parser().parse_args()
    args.points = args.points if args.points.is_absolute() else REPO / args.points
    args.manifest = args.manifest if args.manifest.is_absolute() else REPO / args.manifest
    args.json_output = args.json_output if args.json_output.is_absolute() else REPO / args.json_output

    failures: list[str] = []
    manifest_payload: dict[str, Any] = {}
    if not args.points.exists():
        failures.append(f"missing points CSV: {args.points}")
        rows: list[dict[str, str]] = []
    else:
        rows = load_rows(args.points)

    manifest_failures, manifest_payload = check_manifest(args.manifest)
    failures.extend(manifest_failures)
    row_failures, row_summary = check_rows(rows, args)
    failures.extend(row_failures)

    report = {
        "schema_version": 1,
        "contract_id": "pp_photonjet_stitch_fig5",
        "status": "pass" if not failures else "fail",
        "points_csv": str(args.points),
        "manifest": str(args.manifest),
        "manifest_campaign_tag": manifest_payload.get("campaign_tag"),
        "manifest_histogram": manifest_payload.get("histogram"),
        "reference": {
            "ppg12_remote_root": "/sphenix/user/shuhangli/ppg12/plotting/photon_max_pT_uncut.root",
            "ppg12_objects": [
                "h_max_photon_pT_photon5",
                "h_max_photon_pT_photon10",
                "h_max_photon_pT_photon20",
            ],
            "ppg12_slimtree_entries": EXPECTED_DENOMINATORS,
        },
        "required_contract": {
            "fill_stage": "pre-firstEventCuts/vz source-stage diagnostic",
            "histogram": "SIM/h_ppPhotonStitch_ppg12Fig5PreVz_maxPhotonPt_kept",
            "sample_windows": {
                "photon5": "10 <= pT < 14 GeV",
                "photon10": "14 <= pT < 22 GeV",
                "photon20": "22 <= pT < 40 GeV",
            },
            "normalization": "XSEC[sample] / PPG12 slimtree.num_entries",
            "excluded": [
                "period/SI-DI aggregate",
                "vertex/luminosity/mix/reweight machinery",
                "global fitted scale",
            ],
        },
        "thresholds": {
            "max_abs_current_over_sdcc_minus_one": args.max_abs_current_over_sdcc_minus_one,
            "max_rms_current_over_sdcc_minus_one": args.max_rms_current_over_sdcc_minus_one,
            "expected_bins": args.expected_bins,
        },
        "summary": row_summary,
        "failures": failures,
    }
    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())

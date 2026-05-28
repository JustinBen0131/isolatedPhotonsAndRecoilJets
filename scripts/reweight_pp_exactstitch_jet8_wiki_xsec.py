#!/usr/bin/env python3
"""Recompute pp exact-stitch compact QA with the wiki jet8 cross section.

This is intentionally a downstream normalization variation.  It reuses the
in-situ RecoilJets raw histogram counts already pulled locally and changes only
the run28_jet8 xsec from the historical local value to the wiki value.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


JET8_SAMPLE = "run28_jet8"
JET8_WIKI_XSEC_PB = 1.3013e7
JET_SAMPLES = ["run28_jet8", "run28_jet12", "run28_jet20", "run28_jet30", "run28_jet40"]
JET_BOUNDARIES = [
    ("run28_jet8", "run28_jet12", 14.0),
    ("run28_jet12", "run28_jet20", 21.0),
    ("run28_jet20", "run28_jet30", 32.0),
    ("run28_jet30", "run28_jet40", 42.0),
]


def fnum(x: float) -> str:
    return f"{x:.12g}"


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f))


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def update_contract_points(in_csv: Path, out_csv: Path) -> None:
    rows = load_rows(in_csv)
    for row in rows:
        if row["sample"] != JET8_SAMPLE:
            continue
        old_xsec = float(row["xsec_pb"])
        if old_xsec <= 0:
            raise RuntimeError("jet8 old xsec is non-positive")
        scale = JET8_WIKI_XSEC_PB / old_xsec
        row["xsec_pb"] = fnum(JET8_WIKI_XSEC_PB)
        row["density_pb_per_gev"] = fnum(float(row["density_pb_per_gev"]) * scale)
        row["density_err_pb_per_gev"] = fnum(float(row["density_err_pb_per_gev"]) * scale)
    write_rows(out_csv, rows)


def update_contract_summary(in_json: Path, out_json: Path) -> None:
    data = json.loads(in_json.read_text())
    for sample in data["samples"]:
        if sample["sample"] != JET8_SAMPLE:
            continue
        sample["xsec_pb"] = JET8_WIKI_XSEC_PB
        sample["metadata"]["xsec_pb"] = JET8_WIKI_XSEC_PB
        sample["xsec_note"] = "Updated from sPHENIX Jet Structure Topical Group wiki Pythia simulations table."
    data["normalization_variant"] = "run28_jet8 xsec set to wiki value 1.3013e7 pb; raw in-situ counts unchanged"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(data, indent=2) + "\n")


def rows_by_sample(rows: list[dict[str, str]]) -> dict[str, dict[float, dict[str, str]]]:
    out: dict[str, dict[float, dict[str, str]]] = {}
    for row in rows:
        out.setdefault(row["sample"], {})[float(row["bin_center"])] = row
    return out


def all_density(row: dict[str, str]) -> float:
    return (
        float(row["raw_all_events"])
        * float(row["xsec_pb"])
        / float(row["events_processed_metadata"])
        / (float(row["bin_high"]) - float(row["bin_low"]))
    )


def recompute_boundaries(points_csv: Path, out_csv: Path, out_json: Path) -> None:
    rows = [row for row in load_rows(points_csv) if row["group"] == "jet"]
    by_sample = rows_by_sample(rows)
    out_rows: list[dict[str, object]] = []
    for left, right, boundary in JET_BOUNDARIES:
        left_bins = by_sample[left]
        right_bins = by_sample[right]
        bin_width = float(next(iter(left_bins.values()))["bin_high"]) - float(next(iter(left_bins.values()))["bin_low"])
        left_center = boundary - 0.5 * bin_width
        right_center = boundary + 0.5 * bin_width
        l_at_left = left_bins[left_center]
        r_at_left = right_bins[left_center]
        l_at_right = left_bins[right_center]
        r_at_right = right_bins[right_center]
        left_display = float(l_at_left["density_pb_per_gev"])
        right_display = float(r_at_right["density_pb_per_gev"])
        left_center_l = all_density(l_at_left)
        left_center_r = all_density(r_at_left)
        right_center_l = all_density(l_at_right)
        right_center_r = all_density(r_at_right)
        overlap_left = left_center_l / left_center_r if left_center_r > 0 else math.nan
        overlap_right = right_center_l / right_center_r if right_center_r > 0 else math.nan
        devs = [abs(overlap_left - 1.0), abs(overlap_right - 1.0)]
        max_dev = max(d for d in devs if math.isfinite(d))
        display_ratio = left_display / right_display if right_display > 0 else math.nan
        out_rows.append(
            {
                "group": "jet",
                "boundary_GeV": boundary,
                "left_sample": left,
                "right_sample": right,
                "window_left_GeV": f"[{float(l_at_left['stitch_window_low'])},{float(l_at_left['stitch_window_high'])})",
                "window_right_GeV": f"[{float(r_at_right['stitch_window_low'])},{float(r_at_right['stitch_window_high'])})",
                "bin_width_GeV": bin_width,
                "left_edge_center_GeV": left_center,
                "right_edge_center_GeV": right_center,
                "left_display_density_pb_per_GeV": left_display,
                "right_display_density_pb_per_GeV": right_display,
                "display_adjacent_left_over_right": display_ratio,
                "left_center_all_density_left_sample": left_center_l,
                "left_center_all_density_right_sample": left_center_r,
                "overlap_left_center_left_over_right": overlap_left,
                "right_center_all_density_left_sample": right_center_l,
                "right_center_all_density_right_sample": right_center_r,
                "overlap_right_center_left_over_right": overlap_right,
                "max_same_bin_overlap_fractional_deviation": max_dev,
                "left_xsec_pb": float(l_at_left["xsec_pb"]),
                "right_xsec_pb": float(r_at_right["xsec_pb"]),
                "left_events_processed": float(l_at_left["events_processed_metadata"]),
                "right_events_processed": float(r_at_right["events_processed_metadata"]),
                "no_window_gap_or_overlap": True,
                "pass_same_bin_overlap_within_15pct": max_dev <= 0.15,
                "pass_visible_falling_spectrum_continuity": display_ratio > 1.0,
                "pass_boundary": max_dev <= 0.15 and display_ratio > 1.0,
            }
        )

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)
    out_json.write_text(
        json.dumps(
            {
                "schema": "PP_CURRENTIAN_EXACTSTITCH_JET8_WIKI_XSEC_BOUNDARY_QA_V1",
                "normalization_variant": "run28_jet8 xsec set to 1.3013e7 pb",
                "boundaries": out_rows,
            },
            indent=2,
        )
        + "\n"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qa-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    in_points = args.qa_dir / "pp_currentian_exactstitch_contract_points.csv"
    in_summary = args.qa_dir / "pp_currentian_exactstitch_contract_summary.json"
    out_points = args.out_dir / "pp_currentian_exactstitch_contract_points_jet8WikiXsec.csv"
    out_summary = args.out_dir / "pp_currentian_exactstitch_contract_summary_jet8WikiXsec.json"
    out_boundary_csv = args.out_dir / "pp_currentian_exactstitch_boundary_continuity_qa_jet8WikiXsec.csv"
    out_boundary_json = args.out_dir / "pp_currentian_exactstitch_boundary_continuity_qa_jet8WikiXsec.json"

    update_contract_points(in_points, out_points)
    update_contract_summary(in_summary, out_summary)
    recompute_boundaries(out_points, out_boundary_csv, out_boundary_json)
    print(out_points)
    print(out_summary)
    print(out_boundary_csv)
    print(out_boundary_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

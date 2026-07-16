#!/usr/bin/env python3
"""Build THE-94 inclusive-jet Fig.6 canonical jet8 policy plot points.

This helper is deliberately plotting-only.  It starts from the audited
no-scale THE-94 vs true-period-combined PPG12 Fig.6 table and preserves the
THE-94 nominal values exactly.  In particular, jet8 is not multiplied by the
historical PPG12 hidden-source factor; the point of this artifact is to show
the smooth RecoilJets canonical jet8->jet12 stitch under the visible
``xsec/jet50 * lumi * mix * vertex-weight`` contract.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from statistics import median


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
THE94_BASE = (
    REPO
    / "dataOutput/ppg12Parity/the94_ppg12_inclusivejet_fig6_fixed_20260702_204032"
)
IN_POINTS = (
    THE94_BASE
    / "no_scale_overlay/inclusivejet_fig6_true_period_combined_no_scale_data_over_fit_points.csv"
)
OUT_DIR = THE94_BASE / "canonical_jet8_policy_overlay"
OUT_POINTS = OUT_DIR / "inclusivejet_fig6_canonical_jet8_policy_points.csv"
OUT_MANIFEST = OUT_DIR / "inclusivejet_fig6_canonical_jet8_policy_points_manifest.json"
OUT_NOTE = OUT_DIR / "inclusivejet_fig6_canonical_jet8_policy_note.txt"

SAMPLES = ["jet8", "jet12", "jet20", "jet30", "jet40"]


def as_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def build_rows() -> list[dict[str, object]]:
    rows = read_rows(IN_POINTS)
    out: list[dict[str, object]] = []
    for row in rows:
        sample = row["sample"]
        treatment = (
            "THE94_canonical_visible_formula_no_hidden_PPG12_jet8_factor"
            if sample == "jet8"
            else "no_scale_true_period_combined_PPG12_parity"
        )
        reference_role = (
            "historical_PPG12_jet8_reference_exception"
            if sample == "jet8"
            else "true_period_combined_PPG12_reference"
        )
        out.append(
            {
                "sample": sample,
                "bin_low": as_float(row, "bin_low"),
                "bin_high": as_float(row, "bin_high"),
                "bin_center": as_float(row, "bin_center"),
                "ppg12_reference_value": as_float(row, "ppg12_true_period_value"),
                "ppg12_reference_error": as_float(row, "ppg12_true_period_error"),
                "the94_canonical_value": as_float(row, "current_kept_value"),
                "the94_canonical_error": as_float(row, "current_kept_error"),
                "ppg12_fit_value": as_float(row, "ppg12_fit_value"),
                "ppg12_reference_over_fit": as_float(row, "ppg12_true_period_over_fit"),
                "ppg12_reference_over_fit_error": as_float(row, "ppg12_true_period_over_fit_error"),
                "the94_canonical_over_fit": as_float(row, "current_kept_over_fit"),
                "the94_canonical_over_fit_error": as_float(row, "current_kept_over_fit_error"),
                "the94_canonical_over_ppg12_reference": as_float(
                    row, "current_over_ppg12_true_period"
                ),
                "current_display_scale_factor": 1.0,
                "sample_treatment": treatment,
                "ppg12_reference_role": reference_role,
            }
        )
    out.sort(key=lambda r: (SAMPLES.index(str(r["sample"])), float(r["bin_center"])))
    return out


def sample_summary(rows: list[dict[str, object]]) -> dict[str, object]:
    summary: dict[str, object] = {}
    for sample in SAMPLES:
        sample_rows = [r for r in rows if r["sample"] == sample]
        ratios = [float(r["the94_canonical_over_ppg12_reference"]) for r in sample_rows]
        data_fit = [float(r["the94_canonical_over_fit"]) for r in sample_rows]
        ppg12_data_fit = [float(r["ppg12_reference_over_fit"]) for r in sample_rows]
        summary[sample] = {
            "n_bins": len(sample_rows),
            "median_THE94_over_PPG12_reference": median(ratios),
            "min_THE94_over_PPG12_reference": min(ratios),
            "max_THE94_over_PPG12_reference": max(ratios),
            "median_THE94_over_fit": median(data_fit),
            "median_PPG12_reference_over_fit": median(ppg12_data_fit),
            "treatment": sample_rows[0]["sample_treatment"] if sample_rows else "none",
        }
    return summary


def boundary_summary(rows: list[dict[str, object]]) -> dict[str, object]:
    jet8 = [r for r in rows if r["sample"] == "jet8"]
    jet12 = [r for r in rows if r["sample"] == "jet12"]
    last_jet8 = max(jet8, key=lambda r: float(r["bin_center"]))
    first_jet12 = min(jet12, key=lambda r: float(r["bin_center"]))
    return {
        "last_jet8_bin": [last_jet8["bin_low"], last_jet8["bin_high"]],
        "first_jet12_bin": [first_jet12["bin_low"], first_jet12["bin_high"]],
        "THE94_last_jet8_over_fit": last_jet8["the94_canonical_over_fit"],
        "THE94_first_jet12_over_fit": first_jet12["the94_canonical_over_fit"],
        "THE94_boundary_ratio_last_jet8_to_first_jet12": float(
            last_jet8["the94_canonical_over_fit"]
        )
        / float(first_jet12["the94_canonical_over_fit"]),
        "PPG12_last_jet8_over_fit": last_jet8["ppg12_reference_over_fit"],
        "PPG12_first_jet12_over_fit": first_jet12["ppg12_reference_over_fit"],
        "PPG12_boundary_ratio_last_jet8_to_first_jet12": float(
            last_jet8["ppg12_reference_over_fit"]
        )
        / float(first_jet12["ppg12_reference_over_fit"]),
    }


def write_rows(rows: list[dict[str, object]]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with OUT_POINTS.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    rows = build_rows()
    write_rows(rows)
    manifest = {
        "status": "ok_canonical_jet8_policy_points",
        "interpretation": "canonical RecoilJets jet8 policy comparison; not PPG12-verbatim jet8 parity",
        "input_points_csv": str(IN_POINTS),
        "points_csv": str(OUT_POINTS),
        "current_display_scale_factor": 1.0,
        "current_jet8_policy": (
            "Use the corrected THE-94 nominal jet8 output directly: visible "
            "jet8 xsec/jet50 factor, period lumi, SI/DI mix, and truth-vertex "
            "weights; do not apply the unresolved PPG12 hidden 0.4274 jet8 factor."
        ),
        "ppg12_jet8_status": (
            "Saved PPG12 jet8 true-period reference is retained as a historical "
            "source/normalization exception.  Jet12/20/30/40 are no-scale "
            "true-period-combined PPG12 parity references."
        ),
        "sample_summary": sample_summary(rows),
        "boundary_summary": boundary_summary(rows),
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    OUT_NOTE.write_text(
        "\n".join(
            [
                "THE-94 canonical jet8 policy plot input.",
                "",
                "This artifact preserves the corrected THE-94 nominal jet8 output with no display scale.",
                "It intentionally does not force jet8 to the saved PPG12 jet8 reference.",
                "The saved PPG12 jet8 points are treated as a historical source/normalization exception.",
                "Jet12/20/30/40 remain the no-scale true-period-combined PPG12 parity proof.",
                "",
                "Boundary check:",
                json.dumps(manifest["boundary_summary"], indent=2, sort_keys=True),
                "",
            ]
        )
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

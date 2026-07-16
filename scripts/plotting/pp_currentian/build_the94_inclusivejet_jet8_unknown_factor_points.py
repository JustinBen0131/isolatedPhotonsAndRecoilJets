#!/usr/bin/env python3
"""Build THE-94 inclusive-jet Fig.6 points with jet8 hidden-factor alignment.

This is a diagnostic/presentation table only.  It starts from the audited
true-period-combined no-scale THE-94 vs PPG12 Fig.6 points and multiplies only
the displayed THE-94 jet8 values by the unresolved PPG12 hidden-source factor.
The nominal THE-94 values are preserved in separate columns.
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
OUT_DIR = THE94_BASE / "jet8_unknown_factor_overlay"
OUT_POINTS = OUT_DIR / "inclusivejet_fig6_jet8_unknown_factor_points.csv"
OUT_MANIFEST = OUT_DIR / "inclusivejet_fig6_jet8_unknown_factor_points_manifest.json"
OUT_NOTE = OUT_DIR / "inclusivejet_fig6_jet8_unknown_factor_note.txt"

SAMPLES = ["jet8", "jet12", "jet20", "jet30", "jet40"]

# Audited median PPG12/current ratio for the true-period jet8 bins.  This is
# not a nominal physics correction; it is the unresolved PPG12 saved-reference
# normalization factor isolated in THE-94 jet8 provenance.
JET8_UNKNOWN_FACTOR = 0.42738368246379005


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
        display_scale = JET8_UNKNOWN_FACTOR if sample == "jet8" else 1.0
        nominal_value = as_float(row, "current_kept_value")
        nominal_error = as_float(row, "current_kept_error")
        display_value = nominal_value * display_scale
        display_error = nominal_error * display_scale
        fit_value = as_float(row, "ppg12_fit_value")
        ppg12_value = as_float(row, "ppg12_true_period_value")
        ppg12_error = as_float(row, "ppg12_true_period_error")
        out.append(
            {
                "sample": sample,
                "bin_low": as_float(row, "bin_low"),
                "bin_high": as_float(row, "bin_high"),
                "bin_center": as_float(row, "bin_center"),
                "ppg12_reference_value": ppg12_value,
                "ppg12_reference_error": ppg12_error,
                "current_nominal_value": nominal_value,
                "current_nominal_error": nominal_error,
                "current_display_value": display_value,
                "current_display_error": display_error,
                "ppg12_fit_value": fit_value,
                "ppg12_reference_over_fit": ppg12_value / fit_value,
                "ppg12_reference_over_fit_error": ppg12_error / fit_value,
                "current_display_over_fit": display_value / fit_value,
                "current_display_over_fit_error": display_error / fit_value,
                "current_nominal_over_ppg12_reference": as_float(
                    row, "current_over_ppg12_true_period"
                ),
                "current_display_over_ppg12_reference": display_value / ppg12_value,
                "current_display_scale_factor": display_scale,
                "sample_treatment": (
                    "display_only_jet8_times_unresolved_ppg12_factor"
                    if sample == "jet8"
                    else "nominal_no_scale"
                ),
            }
        )
    out.sort(key=lambda r: (SAMPLES.index(str(r["sample"])), float(r["bin_center"])))
    return out


def sample_summary(rows: list[dict[str, object]]) -> dict[str, object]:
    summary: dict[str, object] = {}
    for sample in SAMPLES:
        sample_rows = [r for r in rows if r["sample"] == sample]
        nominal = [float(r["current_nominal_over_ppg12_reference"]) for r in sample_rows]
        display = [float(r["current_display_over_ppg12_reference"]) for r in sample_rows]
        summary[sample] = {
            "n_bins": len(sample_rows),
            "median_nominal_current_over_ppg12": median(nominal),
            "median_display_current_over_ppg12": median(display),
            "min_display_current_over_ppg12": min(display),
            "max_display_current_over_ppg12": max(display),
            "display_scale_factor": sample_rows[0]["current_display_scale_factor"],
            "treatment": sample_rows[0]["sample_treatment"],
        }
    return summary


def write_rows(rows: list[dict[str, object]]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with OUT_POINTS.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    rows = build_rows()
    write_rows(rows)
    manifest = {
        "status": "ok_jet8_unknown_factor_points",
        "interpretation": (
            "Diagnostic historical-source-aligned display: THE-94 jet8 is "
            "multiplied by the unresolved PPG12 hidden factor while jet12-40 "
            "remain nominal no-scale."
        ),
        "input_points_csv": str(IN_POINTS),
        "points_csv": str(OUT_POINTS),
        "jet8_display_factor": JET8_UNKNOWN_FACTOR,
        "jet8_factor_origin": (
            "Audited median PPG12/current true-period jet8 ratio; exact PPG12 "
            "production source for this factor remains unknown."
        ),
        "not_a_nominal_correction": True,
        "sample_summary": sample_summary(rows),
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    OUT_NOTE.write_text(
        "\n".join(
            [
                "THE-94 inclusive-jet Fig.6 jet8 unknown-factor diagnostic.",
                "",
                f"Input points: {IN_POINTS}",
                f"Output points: {OUT_POINTS}",
                f"Jet8 display factor: {JET8_UNKNOWN_FACTOR:.17g}",
                "",
                "Only THE-94 jet8 display values/errors are multiplied by this factor.",
                "THE-94 jet12/20/30/40 are not scaled.",
                "Nominal THE-94 jet8 values are preserved in the CSV.",
                "Interpretation: this proves the residual jet8 mismatch is a flat scalar alignment; it is not a nominal no-scale parity claim.",
                "",
            ]
        )
    )
    print(OUT_POINTS)
    print(OUT_MANIFEST)
    print(OUT_NOTE)
    print(json.dumps(manifest["sample_summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

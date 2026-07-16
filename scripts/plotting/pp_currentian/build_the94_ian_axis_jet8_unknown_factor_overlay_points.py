#!/usr/bin/env python3
"""Build IAN-axis inclusive-jet overlay points with current jet8 displayed.

This is the same plot grammar as the THE-94 IAN-axis jet12-40 overlay, but it
also draws current jet8 after applying the unresolved jet8 alignment needed to
match the PPG12 SDCC source points on the IAN/source-scope canvas.

The nominal THE-94 values are preserved in source columns.  This is a display
diagnostic, not a nominal physics-output correction.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from statistics import median


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
THE94_BASE = (
    REPO
    / "dataOutput/ppg12Parity/the94_ppg12_inclusivejet_fig6_fixed_20260702_204032"
)
PPG12_IAN_CSV = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024/"
    "strict_stitched_inclusivejet/jet_sdcc_over_current_shape_overlay_ppg12_jet8_xsecfix_points.csv"
)
THE94_TRUE_PERIOD_CSV = (
    THE94_BASE
    / "no_scale_overlay/inclusivejet_fig6_true_period_combined_no_scale_data_over_fit_points.csv"
)
OUT_DIR = THE94_BASE / "ian_axis_jet8_unknown_factor_overlay"
OUT_POINTS = OUT_DIR / "inclusivejet_fig6_ian_axis_ppg12_the94_jet8_unknown_factor_points.csv"
OUT_MANIFEST = OUT_DIR / "inclusivejet_fig6_ian_axis_ppg12_the94_jet8_unknown_factor_points_manifest.json"

SAMPLES = ["jet8", "jet12", "jet20", "jet30", "jet40"]
NONJET8_SAMPLES = {"jet12", "jet20", "jet30", "jet40"}
JET8_TRUE_PERIOD_UNKNOWN_FACTOR = 0.42738368246379005


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def key(row: dict[str, str]) -> tuple[str, float]:
    return row["sample"], round(f(row, "bin_center"), 6)


def finite_or_blank(value: float) -> float | str:
    return value if math.isfinite(value) else ""


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ppg12_rows = {
        key(row): row
        for row in read_rows(PPG12_IAN_CSV)
        if row["sample"] in SAMPLES and f(row, "bin_center") <= 50.0
    }
    the94_rows = {
        key(row): row
        for row in read_rows(THE94_TRUE_PERIOD_CSV)
        if row["sample"] in SAMPLES and f(row, "bin_center") <= 50.0
    }
    common = sorted(
        set(ppg12_rows) & set(the94_rows),
        key=lambda item: (SAMPLES.index(item[0]), item[1]),
    )
    if not common:
        raise RuntimeError("no common rows between PPG12 IAN and THE-94 tables")

    nonjet8_display_ratios = []
    jet8_display_ratios = []
    for item in common:
        sample, _ = item
        ppg12 = ppg12_rows[item]
        current = the94_rows[item]
        cur_value = f(current, "current_kept_value")
        if cur_value <= 0.0:
            continue
        ratio = f(ppg12, "ppg12_sdcc_value") / cur_value
        if sample in NONJET8_SAMPLES:
            nonjet8_display_ratios.append(ratio)
        elif sample == "jet8":
            jet8_display_ratios.append(ratio)
    if not nonjet8_display_ratios or not jet8_display_ratios:
        raise RuntimeError("cannot derive display factors")

    common_nonjet8_factor = median(nonjet8_display_ratios)
    jet8_ian_axis_factor = median(jet8_display_ratios)

    out_rows: list[dict[str, object]] = []
    for item in common:
        sample, _ = item
        ppg12 = ppg12_rows[item]
        current = the94_rows[item]
        fit = f(ppg12, "ppg12_fit_value")
        nominal_value = f(current, "current_kept_value")
        nominal_error = f(current, "current_kept_error")
        display_factor = jet8_ian_axis_factor if sample == "jet8" else common_nonjet8_factor
        display_value = nominal_value * display_factor
        display_error = nominal_error * display_factor
        out_rows.append(
            {
                "sample": sample,
                "bin_low": f(ppg12, "bin_low"),
                "bin_high": f(ppg12, "bin_high"),
                "bin_center": f(ppg12, "bin_center"),
                "ppg12_ian_value": f(ppg12, "ppg12_sdcc_value"),
                "ppg12_ian_error": f(ppg12, "ppg12_sdcc_error"),
                "ppg12_fit_value": fit,
                "ppg12_ian_over_fit": f(ppg12, "ppg12_sdcc_over_fit"),
                "ppg12_ian_over_fit_error": f(ppg12, "ppg12_sdcc_error") / fit,
                "the94_nominal_value": nominal_value,
                "the94_nominal_error": nominal_error,
                "the94_nominal_over_ppg12_ian": nominal_value / f(ppg12, "ppg12_sdcc_value"),
                "the94_display_value": display_value,
                "the94_display_error": display_error,
                "the94_display_over_fit": display_value / fit,
                "the94_display_over_fit_error": display_error / fit,
                "the94_display_over_ppg12_ian": display_value / f(ppg12, "ppg12_sdcc_value"),
                "draw_the94": 1,
                "source_scope_display_factor": finite_or_blank(display_factor),
                "jet8_true_period_unknown_factor": (
                    JET8_TRUE_PERIOD_UNKNOWN_FACTOR if sample == "jet8" else ""
                ),
                "treatment": (
                    "jet8_display_aligned_to_ppg12_ian_unknown_factor"
                    if sample == "jet8"
                    else "common_nonjet8_source_scope_projection"
                ),
            }
        )

    with OUT_POINTS.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)

    summary: dict[str, object] = {}
    for sample in SAMPLES:
        rows = [r for r in out_rows if r["sample"] == sample]
        display = [float(r["the94_display_over_ppg12_ian"]) for r in rows]
        summary[sample] = {
            "n_bins": len(rows),
            "median_display_the94_over_ppg12_ian": median(display),
            "min_display_the94_over_ppg12_ian": min(display),
            "max_display_the94_over_ppg12_ian": max(display),
            "display_factor": rows[0]["source_scope_display_factor"] if rows else None,
            "treatment": rows[0]["treatment"] if rows else None,
        }

    manifest = {
        "status": "ok_ian_axis_jet8_unknown_factor_points",
        "interpretation": (
            "PPG12 IAN/source-scope Fig.6 axis and fit reproduction with THE-94 "
            "current output overlaid for jet8/12/20/30/40. Current jet8 is "
            "display-aligned to PPG12 SDCC using the unresolved jet8 factor."
        ),
        "ppg12_ian_csv": str(PPG12_IAN_CSV),
        "the94_true_period_csv": str(THE94_TRUE_PERIOD_CSV),
        "points_csv": str(OUT_POINTS),
        "common_nonjet8_source_scope_display_factor": common_nonjet8_factor,
        "jet8_ian_axis_display_factor": jet8_ian_axis_factor,
        "jet8_true_period_unknown_factor": JET8_TRUE_PERIOD_UNKNOWN_FACTOR,
        "jet8_note": (
            "The IAN-axis jet8 display factor is median(PPG12 IAN source / "
            "nominal THE-94 true-period current) over jet8 bins. This is the "
            "source-scope version of the unresolved jet8 alignment, not a "
            "nominal production correction."
        ),
        "sample_summary": summary,
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(OUT_POINTS)
    print(OUT_MANIFEST)
    print(f"common_nonjet8_source_scope_display_factor={common_nonjet8_factor:.12g}")
    print(f"jet8_ian_axis_display_factor={jet8_ian_axis_factor:.12g}")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

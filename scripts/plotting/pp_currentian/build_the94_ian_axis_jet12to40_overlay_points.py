#!/usr/bin/env python3
"""Build PPG12-IAN-axis inclusive-jet overlay points for THE-94.

The output table is for one presentation plot only.  It keeps the PPG12 Fig. 6
IAN/source-scope values and fit exactly as the displayed reference.  THE-94 is
drawn only for jet12/20/30/40, projected onto the IAN source-scope with one
common non-jet8 display factor.  THE-94 jet8 is intentionally not drawn.

Nominal THE-94 values are retained in the CSV.  The display factor is not a
physics-output correction and is recorded in the manifest.
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
OUT_DIR = THE94_BASE / "ian_axis_jet12to40_overlay"
OUT_POINTS = OUT_DIR / "inclusivejet_fig6_ian_axis_ppg12_jet8to40_the94_jet12to40_points.csv"
OUT_MANIFEST = OUT_DIR / "inclusivejet_fig6_ian_axis_ppg12_jet8to40_the94_jet12to40_points_manifest.json"

SAMPLES = ["jet8", "jet12", "jet20", "jet30", "jet40"]
THE94_DRAW_SAMPLES = {"jet12", "jet20", "jet30", "jet40"}


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

    display_ratios = []
    for item in common:
        sample, _ = item
        if sample not in THE94_DRAW_SAMPLES:
            continue
        ppg12 = ppg12_rows[item]
        current = the94_rows[item]
        cur_value = f(current, "current_kept_value")
        if cur_value > 0.0:
            display_ratios.append(f(ppg12, "ppg12_sdcc_value") / cur_value)
    if not display_ratios:
        raise RuntimeError("cannot derive non-jet8 source-scope display factor")
    common_display_factor = median(display_ratios)

    out_rows: list[dict[str, object]] = []
    for item in common:
        sample, _ = item
        ppg12 = ppg12_rows[item]
        current = the94_rows[item]
        draw_current = sample in THE94_DRAW_SAMPLES
        fit = f(ppg12, "ppg12_fit_value")
        nominal_value = f(current, "current_kept_value")
        nominal_error = f(current, "current_kept_error")
        display_factor = common_display_factor if draw_current else float("nan")
        display_value = nominal_value * display_factor if draw_current else float("nan")
        display_error = nominal_error * display_factor if draw_current else float("nan")
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
                "the94_display_value": finite_or_blank(display_value),
                "the94_display_error": finite_or_blank(display_error),
                "the94_display_over_fit": finite_or_blank(display_value / fit if draw_current else float("nan")),
                "the94_display_over_fit_error": finite_or_blank(display_error / fit if draw_current else float("nan")),
                "the94_display_over_ppg12_ian": finite_or_blank(
                    display_value / f(ppg12, "ppg12_sdcc_value") if draw_current else float("nan")
                ),
                "draw_the94": int(draw_current),
                "source_scope_display_factor": finite_or_blank(display_factor),
                "treatment": (
                    "common_nonjet8_source_scope_projection"
                    if draw_current
                    else "ppg12_reference_only_no_the94_jet8"
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
        display = [
            float(r["the94_display_over_ppg12_ian"])
            for r in rows
            if r["draw_the94"] and r["the94_display_over_ppg12_ian"] != ""
        ]
        summary[sample] = {
            "n_bins": len(rows),
            "draw_the94": sample in THE94_DRAW_SAMPLES,
            "median_display_the94_over_ppg12_ian": median(display) if display else None,
            "min_display_the94_over_ppg12_ian": min(display) if display else None,
            "max_display_the94_over_ppg12_ian": max(display) if display else None,
        }

    manifest = {
        "status": "ok_ian_axis_jet12to40_points",
        "interpretation": (
            "PPG12 IAN/source-scope Fig.6 axis and fit reproduction with THE-94 "
            "overlaid only for jet12/20/30/40. THE-94 jet8 is not drawn."
        ),
        "ppg12_ian_csv": str(PPG12_IAN_CSV),
        "the94_true_period_csv": str(THE94_TRUE_PERIOD_CSV),
        "points_csv": str(OUT_POINTS),
        "the94_draw_samples": sorted(THE94_DRAW_SAMPLES),
        "jet8_treatment": "PPG12 reference is shown; THE-94 jet8 is not drawn",
        "source_scope_display_factor_common_nonjet8": common_display_factor,
        "display_factor_note": (
            "Single common factor PPG12_IAN_source / THE94_true_period_current, "
            "derived from all jet12/20/30/40 bins. This is a display projection "
            "onto the IAN source-scope, not a physics-output correction."
        ),
        "sample_summary": summary,
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(OUT_POINTS)
    print(OUT_MANIFEST)
    print(f"common_nonjet8_source_scope_display_factor={common_display_factor:.12g}")


if __name__ == "__main__":
    main()

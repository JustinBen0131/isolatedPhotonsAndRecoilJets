#!/usr/bin/env python3
"""Build IAN-source-scope overlay points for THE-94 inclusive-jet Fig.6.

The PPG12 IAN/no-suffix source and the THE-94 true-period campaign are different
normalization targets.  This builder preserves nominal THE-94 values in source
columns and adds separate display columns that project the THE-94 current points
onto the IAN source-scope visual target using per-sample median PPG12/current
factors.  It is for plotting only; it does not alter campaign ROOT outputs.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from statistics import median


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
BASE = REPO / "dataOutput/ppg12Parity/the94_ppg12_inclusivejet_fig6_fixed_20260702_204032"
PPG12_SOURCE_CSV = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024/strict_stitched_inclusivejet/"
    / "jet_sdcc_over_current_shape_overlay_ppg12_jet8_xsecfix_points.csv"
)
THE94_CURRENT_VS_SOURCE_CSV = (
    BASE
    / "contract_canary/fig6_inclusive/fig6_inclusive_current_kept_vs_ppg12_source_points.csv"
)
OUT_DIR = BASE / "ian_source_overlay_the94_current"
OUT_POINTS = OUT_DIR / "inclusivejet_fig6_ian_source_overlay_the94_current_source_aligned_points.csv"
OUT_MANIFEST = OUT_DIR / "inclusivejet_fig6_ian_source_overlay_the94_current_source_aligned_points_manifest.json"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def key(row: dict[str, str]) -> tuple[str, float]:
    return (row["sample"], round(f(row, "bin_center"), 6))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ppg12_rows = {key(row): row for row in read_rows(PPG12_SOURCE_CSV) if f(row, "bin_center") <= 50.0}
    current_rows = {key(row): row for row in read_rows(THE94_CURRENT_VS_SOURCE_CSV) if f(row, "bin_center") <= 50.0}
    common_keys = sorted(
        set(ppg12_rows) & set(current_rows),
        key=lambda item: (["jet8", "jet12", "jet20", "jet30", "jet40"].index(item[0]), item[1]),
    )
    if not common_keys:
        raise RuntimeError("no common rows between PPG12 source and THE-94 current tables")

    factors: dict[str, float] = {}
    for sample in ["jet8", "jet12", "jet20", "jet30", "jet40"]:
        ratios = []
        for item in common_keys:
            if item[0] != sample:
                continue
            p = ppg12_rows[item]
            c = current_rows[item]
            cur = f(c, "current_kept_value")
            if cur > 0:
                ratios.append(f(p, "ppg12_sdcc_value") / cur)
        if ratios:
            factors[sample] = median(ratios)

    out_rows: list[dict[str, object]] = []
    for item in common_keys:
        sample, _ = item
        p = ppg12_rows[item]
        c = current_rows[item]
        factor = factors[sample]
        current_nominal = f(c, "current_kept_value")
        current_nominal_err = f(c, "current_kept_error")
        current_display = current_nominal * factor
        current_display_err = current_nominal_err * factor
        fit = f(p, "ppg12_fit_value")
        out_rows.append(
            {
                "sample": sample,
                "bin_low": f(p, "bin_low"),
                "bin_high": f(p, "bin_high"),
                "bin_center": f(p, "bin_center"),
                "ppg12_source_value": f(p, "ppg12_sdcc_value"),
                "ppg12_source_error": f(p, "ppg12_sdcc_error"),
                "ppg12_fit_value": fit,
                "ppg12_source_over_fit": f(p, "ppg12_sdcc_over_fit"),
                "ppg12_source_over_fit_error": f(p, "ppg12_sdcc_error") / fit,
                "current_nominal_value": current_nominal,
                "current_nominal_error": current_nominal_err,
                "current_nominal_over_ppg12_source": f(c, "current_kept_over_ppg12"),
                "current_display_value": current_display,
                "current_display_error": current_display_err,
                "current_display_over_fit": current_display / fit,
                "current_display_over_fit_error": current_display_err / fit,
                "current_display_over_ppg12_source": current_display / f(p, "ppg12_sdcc_value"),
                "source_scope_display_factor": factor,
                "treatment": "source_scope_display_alignment",
                "current_root": c["current_root"],
            }
        )

    fields = list(out_rows[0].keys())
    with OUT_POINTS.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(out_rows)

    summary = {}
    for sample, factor in factors.items():
        rows = [r for r in out_rows if r["sample"] == sample]
        summary[sample] = {
            "n_bins": len(rows),
            "display_factor": factor,
            "median_nominal_current_over_ppg12_source": median(
                float(r["current_nominal_over_ppg12_source"]) for r in rows
            ),
            "median_display_current_over_ppg12_source": median(
                float(r["current_display_over_ppg12_source"]) for r in rows
            ),
        }

    manifest = {
        "status": "ok_source_scope_overlay_points",
        "interpretation": "IAN-source-scope display comparison; not nominal no-scale parity",
        "ppg12_source_csv": str(PPG12_SOURCE_CSV),
        "the94_current_vs_source_csv": str(THE94_CURRENT_VS_SOURCE_CSV),
        "points_csv": str(OUT_POINTS),
        "display_treatment": (
            "THE-94 current values are preserved in nominal columns and displayed after "
            "per-sample median source-scope alignment to the PPG12 IAN/no-suffix source."
        ),
        "sample_summary": summary,
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(OUT_POINTS)
    print(OUT_MANIFEST)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Make a wide two-panel PPG12/current stitched spectra ratio summary.

The top panel summarizes the corrected photon+jet stitch comparison using the
strict no-scale PPG12-source/current point table. The bottom panel uses the
explicit source-exposure-corrected inclusive-jet table: the common jet12-40
source-count convention is removed, while the jet8-only residual is left visible.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt


BASE = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12Parity/"
    "the76_ppg12_parity_full_20260701_003024"
)
PHOTON_POINTS = (
    BASE
    / "strict_stitched_photonjet"
    / "photon_strict_july1_over_ppg12_source_points.csv"
)
JET_POINTS = (
    BASE
    / "strict_stitched_inclusivejet"
    / "jet_source_scope_common_exposure_corrected_points.csv"
)
OUT_DIR = BASE / "summary_stitched_ratio"
OUT_PNG = OUT_DIR / "ppg12_over_current_stitched_spectra_two_panel_summary_fixed_photon_corrected_inclusive.png"
OUT_CSV = OUT_DIR / "ppg12_over_current_stitched_spectra_two_panel_summary_fixed_photon_corrected_inclusive_points.csv"
OUT_MANIFEST = OUT_DIR / "ppg12_over_current_stitched_spectra_two_panel_summary_fixed_photon_corrected_inclusive_manifest.json"

PHOTON_COLORS = {
    "photon5": "#e83e8c",
    "photon10": "#2ca02c",
    "photon20": "#1da1f2",
}
JET_COLORS = {
    "jet8": "#e83e8c",
    "jet12": "#2ca02c",
    "jet20": "#1da1f2",
    "jet30": "#ff6f00",
    "jet40": "#d65ad1",
}


def _safe_float(row: dict[str, str], key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return float("nan")


def _ratio_error(num: float, num_err: float, den: float, den_err: float) -> float:
    if not all(math.isfinite(v) for v in (num, num_err, den, den_err)):
        return float("nan")
    if num <= 0 or den <= 0:
        return float("nan")
    ratio = num / den
    return ratio * math.hypot(num_err / num, den_err / den)


def read_photon() -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    with PHOTON_POINTS.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            sdcc = _safe_float(row, "ppg12_source_value")
            sdcc_err = _safe_float(row, "ppg12_source_error")
            current = _safe_float(row, "july1_value")
            current_err = _safe_float(row, "july1_error")
            if current <= 0:
                continue
            ratio = sdcc / current
            rows.append(
                {
                    "panel": "photon+jet",
                    "sample": row["sample"],
                    "bin_low": _safe_float(row, "bin_low"),
                    "bin_high": _safe_float(row, "bin_high"),
                    "bin_center": _safe_float(row, "bin_center"),
                    "ppg12_value": sdcc,
                    "ppg12_error": sdcc_err,
                    "current_value": current,
                    "current_error": current_err,
                    "ratio_ppg12_over_current": ratio,
                    "ratio_error": _ratio_error(sdcc, sdcc_err, current, current_err),
                    "current_definition": "July1 photon+jet strict no-scale output; PPG12 photon-slice ownership and xsec/N convention",
                }
            )
    return rows


def read_jet() -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    with JET_POINTS.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            sdcc = _safe_float(row, "ppg12_source_value")
            sdcc_err = _safe_float(row, "ppg12_source_error")
            current = _safe_float(row, "current_common_exposure_value")
            current_err = _safe_float(row, "current_common_exposure_error")
            if current <= 0:
                continue
            ratio = sdcc / current
            rows.append(
                {
                    "panel": "inclusive jet",
                    "sample": row["sample"],
                    "bin_low": _safe_float(row, "bin_low"),
                    "bin_high": _safe_float(row, "bin_high"),
                    "bin_center": _safe_float(row, "bin_center"),
                    "ppg12_value": sdcc,
                    "ppg12_error": sdcc_err,
                    "current_value": current,
                    "current_error": current_err,
                    "ratio_ppg12_over_current": ratio,
                    "ratio_error": _ratio_error(sdcc, sdcc_err, current, current_err),
                    "current_definition": "July1 inclusive-jet output after explicit common source-exposure correction from jet12-40; jet8 residual left visible",
                }
            )
    return rows


def write_points(rows: list[dict[str, float | str]]) -> None:
    fieldnames = [
        "panel",
        "sample",
        "bin_low",
        "bin_high",
        "bin_center",
        "ppg12_value",
        "ppg12_error",
        "current_value",
        "current_error",
        "ratio_ppg12_over_current",
        "ratio_error",
        "current_definition",
    ]
    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def draw_panel(
    ax,
    rows: list[dict[str, float | str]],
    colors: dict[str, str],
    title: str,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> None:
    ax.axhline(1.0, color="0.45", lw=1.5, ls=(0, (5, 5)), zorder=0)
    for sample, color in colors.items():
        sample_rows = [r for r in rows if r["sample"] == sample]
        if not sample_rows:
            continue
        x = [float(r["bin_center"]) for r in sample_rows]
        y = [float(r["ratio_ppg12_over_current"]) for r in sample_rows]
        ey = [float(r["ratio_error"]) for r in sample_rows]
        ax.errorbar(
            x,
            y,
            yerr=ey,
            fmt="o",
            ms=6.5,
            lw=1.2,
            elinewidth=1.1,
            capsize=0,
            color=color,
            markeredgecolor=color,
            markerfacecolor=color,
            label=sample,
        )
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_ylabel("PPG12 / Current", fontsize=24)
    ax.set_title(title, loc="left", fontsize=26, fontweight="bold", pad=8)
    ax.tick_params(axis="both", which="major", labelsize=20, direction="in", top=True, right=True, length=8, width=1.4)
    ax.tick_params(axis="both", which="minor", direction="in", top=True, right=True, length=4, width=1.0)
    ax.minorticks_on()
    ax.legend(
        loc="upper right",
        ncol=len(colors),
        frameon=False,
        fontsize=17,
        handletextpad=0.4,
        columnspacing=1.0,
        borderpad=0.1,
    )
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    photon = read_photon()
    jet = read_jet()
    all_rows = photon + jet
    write_points(all_rows)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.unicode_minus": False,
        }
    )
    fig, axes = plt.subplots(2, 1, figsize=(16, 9), sharex=False, sharey=False)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.93, bottom=0.105, hspace=0.28)

    draw_panel(axes[0], photon, PHOTON_COLORS, "photon+jet stitched spectra", (10, 40), (0.94, 1.06))
    draw_panel(axes[1], jet, JET_COLORS, "inclusive jet stitched spectra, source-count corrected", (9, 50), (0.86, 1.08))
    axes[0].tick_params(labelbottom=False)
    axes[1].set_xlabel("Leading object $p_T$ or $E_T$ [GeV]", fontsize=26)

    fig.savefig(OUT_PNG, dpi=160)
    plt.close(fig)

    ratios = [float(r["ratio_ppg12_over_current"]) for r in all_rows]
    manifest = {
        "status": "ok",
        "plot_png": str(OUT_PNG),
        "points_csv": str(OUT_CSV),
        "photon_points_csv": str(PHOTON_POINTS),
        "inclusive_jet_points_csv": str(JET_POINTS),
        "ratio_definition": "PPG12 SDCC stitched spectrum divided by current analysis stitched spectrum",
        "photon_current_definition": "July1 photon+jet strict no-scale output from photon_strict_july1_over_ppg12_source_points.csv",
        "inclusive_jet_current_definition": "July1 inclusive-jet source-exposure corrected output from jet_source_scope_common_exposure_corrected_points.csv; correction derived from jet12-40 only, with jet8 residual retained",
        "canvas_pixels": [2560, 1440],
        "panels": ["photon+jet top", "inclusive jet bottom"],
        "y_ranges": {"photon+jet": [0.94, 1.06], "inclusive jet": [0.86, 1.08]},
        "n_points": len(all_rows),
        "ratio_min": min(ratios),
        "ratio_max": max(ratios),
        "important_caveat": "Inclusive jet is corrected only for the common source-count convention; jet8 remains visibly non-closed and is not used to derive the correction.",
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(OUT_PNG)
    print(OUT_CSV)
    print(OUT_MANIFEST)


if __name__ == "__main__":
    main()

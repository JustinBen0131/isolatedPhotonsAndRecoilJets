#!/usr/bin/env python3
"""Make the slide-style stitched spectra ratio closure summary.

Top panel:
  Latest strict pre-vz photon+jet stitch check, unchanged.

Bottom panel:
  THE-94 inclusive-jet IAN-axis display comparison. Jet8 uses the explicit
  unresolved alignment factor recorded in the input CSV; jet12-40 use the
  common source-frame projection from the same overlay product.

This is a plotting product only. It does not modify nominal RecoilJets output
or promote the jet8 factor as a physics correction.
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import median

import matplotlib.pyplot as plt


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")

PHOTON_POINTS = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig5_prevz_strict_20260705_224243"
    / "strict_stitched_photonjet_prevz"
    / "photon_data_over_fit_sdcc_vs_current_overlay_prevz_strict_points.csv"
)

THE94_BASE = (
    REPO
    / "dataOutput/ppg12Parity/the94_ppg12_inclusivejet_fig6_fixed_20260702_204032"
)

INCLUSIVE_POINTS = (
    THE94_BASE
    / "ian_axis_jet8_unknown_factor_overlay"
    / "inclusivejet_fig6_ian_axis_ppg12_the94_jet8_unknown_factor_points.csv"
)

OUT_DIR = THE94_BASE / "summary_stitched_ratio_jet8_aligned"
OUT_PNG = OUT_DIR / "stitched_spectra_ratio_closure_photon_prevz_inclusive_jet8_aligned.png"
OUT_CSV = OUT_DIR / "stitched_spectra_ratio_closure_photon_prevz_inclusive_jet8_aligned_points.csv"
OUT_MANIFEST = OUT_DIR / "stitched_spectra_ratio_closure_photon_prevz_inclusive_jet8_aligned_manifest.json"

COLORS = {
    "photon5": "#e83e8c",
    "photon10": "#2ca02c",
    "photon20": "#1da1f2",
    "jet8": "#e83e8c",
    "jet12": "#2ca02c",
    "jet20": "#1da1f2",
    "jet30": "#ff6f00",
    "jet40": "#d65ad1",
}


def _f(row: dict[str, str], key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return math.nan


def _ratio_error(num: float, num_err: float, den: float, den_err: float) -> float:
    if not all(math.isfinite(v) for v in (num, num_err, den, den_err)):
        return math.nan
    if num <= 0 or den <= 0:
        return math.nan
    return (num / den) * math.hypot(num_err / num, den_err / den)


def read_photon_rows() -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    with PHOTON_POINTS.open() as f:
        for row in csv.DictReader(f):
            ppg12 = _f(row, "ppg12_sdcc_value")
            ppg12_err = _f(row, "ppg12_sdcc_error")
            current = _f(row, "current_value")
            current_err = _f(row, "current_error")
            if not (ppg12 > 0 and current > 0):
                continue
            rows.append(
                {
                    "panel": "photon+jet",
                    "sample": row["sample"],
                    "bin_low": _f(row, "bin_low"),
                    "bin_high": _f(row, "bin_high"),
                    "bin_center": _f(row, "bin_center"),
                    "ppg12_value": ppg12,
                    "ppg12_error": ppg12_err,
                    "current_value": current,
                    "current_error": current_err,
                    "ratio_ppg12_over_current": ppg12 / current,
                    "ratio_error": _ratio_error(ppg12, ppg12_err, current, current_err),
                    "treatment": "strict_prevz_photon_top_panel_unchanged",
                    "source_csv": str(PHOTON_POINTS),
                }
            )
    return rows


def read_inclusive_rows() -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    with INCLUSIVE_POINTS.open() as f:
        for row in csv.DictReader(f):
            if row.get("draw_the94") not in {"1", "true", "True"}:
                continue
            ppg12 = _f(row, "ppg12_ian_value")
            ppg12_err = _f(row, "ppg12_ian_error")
            current = _f(row, "the94_display_value")
            current_err = _f(row, "the94_display_error")
            if not (ppg12 > 0 and current > 0):
                continue
            rows.append(
                {
                    "panel": "inclusive jet",
                    "sample": row["sample"],
                    "bin_low": _f(row, "bin_low"),
                    "bin_high": _f(row, "bin_high"),
                    "bin_center": _f(row, "bin_center"),
                    "ppg12_value": ppg12,
                    "ppg12_error": ppg12_err,
                    "current_value": current,
                    "current_error": current_err,
                    "ratio_ppg12_over_current": ppg12 / current,
                    "ratio_error": _ratio_error(ppg12, ppg12_err, current, current_err),
                    "treatment": row.get("treatment", ""),
                    "source_csv": str(INCLUSIVE_POINTS),
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
        "treatment",
        "source_csv",
    ]
    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def draw_panel(
    ax,
    rows: list[dict[str, float | str]],
    samples: list[str],
    title: str,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> None:
    ax.axhline(1.0, color="0.50", lw=1.0, ls=(0, (5, 5)), zorder=0)
    for sample in samples:
        pts = [r for r in rows if r["sample"] == sample]
        if not pts:
            continue
        ax.errorbar(
            [float(r["bin_center"]) for r in pts],
            [float(r["ratio_ppg12_over_current"]) for r in pts],
            yerr=[float(r["ratio_error"]) for r in pts],
            fmt="o",
            ms=4.6,
            lw=1.0,
            elinewidth=1.0,
            capsize=0,
            color=COLORS[sample],
            markeredgecolor=COLORS[sample],
            markerfacecolor=COLORS[sample],
            label=sample,
        )
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_ylabel("PPG12 / Current", fontsize=16)
    ax.set_title(title, loc="left", fontsize=20, fontweight="bold", pad=6)
    ax.tick_params(axis="both", which="major", labelsize=12, direction="in", top=True, right=True, length=6, width=0.9)
    ax.tick_params(axis="both", which="minor", direction="in", top=True, right=True, length=3, width=0.8)
    ax.minorticks_on()
    ax.legend(
        loc="upper right",
        ncol=len(samples),
        frameon=False,
        fontsize=14.5,
        handletextpad=0.35,
        columnspacing=0.8,
        borderpad=0.1,
    )
    for spine in ax.spines.values():
        spine.set_linewidth(0.9)


def _summary(rows: list[dict[str, float | str]]) -> dict[str, dict[str, float | int]]:
    by_sample: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        by_sample[str(row["sample"])].append(float(row["ratio_ppg12_over_current"]))
    return {
        sample: {
            "n_bins": len(vals),
            "min": min(vals),
            "median": median(vals),
            "max": max(vals),
        }
        for sample, vals in sorted(by_sample.items())
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    photon = read_photon_rows()
    inclusive = read_inclusive_rows()
    rows = photon + inclusive
    write_points(rows)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.unicode_minus": False,
        }
    )

    fig = plt.figure(figsize=(13.333, 7.5), dpi=180)
    fig.patch.set_facecolor("white")

    ax_top = fig.add_axes([0.105, 0.565, 0.865, 0.385])
    ax_bottom = fig.add_axes([0.105, 0.105, 0.865, 0.385])
    draw_panel(
        ax_top,
        photon,
        ["photon5", "photon10", "photon20"],
        "photon+jet stitched spectra, strict pre-vz diagnostic",
        (10, 40),
        (0.94, 1.06),
    )
    draw_panel(
        ax_bottom,
        inclusive,
        ["jet8", "jet12", "jet20", "jet30", "jet40"],
        "inclusive jet stitched spectra, jet8 display-aligned",
        (9, 50),
        (0.94, 1.06),
    )
    ax_top.tick_params(labelbottom=False)
    ax_bottom.set_xlabel("Leading object $p_T$ or $E_T$ [GeV]", fontsize=16)

    fig.savefig(OUT_PNG)
    plt.close(fig)

    manifest = {
        "status": "ok",
        "plot_png": str(OUT_PNG),
        "points_csv": str(OUT_CSV),
        "photon_input_csv": str(PHOTON_POINTS),
        "inclusive_input_csv": str(INCLUSIVE_POINTS),
        "ratio_definition": "PPG12 / Current",
        "layout": "two ratio panels only; no main slide title or page number; enlarged axes fill the PNG canvas",
        "top_panel": "latest strict pre-vz photon+jet source-stage output; unchanged from photon_data_over_fit_sdcc_vs_current_overlay_prevz_strict_points.csv",
        "bottom_panel": "THE-94 inclusive IAN-axis display comparison; jet8 explicitly display-aligned to PPG12, jet12-40 use common source-frame projection",
        "jet8_caveat": "Jet8 alignment is a plotting/display diagnostic for the unresolved PPG12 jet8 normalization factor, not a nominal no-scale production correction.",
        "inclusive_display_factors": {
            "nonjet8_source_scope_display_factor": 0.3324931587057557,
            "jet8_ian_axis_display_factor": 0.29252707808847955,
            "jet8_true_period_unknown_factor": 0.42738368246379005,
        },
        "ratio_summary": {
            "photon": _summary(photon),
            "inclusive": _summary(inclusive),
        },
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")

    print(OUT_PNG)
    print(OUT_CSV)
    print(OUT_MANIFEST)


if __name__ == "__main__":
    main()

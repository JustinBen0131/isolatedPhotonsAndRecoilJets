#!/usr/bin/env python3
"""Render a THE-94 inclusive-jet Fig.6 data/fit presentation overlay.

This intentionally does not modify the nominal THE-94 campaign output.  It
starts from the corrected no-scale, true-period THE-94 Fig.6 points and applies
one plot-only jet8 historical-source alignment factor to the *displayed* current
jet8 points.  The factor is the median PPG12/current ratio over the jet8 bins in
the already-audited no-scale table, and is recorded in the output CSV/manifest.
Jet12/20/30/40 are left unscaled.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
THE94_BASE = (
    REPO
    / "dataOutput/ppg12Parity/the94_ppg12_inclusivejet_fig6_fixed_20260702_204032"
)
NO_SCALE_DIR = THE94_BASE / "no_scale_overlay"
IN_POINTS_CSV = NO_SCALE_DIR / "inclusivejet_fig6_true_period_combined_no_scale_data_over_fit_points.csv"

OUT_DIR = THE94_BASE / "presentation_data_over_fit_overlay"
OUT_PNG = OUT_DIR / "inclusivejet_fig6_data_over_fit_overlay_jet8_historical_aligned.png"
OUT_POINTS = OUT_DIR / "inclusivejet_fig6_data_over_fit_overlay_jet8_historical_aligned_points.csv"
OUT_MANIFEST = OUT_DIR / "inclusivejet_fig6_data_over_fit_overlay_jet8_historical_aligned_manifest.json"
OUT_NOTE = OUT_DIR / "inclusivejet_fig6_data_over_fit_overlay_jet8_historical_aligned_note.txt"

SAMPLES = ["jet8", "jet12", "jet20", "jet30", "jet40"]
COLORS = {
    "jet8": "#d62aa0",
    "jet12": "#239b35",
    "jet20": "#169ce8",
    "jet30": "#ff6f00",
    "jet40": "#c12ac7",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def ratio_error(err: float, fit: float) -> float:
    if fit <= 0 or not all(math.isfinite(v) for v in (err, fit)):
        return float("nan")
    return err / fit


def derive_jet8_alignment(rows: list[dict[str, str]]) -> float:
    ratios = [
        f(row, "ppg12_true_period_value") / f(row, "current_kept_value")
        for row in rows
        if row["sample"] == "jet8" and f(row, "current_kept_value") > 0
    ]
    if not ratios:
        raise RuntimeError("cannot derive jet8 alignment factor: no jet8 rows")
    return float(np.median(ratios))


def build_fit(rows: list[dict[str, str]]) -> tuple[np.ndarray, np.ndarray]:
    x = np.array([f(row, "bin_center") for row in rows], dtype=float)
    y = np.array([f(row, "ppg12_true_period_value") for row in rows], dtype=float)
    mask = (x >= 9.0) & (x <= 50.0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    coeff = np.polyfit(np.log(x[mask]), np.log(y[mask]), deg=4)
    grid = np.linspace(9.0, 50.0, 600)
    curve = np.exp(np.polyval(coeff, np.log(grid)))
    return grid, curve


def build_rows() -> tuple[list[dict[str, object]], np.ndarray, np.ndarray, float]:
    source_rows = read_csv(IN_POINTS_CSV)
    jet8_factor = derive_jet8_alignment(source_rows)
    fit_grid, fit_curve = build_fit(source_rows)
    fit_by_x = {
        round(f(row, "bin_center"), 6): f(row, "ppg12_fit_value")
        for row in source_rows
    }

    rows: list[dict[str, object]] = []
    for row in source_rows:
        sample = row["sample"]
        center = round(f(row, "bin_center"), 6)
        fit_value = fit_by_x[center]
        ppg12 = f(row, "ppg12_true_period_value")
        ppg12_err = f(row, "ppg12_true_period_error")
        current_nominal = f(row, "current_kept_value")
        current_nominal_err = f(row, "current_kept_error")
        display_factor = jet8_factor if sample == "jet8" else 1.0
        current_display = current_nominal * display_factor
        current_display_err = current_nominal_err * display_factor
        treatment = (
            "plot_only_jet8_historical_source_alignment"
            if sample == "jet8"
            else "nominal_no_scale"
        )
        rows.append(
            {
                "sample": sample,
                "bin_low": f(row, "bin_low"),
                "bin_high": f(row, "bin_high"),
                "bin_center": center,
                "ppg12_true_period_value": ppg12,
                "ppg12_true_period_error": ppg12_err,
                "current_nominal_value": current_nominal,
                "current_nominal_error": current_nominal_err,
                "current_display_value": current_display,
                "current_display_error": current_display_err,
                "ppg12_fit_value": fit_value,
                "ppg12_true_period_over_fit": ppg12 / fit_value,
                "ppg12_true_period_over_fit_error": ratio_error(ppg12_err, fit_value),
                "current_display_over_fit": current_display / fit_value,
                "current_display_over_fit_error": ratio_error(current_display_err, fit_value),
                "current_nominal_over_ppg12_true_period": f(row, "current_over_ppg12_true_period"),
                "display_current_over_ppg12_true_period": current_display / ppg12,
                "jet8_display_alignment_factor": display_factor,
                "treatment": treatment,
            }
        )

    rows.sort(key=lambda r: (SAMPLES.index(str(r["sample"])), float(r["bin_center"])))
    return rows, fit_grid, fit_curve, jet8_factor


def write_points(rows: list[dict[str, object]]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample",
        "bin_low",
        "bin_high",
        "bin_center",
        "ppg12_true_period_value",
        "ppg12_true_period_error",
        "current_nominal_value",
        "current_nominal_error",
        "current_display_value",
        "current_display_error",
        "ppg12_fit_value",
        "ppg12_true_period_over_fit",
        "ppg12_true_period_over_fit_error",
        "current_display_over_fit",
        "current_display_over_fit_error",
        "current_nominal_over_ppg12_true_period",
        "display_current_over_ppg12_true_period",
        "jet8_display_alignment_factor",
        "treatment",
    ]
    with OUT_POINTS.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sample_summary(rows: list[dict[str, object]]) -> dict[str, object]:
    out: dict[str, object] = {}
    for sample in SAMPLES:
        sample_rows = [r for r in rows if r["sample"] == sample]
        nominal = [float(r["current_nominal_over_ppg12_true_period"]) for r in sample_rows]
        display = [float(r["display_current_over_ppg12_true_period"]) for r in sample_rows]
        out[sample] = {
            "n_bins": len(sample_rows),
            "median_nominal_current_over_ppg12": float(np.median(nominal)),
            "median_display_current_over_ppg12": float(np.median(display)),
            "min_display_current_over_ppg12": float(min(display)),
            "max_display_current_over_ppg12": float(max(display)),
            "treatment": str(sample_rows[0]["treatment"]) if sample_rows else "none",
        }
    return out


def render(rows: list[dict[str, object]], fit_grid: np.ndarray, fit_curve: np.ndarray) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.45,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 7.5,
            "ytick.major.size": 7.5,
            "xtick.minor.size": 4,
            "ytick.minor.size": 4,
        }
    )

    fig = plt.figure(figsize=(8.0, 8.89), dpi=200)
    top = fig.add_axes([0.13, 0.40, 0.79, 0.56])
    bot = fig.add_axes([0.13, 0.09, 0.79, 0.30], sharex=top)
    top.set_yscale("log")
    top.set_xlim(8.8, 50.2)
    top.set_ylim(5.0e3, 2.0e12)
    bot.set_ylim(0.84, 1.20)

    for sample in SAMPLES:
        pts = [r for r in rows if r["sample"] == sample]
        color = COLORS[sample]
        x = [float(r["bin_center"]) for r in pts]
        top.errorbar(
            x,
            [float(r["current_display_value"]) for r in pts],
            yerr=[float(r["current_display_error"]) for r in pts],
            fmt="o",
            ms=4.2,
            mfc=color,
            mec=color,
            mew=0.7,
            ecolor=color,
            elinewidth=0.55,
            linestyle="none",
            zorder=3,
        )
        top.errorbar(
            x,
            [float(r["ppg12_true_period_value"]) for r in pts],
            yerr=[float(r["ppg12_true_period_error"]) for r in pts],
            fmt="o",
            ms=5.6,
            mfc="none",
            mec=color,
            mew=1.25,
            ecolor=color,
            elinewidth=0.65,
            linestyle="none",
            zorder=4,
        )
        bot.errorbar(
            x,
            [float(r["current_display_over_fit"]) for r in pts],
            yerr=[float(r["current_display_over_fit_error"]) for r in pts],
            fmt="o",
            ms=3.6,
            mfc="black",
            mec="black",
            mew=0.6,
            ecolor="black",
            elinewidth=0.48,
            linestyle="none",
            zorder=3,
        )
        bot.errorbar(
            x,
            [float(r["ppg12_true_period_over_fit"]) for r in pts],
            yerr=[float(r["ppg12_true_period_over_fit_error"]) for r in pts],
            fmt="o",
            ms=5.0,
            mfc="none",
            mec="black",
            mew=1.05,
            ecolor="black",
            elinewidth=0.55,
            linestyle="none",
            zorder=4,
        )

    top.plot(fit_grid, fit_curve, color="red", lw=1.6, zorder=2)
    bot.axhline(1.0, color="0.45", lw=1.0, ls=(0, (5, 5)))

    top.set_ylabel("counts", fontsize=24)
    bot.set_ylabel("MC / Fit", fontsize=22)
    bot.set_xlabel(r"Leading $p_T^\mathrm{jet}$ [GeV]", fontsize=25)
    top.tick_params(which="both", labelbottom=False, labelsize=18)
    bot.tick_params(which="both", labelsize=18)
    top.minorticks_on()
    bot.minorticks_on()
    bot.set_yticks([0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20])

    top.text(
        0.55,
        0.97,
        r"$\bf{\it{sPHENIX}}$ Internal" + "\n" + r"$p$+$p$ $\sqrt{s}=200$ GeV" + "\nPYTHIA8",
        transform=top.transAxes,
        fontsize=21,
        va="top",
    )
    top.text(
        0.06,
        0.16,
        "jet8 current: plot-only\nhistorical-source aligned",
        transform=top.transAxes,
        fontsize=15.5,
        va="bottom",
    )

    sample_handles = [
        Line2D([0], [0], marker="o", color=COLORS[s], mfc=COLORS[s], lw=0, ms=8, label=s)
        for s in SAMPLES
    ]
    source_handles = [
        Line2D([0], [0], marker="o", color="black", mfc="none", lw=0, ms=8, label="PPG12 source"),
        Line2D([0], [0], marker="o", color="black", mfc="black", lw=0, ms=8, label="Current output"),
        Line2D([0], [0], color="red", lw=1.8, label="PPG12 fit guide"),
    ]
    leg1 = top.legend(
        handles=sample_handles,
        title="sample",
        frameon=False,
        fontsize=15,
        title_fontsize=15,
        loc="lower left",
        bbox_to_anchor=(0.07, 0.30),
        handletextpad=0.6,
        labelspacing=0.48,
    )
    top.add_artist(leg1)
    top.legend(
        handles=source_handles,
        title="source",
        frameon=False,
        fontsize=14,
        title_fontsize=14,
        loc="lower left",
        bbox_to_anchor=(0.37, 0.30),
        handletextpad=0.7,
        labelspacing=0.55,
    )

    fig.savefig(OUT_PNG)
    plt.close(fig)


def write_manifest(rows: list[dict[str, object]], jet8_factor: float) -> None:
    manifest = {
        "status": "ok_presentation_data_over_fit_overlay",
        "interpretation": "historical-source-aligned comparison; not nominal no-scale jet8 parity",
        "png": str(OUT_PNG),
        "points_csv": str(OUT_POINTS),
        "note": str(OUT_NOTE),
        "input_points_csv": str(IN_POINTS_CSV),
        "ppg12_reference": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet_bdt_nom.root::h_max_truth_jet_pT, extracted read-only as the true period-combined Fig.6 reference",
        "current_output": str(THE94_BASE / "final_roots/inclusivejet"),
        "fit_source": "same smooth degree-4 log-polynomial guide used by the THE-94 true-period no-scale data/fit renderer, computed from the plotted PPG12 true-period points",
        "jet8_treatment": {
            "description": "Plot-only historical-source alignment of current jet8 display points to the saved PPG12 jet8 reference; nominal current jet8 values remain unchanged in source columns.",
            "factor_applied_to_current_jet8_display_columns": jet8_factor,
            "derived_from": "median(PPG12 true-period value / nominal THE-94 current value) over the five jet8 bins in the audited no-scale points CSV",
            "reason": "THE-94 jet8 provenance audit isolated a flat saved-PPG12 jet8 source/normalization exception of about 0.4274; exact historical recipe is not reconstructable from available PPG12 configs/logs.",
            "scope": "this figure only",
        },
        "sample_summary": sample_summary(rows),
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def write_note(jet8_factor: float) -> None:
    text = (
        "THE-94 inclusive-jet Fig.6 presentation data/fit overlay\n"
        "\n"
        f"Input points: {IN_POINTS_CSV}\n"
        f"Output PNG: {OUT_PNG}\n"
        f"Output CSV: {OUT_POINTS}\n"
        "\n"
        "Jet8 treatment: current jet8 display values and errors are multiplied by "
        f"{jet8_factor:.12g}. This is the median PPG12/current ratio over the "
        "audited no-scale jet8 bins and is used only to align the historical "
        "PPG12 jet8 source for this visual comparison. Nominal THE-94 output and "
        "jet12/20/30/40 values are not changed.\n"
        "\n"
        "Interpretation: historical-source-aligned comparison / presentation-mode "
        "corrected comparison, not nominal no-scale parity for jet8.\n"
    )
    OUT_NOTE.write_text(text)


def main() -> None:
    rows, fit_grid, fit_curve, jet8_factor = build_rows()
    write_points(rows)
    render(rows, fit_grid, fit_curve)
    write_manifest(rows, jet8_factor)
    write_note(jet8_factor)
    print(OUT_PNG)
    print(OUT_POINTS)
    print(OUT_MANIFEST)
    print(OUT_NOTE)


if __name__ == "__main__":
    main()

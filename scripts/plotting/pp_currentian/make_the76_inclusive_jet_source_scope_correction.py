#!/usr/bin/env python3
"""Build the THE-76 inclusive-jet source-scope correction artifacts.

PPG12 Fig. 6 is a weighted-count source diagnostic, not a cross-section
normalised spectrum. The July 1 current sample has a larger effective
period/component event population, so the raw all-component counts are high by
a nearly constant factor for jet12-40. This helper applies that common
non-jet8 source-exposure factor explicitly and leaves the jet8 residual visible.
"""

from __future__ import annotations

import csv
import json
import math
import statistics
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


BASE = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12Parity/"
    "the76_ppg12_parity_full_20260701_003024"
)
IN_CSV = BASE / "strict_stitched_inclusivejet" / "jet_component_exposure_audit_points.csv"
OUT_DIR = BASE / "strict_stitched_inclusivejet"
OUT_CSV = OUT_DIR / "jet_source_scope_common_exposure_corrected_points.csv"
OUT_PNG = OUT_DIR / "jet_data_over_fit_sdcc_vs_current_overlay_root_ppg12_style_source_scope_corrected.png"
OUT_MANIFEST = OUT_DIR / "jet_data_over_fit_sdcc_vs_current_overlay_root_ppg12_style_source_scope_corrected_manifest.json"

SAMPLES = ["jet8", "jet12", "jet20", "jet30", "jet40"]
COLORS = {
    "jet8": "#e83e8c",
    "jet12": "#2ca02c",
    "jet20": "#1da1f2",
    "jet30": "#ff6f00",
    "jet40": "#d65ad1",
}


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def ratio_error(num: float, num_err: float, den: float, den_err: float) -> float:
    if not all(math.isfinite(v) for v in (num, num_err, den, den_err)):
        return float("nan")
    if num <= 0 or den <= 0:
        return float("nan")
    r = num / den
    return r * math.hypot(num_err / num, den_err / den)


def load_all_component_rows() -> list[dict[str, str]]:
    with IN_CSV.open(newline="") as handle:
        rows = [r for r in csv.DictReader(handle) if r["mode"] == "all_components"]
    rows.sort(key=lambda r: (SAMPLES.index(r["sample"]), f(r, "bin_center")))
    return rows


def write_csv(rows: list[dict[str, object]]) -> None:
    fields = [
        "sample",
        "bin_low",
        "bin_high",
        "bin_center",
        "ppg12_source_value",
        "ppg12_source_error",
        "ppg12_fit_value",
        "ppg12_source_over_fit",
        "current_unscaled_value",
        "current_unscaled_error",
        "current_unscaled_over_ppg12",
        "common_source_exposure_scale",
        "current_common_exposure_value",
        "current_common_exposure_error",
        "current_common_exposure_over_ppg12",
        "current_common_exposure_over_fit",
        "current_common_exposure_over_ppg12_error",
        "correction_role",
    ]
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot(rows: list[dict[str, object]], exposure_scale: float) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.25,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )

    fig = plt.figure(figsize=(8.0, 8.89), dpi=200)
    top = fig.add_axes([0.13, 0.40, 0.79, 0.56])
    bot = fig.add_axes([0.13, 0.09, 0.79, 0.30], sharex=top)
    top.set_yscale("log")
    top.set_xlim(9, 50)
    top.set_ylim(5e3, 2e12)
    bot.set_ylim(0.85, 1.18)

    for sample in SAMPLES:
        pts = [r for r in rows if r["sample"] == sample]
        if not pts:
            continue
        color = COLORS[sample]
        x = [float(r["bin_center"]) for r in pts]
        top.errorbar(
            x,
            [float(r["ppg12_source_value"]) for r in pts],
            yerr=[float(r["ppg12_source_error"]) for r in pts],
            fmt="o",
            ms=4.0,
            mfc="white",
            mec=color,
            mew=1.05,
            ecolor=color,
            elinewidth=0.55,
            linestyle="none",
            zorder=3,
        )
        top.errorbar(
            x,
            [float(r["current_common_exposure_value"]) for r in pts],
            yerr=[float(r["current_common_exposure_error"]) for r in pts],
            fmt="o",
            ms=3.5,
            mfc=color,
            mec=color,
            mew=0.65,
            ecolor=color,
            elinewidth=0.5,
            linestyle="none",
            zorder=4,
        )
        bot.errorbar(
            x,
            [float(r["current_common_exposure_over_fit"]) for r in pts],
            yerr=[
                float(r["current_common_exposure_error"]) / float(r["ppg12_fit_value"])
                if float(r["ppg12_fit_value"]) > 0
                else float("nan")
                for r in pts
            ],
            fmt="o",
            ms=3.4,
            mfc=color,
            mec=color,
            mew=0.6,
            ecolor=color,
            elinewidth=0.5,
            linestyle="none",
            zorder=4,
        )
        bot.errorbar(
            x,
            [float(r["ppg12_source_over_fit"]) for r in pts],
            yerr=[
                float(r["ppg12_source_error"]) / float(r["ppg12_fit_value"])
                if float(r["ppg12_fit_value"]) > 0
                else float("nan")
                for r in pts
            ],
            fmt="o",
            ms=3.6,
            mfc="white",
            mec="black",
            mew=0.9,
            ecolor="black",
            elinewidth=0.45,
            linestyle="none",
            alpha=0.85,
            zorder=3,
        )

    fit_rows = sorted(rows, key=lambda r: float(r["bin_center"]))
    top.plot(
        [float(r["bin_center"]) for r in fit_rows],
        [float(r["ppg12_fit_value"]) for r in fit_rows],
        color="red",
        lw=1.35,
    )
    bot.axhline(1.0, color="0.4", lw=1.0, ls=(0, (5, 5)))

    top.set_ylabel("counts", fontsize=23)
    bot.set_ylabel("MC / Fit", fontsize=19)
    bot.set_xlabel(r"Leading $p_T^\mathrm{jet}$ [GeV]", fontsize=24)
    top.tick_params(which="both", labelbottom=False, labelsize=17, length=7)
    bot.tick_params(which="both", labelsize=17, length=7)
    top.minorticks_on()
    bot.minorticks_on()

    top.text(
        0.50,
        0.96,
        r"$\bf{\it{sPHENIX}}$ Internal" + "\n" + r"$p$+$p$ $\sqrt{s}=200$ GeV" + "\nPYTHIA8",
        transform=top.transAxes,
        fontsize=19,
        va="top",
    )
    top.text(
        0.50,
        0.77,
        f"current source-exposure scale = {exposure_scale:.6f}",
        transform=top.transAxes,
        fontsize=12,
        va="top",
    )

    sample_handles = [
        Line2D([0], [0], marker="o", color=COLORS[s], mfc=COLORS[s], lw=0, ms=7, label=s)
        for s in SAMPLES
    ]
    source_handles = [
        Line2D([0], [0], marker="o", color="black", mfc="white", lw=0, ms=7, label="PPG12 SDCC"),
        Line2D([0], [0], marker="o", color="black", mfc="black", lw=0, ms=7, label="Current corrected"),
        Line2D([0], [0], color="red", lw=1.5, label="PPG12 fit"),
    ]
    leg1 = top.legend(
        handles=sample_handles,
        title="sample",
        frameon=False,
        fontsize=12.5,
        title_fontsize=12.5,
        loc="lower left",
        bbox_to_anchor=(0.06, 0.03),
        handletextpad=0.45,
        labelspacing=0.42,
    )
    top.add_artist(leg1)
    top.legend(
        handles=source_handles,
        title="source",
        frameon=False,
        fontsize=12.5,
        title_fontsize=12.5,
        loc="lower left",
        bbox_to_anchor=(0.34, 0.08),
        handletextpad=0.55,
        labelspacing=0.48,
    )

    fig.savefig(OUT_PNG)
    plt.close(fig)


def main() -> None:
    source_rows = load_all_component_rows()
    nonjet8_ratios = [
        f(r, "ratio_july1_over_ppg12_source")
        for r in source_rows
        if r["sample"] != "jet8" and math.isfinite(f(r, "ratio_july1_over_ppg12_source"))
    ]
    if not nonjet8_ratios:
        raise RuntimeError("no non-jet8 ratios available")
    exposure_scale = 1.0 / statistics.median(nonjet8_ratios)

    rows: list[dict[str, object]] = []
    by_sample_ratio: dict[str, list[float]] = {s: [] for s in SAMPLES}
    for r in source_rows:
        sample = r["sample"]
        ppg = f(r, "ppg12_source_value")
        ppg_err = f(r, "ppg12_source_error")
        cur = f(r, "july1_value")
        cur_err = f(r, "july1_error")
        fit = f(r, "ppg12_source_value") / f(r, "ratio_july1_over_ppg12_source")
        # The audit CSV does not carry the fit. Recover it from the SDCC source
        # table by reading the companion strict points file.
        rows.append(
            {
                "sample": sample,
                "bin_low": f(r, "bin_low"),
                "bin_high": f(r, "bin_high"),
                "bin_center": f(r, "bin_center"),
                "ppg12_source_value": ppg,
                "ppg12_source_error": ppg_err,
                "ppg12_fit_value": 0.0,
                "ppg12_source_over_fit": 0.0,
                "current_unscaled_value": cur,
                "current_unscaled_error": cur_err,
                "current_unscaled_over_ppg12": cur / ppg if ppg > 0 else float("nan"),
                "common_source_exposure_scale": exposure_scale,
                "current_common_exposure_value": cur * exposure_scale,
                "current_common_exposure_error": cur_err * exposure_scale,
                "current_common_exposure_over_ppg12": cur * exposure_scale / ppg if ppg > 0 else float("nan"),
                "current_common_exposure_over_fit": 0.0,
                "current_common_exposure_over_ppg12_error": ratio_error(cur * exposure_scale, cur_err * exposure_scale, ppg, ppg_err),
                "correction_role": "jet8_residual_check" if sample == "jet8" else "common_nonjet8_source_exposure",
            }
        )

    source_by_key: dict[tuple[str, float], dict[str, str]] = {}
    source_csv = BASE / "reference_audit/sdcc_pull/live_ppg12_ian_source_points.csv"
    with source_csv.open(newline="") as handle:
        for srow in csv.DictReader(handle):
            if srow["group"] == "jet":
                source_by_key[(srow["sample"], round(float(srow["bin_center"]), 6))] = srow
    for row in rows:
        key = (str(row["sample"]), round(float(row["bin_center"]), 6))
        srow = source_by_key.get(key)
        if not srow:
            raise RuntimeError(f"missing source fit row for {key}")
        fit_value = float(srow["fit_value"])
        row["ppg12_fit_value"] = fit_value
        row["ppg12_source_over_fit"] = float(srow["fit_ratio"])
        row["current_common_exposure_over_fit"] = (
            float(row["current_common_exposure_value"]) / fit_value if fit_value > 0 else float("nan")
        )
        by_sample_ratio[str(row["sample"])].append(float(row["current_common_exposure_over_ppg12"]))

    write_csv(rows)
    plot(rows, exposure_scale)

    summary = {}
    for sample, vals in by_sample_ratio.items():
        vals = [v for v in vals if math.isfinite(v)]
        if vals:
            summary[sample] = {
                "n": len(vals),
                "median_current_corrected_over_ppg12": statistics.median(vals),
                "min": min(vals),
                "max": max(vals),
            }

    OUT_MANIFEST.write_text(
        json.dumps(
            {
                "status": "ok_common_source_exposure_correction_with_jet8_residual",
                "plot_png": str(OUT_PNG),
                "points_csv": str(OUT_CSV),
                "input_component_audit_csv": str(IN_CSV),
                "ppg12_source_csv": str(source_csv),
                "common_source_exposure_scale": exposure_scale,
                "scale_definition": "1 / median(July1 all-components / PPG12 source) over jet12-40 bins only",
                "why_not_arbitrary_plot_scale": (
                    "PPG12 Fig.6 is a weighted-count source diagnostic, so raw counts scale with effective processed event population. "
                    "This factor converts the July1 full period/component population onto the PPG12 no-suffix source-count convention. "
                    "It is written explicitly and jet8 is excluded from deriving it so the jet8-only residual remains visible."
                ),
                "summary_by_sample": summary,
                "jet8_residual_factor_vs_nonjet8": summary.get("jet8", {}).get("median_current_corrected_over_ppg12"),
                "needs_condor_rerun": False,
                "rerun_note": (
                    "No rerun is needed to remove the common jet12-40 offset; that is an offline source-count convention. "
                    "A jet8-only rerun is only justified if downstream analyses require absolute Fig.6 raw-count closure in 9-14 GeV against the older PPG12 no-suffix source."
                ),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(OUT_PNG)
    print(OUT_CSV)
    print(OUT_MANIFEST)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Render split inclusive-jet Fig.6 reference/canonical fit comparison.

Left side: PPG12 IAN/source-scope reference, restricted to jet12-40, with the
IAN PPG12 fit guide.  This avoids treating the saved PPG12 jet8 object as the
authoritative low-pT stitch.

Right side: corrected THE-94 output over the full jet8-40 range, with a smooth
fit guide fitted to the THE-94 points themselves.  The THE-94 points are not
rescaled; this is the canonical RecoilJets jet8 policy view.
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
CURRENT_INPUT_CSV = (
    THE94_BASE
    / "no_scale_overlay/inclusivejet_fig6_true_period_combined_no_scale_data_over_fit_points.csv"
)
PPG12_IAN_INPUT_CSV = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024/"
    "strict_stitched_inclusivejet/jet_sdcc_over_current_shape_overlay_ppg12_jet8_xsecfix_points.csv"
)
OUT_DIR = THE94_BASE / "canonical_jet8_policy_overlay"
OUT_PNG = OUT_DIR / "inclusivejet_fig6_split_ppg12ian_jet12to40_the94_jet8to40_smoothfit.png"
OUT_POINTS = OUT_DIR / "inclusivejet_fig6_split_ppg12ian_jet12to40_the94_jet8to40_smoothfit_points.csv"
OUT_MANIFEST = OUT_DIR / "inclusivejet_fig6_split_ppg12ian_jet12to40_the94_jet8to40_smoothfit_manifest.json"
OUT_NOTE = OUT_DIR / "inclusivejet_fig6_split_ppg12ian_jet12to40_the94_jet8to40_smoothfit_note.txt"

SAMPLES = ["jet8", "jet12", "jet20", "jet30", "jet40"]
COLORS = {
    "jet8": "#d62aa0",
    "jet12": "#239b35",
    "jet20": "#169ce8",
    "jet30": "#ff6f00",
    "jet40": "#c12ac7",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def ff(row: dict[str, str], key: str) -> float:
    return float(row[key])


def fit_log_poly(rows: list[dict[str, str]], value_key: str, degree: int) -> np.ndarray:
    x = np.array([ff(r, "bin_center") for r in rows], dtype=float)
    y = np.array([ff(r, value_key) for r in rows], dtype=float)
    mask = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    return np.polyfit(np.log(x[mask]), np.log(y[mask]), degree)


def eval_log_poly(coeff: np.ndarray, x: np.ndarray | float) -> np.ndarray | float:
    return np.exp(np.polyval(coeff, np.log(x)))


def ppg12_fit_at(row: dict[str, str]) -> float:
    return ff(row, "ppg12_fit_value")


def build_output_rows(
    current_rows: list[dict[str, str]],
    ppg12_ian_rows: list[dict[str, str]],
    the94_coeff: np.ndarray,
) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    ppg12_by_key = {
        (r["sample"], float(r["bin_center"])): r
        for r in ppg12_ian_rows
    }
    for r in current_rows:
        x = ff(r, "bin_center")
        ppg12_fit = ppg12_fit_at(r)
        the94_fit = float(eval_log_poly(the94_coeff, x))
        ppg12_ian = ppg12_by_key.get((r["sample"], x))
        out.append(
            {
                "sample": r["sample"],
                "bin_low": ff(r, "bin_low"),
                "bin_high": ff(r, "bin_high"),
                "bin_center": x,
                "ppg12_ian_value": ff(ppg12_ian, "ppg12_sdcc_value") if ppg12_ian else float("nan"),
                "ppg12_ian_error": ff(ppg12_ian, "ppg12_sdcc_error") if ppg12_ian else float("nan"),
                "ppg12_ian_fit_value": ff(ppg12_ian, "ppg12_fit_value") if ppg12_ian else float("nan"),
                "ppg12_ian_over_fit": ff(ppg12_ian, "ppg12_sdcc_over_fit") if ppg12_ian else float("nan"),
                "ppg12_ian_over_fit_error": (
                    ff(ppg12_ian, "ppg12_sdcc_error") / ff(ppg12_ian, "ppg12_fit_value")
                    if ppg12_ian
                    else float("nan")
                ),
                "ppg12_true_period_value": ff(r, "ppg12_true_period_value"),
                "ppg12_true_period_error": ff(r, "ppg12_true_period_error"),
                "ppg12_true_period_fit_value": ppg12_fit,
                "ppg12_true_period_over_fit": ff(r, "ppg12_true_period_value") / ppg12_fit,
                "the94_canonical_value": ff(r, "current_kept_value"),
                "the94_canonical_error": ff(r, "current_kept_error"),
                "the94_smoothfit_value": the94_fit,
                "the94_canonical_over_the94_smoothfit": ff(r, "current_kept_value") / the94_fit,
                "the94_canonical_over_the94_smoothfit_error": ff(r, "current_kept_error") / the94_fit,
                "the94_current_over_ppg12_true_period_reference": ff(r, "current_over_ppg12_true_period"),
            }
        )
    return out


def sample_rows(rows: list[dict[str, object]] | list[dict[str, str]], sample: str) -> list:
    return [r for r in rows if r["sample"] == sample]


def ratio_stats(rows: list[dict[str, object]], key: str) -> dict[str, float]:
    vals = [float(r[key]) for r in rows]
    return {
        "min": float(min(vals)),
        "max": float(max(vals)),
        "median": float(np.median(vals)),
        "rms_about_1": float(math.sqrt(sum((v - 1.0) ** 2 for v in vals) / len(vals))),
    }


def write_points(rows: list[dict[str, object]]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with OUT_POINTS.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def marker_handles(samples: list[str]) -> list[Line2D]:
    return [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=COLORS[s],
            markeredgecolor=COLORS[s],
            markersize=7,
            label=s,
        )
        for s in samples
    ]


def draw_points(ax, rows, value_key, error_key, samples, *, open_marker: bool) -> None:
    for sample in samples:
        rs = [r for r in rows if r["sample"] == sample]
        if not rs:
            continue
        x = [float(r["bin_center"]) for r in rs]
        y = [float(r[value_key]) for r in rs]
        e = [float(r[error_key]) for r in rs]
        ax.errorbar(
            x,
            y,
            yerr=e,
            fmt="o",
            ms=4.8,
            mfc="none" if open_marker else COLORS[sample],
            mec=COLORS[sample],
            mew=1.2 if open_marker else 0.7,
            ecolor=COLORS[sample],
            elinewidth=0.65,
            linestyle="none",
            zorder=3,
        )


def render(rows: list[dict[str, object]], the94_coeff: np.ndarray) -> None:
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
            "xtick.major.size": 6.5,
            "ytick.major.size": 6.5,
            "xtick.minor.size": 3.5,
            "ytick.minor.size": 3.5,
        }
    )

    left_samples = ["jet12", "jet20", "jet30", "jet40"]
    right_samples = SAMPLES
    left_rows = [r for r in rows if r["sample"] in left_samples]
    right_rows = rows

    fig = plt.figure(figsize=(15.6, 8.4), dpi=180)
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[3.1, 1.25],
        hspace=0.035,
        wspace=0.16,
        left=0.075,
        right=0.985,
        bottom=0.105,
        top=0.94,
    )
    ax_l = fig.add_subplot(gs[0, 0])
    ax_lr = fig.add_subplot(gs[1, 0], sharex=ax_l)
    ax_r = fig.add_subplot(gs[0, 1])
    ax_rr = fig.add_subplot(gs[1, 1], sharex=ax_r)

    for ax in [ax_l, ax_r]:
        ax.set_yscale("log")
        ax.tick_params(labelbottom=False)
    ax_l.set_xlim(13.6, 50.4)
    ax_r.set_xlim(8.6, 50.4)
    ax_l.set_ylim(7.0e5, 1.4e12)
    ax_r.set_ylim(7.0e5, 6.0e12)
    ax_lr.set_ylim(0.84, 1.16)
    ax_rr.set_ylim(0.94, 1.04)

    # Left: PPG12 IAN/source stable reference region.
    draw_points(ax_l, left_rows, "ppg12_ian_value", "ppg12_ian_error", left_samples, open_marker=True)
    x_fit_l = np.linspace(14.0, 50.0, 500)
    ppg12_fit_coeff = fit_log_poly(
        [
            {
                "bin_center": str(r["bin_center"]),
                "ppg12_ian_fit_value": str(r["ppg12_ian_fit_value"]),
            }
            for r in rows
            if np.isfinite(float(r["ppg12_ian_fit_value"]))
        ],
        "ppg12_ian_fit_value",
        5,
    )
    ax_l.plot(x_fit_l, eval_log_poly(ppg12_fit_coeff, x_fit_l), color="red", lw=1.8, label="PPG12 fit")
    draw_points(
        ax_lr,
        left_rows,
        "ppg12_ian_over_fit",
        "ppg12_ian_over_fit_error",
        left_samples,
        open_marker=False,
    )

    # Right: THE-94 full canonical range.
    draw_points(ax_r, right_rows, "the94_canonical_value", "the94_canonical_error", right_samples, open_marker=False)
    x_fit_r = np.linspace(9.0, 50.0, 700)
    ax_r.plot(x_fit_r, eval_log_poly(the94_coeff, x_fit_r), color="red", lw=1.8, label="THE-94 smooth fit")
    draw_points(
        ax_rr,
        right_rows,
        "the94_canonical_over_the94_smoothfit",
        "the94_canonical_over_the94_smoothfit_error",
        right_samples,
        open_marker=False,
    )

    for ax in [ax_lr, ax_rr]:
        ax.axhline(1.0, color="0.45", lw=1.0, ls=(0, (4, 4)))
        ax.set_ylabel("Data / Fit", fontsize=18)
        ax.tick_params(labelsize=13)
        ax.minorticks_on()
    for ax in [ax_l, ax_r]:
        ax.set_ylabel("counts", fontsize=21)
        ax.tick_params(labelsize=14)
        ax.minorticks_on()
    ax_lr.set_xlabel(r"Leading $p_T^{\mathrm{jet}}$ [GeV]", fontsize=20)
    ax_rr.set_xlabel(r"Leading $p_T^{\mathrm{jet}}$ [GeV]", fontsize=20)

    ax_l.text(0.05, 0.92, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax_l.transAxes, fontsize=18)
    ax_l.text(0.05, 0.85, r"$p$+$p$ $\sqrt{s}=200$ GeV", transform=ax_l.transAxes, fontsize=16)
    ax_l.text(0.05, 0.78, "PYTHIA8", transform=ax_l.transAxes, fontsize=16)
    ax_l.text(0.05, 0.12, "PPG12 IAN reference only\njet12-40", transform=ax_l.transAxes, fontsize=15)

    ax_r.text(0.46, 0.92, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax_r.transAxes, fontsize=18)
    ax_r.text(0.46, 0.85, r"$p$+$p$ $\sqrt{s}=200$ GeV", transform=ax_r.transAxes, fontsize=16)
    ax_r.text(0.46, 0.78, "PYTHIA8", transform=ax_r.transAxes, fontsize=16)
    ax_r.text(
        0.05,
        0.12,
        "THE-94 canonical\njet8-40, no hidden jet8 factor",
        transform=ax_r.transAxes,
        fontsize=15,
    )

    left_legend = marker_handles(left_samples)
    right_legend = marker_handles(right_samples)
    ax_l.legend(handles=left_legend + [Line2D([0], [0], color="red", lw=1.8, label="PPG12 fit")],
                loc="upper right", frameon=False, fontsize=13, handlelength=1.2)
    ax_r.legend(handles=right_legend + [Line2D([0], [0], color="red", lw=1.8, label="THE-94 fit")],
                loc="upper right", frameon=False, fontsize=13, handlelength=1.2)

    fig.savefig(OUT_PNG)
    plt.close(fig)


def main() -> None:
    current_rows = read_rows(CURRENT_INPUT_CSV)
    ppg12_ian_rows = read_rows(PPG12_IAN_INPUT_CSV)
    the94_coeff = fit_log_poly(current_rows, "current_kept_value", 5)
    rows = build_output_rows(current_rows, ppg12_ian_rows, the94_coeff)
    write_points(rows)
    render(rows, the94_coeff)

    left_rows = [r for r in rows if r["sample"] in {"jet12", "jet20", "jet30", "jet40"}]
    right_rows = rows
    jet8_rows = [r for r in rows if r["sample"] == "jet8"]
    jet12_rows = [r for r in rows if r["sample"] == "jet12"]
    last_jet8 = max(jet8_rows, key=lambda r: float(r["bin_center"]))
    first_jet12 = min(jet12_rows, key=lambda r: float(r["bin_center"]))
    manifest = {
        "status": "ok_split_ppg12_the94_fit_comparison",
        "current_input_points_csv": str(CURRENT_INPUT_CSV),
        "ppg12_ian_input_points_csv": str(PPG12_IAN_INPUT_CSV),
        "png": str(OUT_PNG),
        "points_csv": str(OUT_POINTS),
        "interpretation": (
            "Left panel: PPG12 IAN/source reference restricted to jet12-40 with "
            "the IAN PPG12 fit guide. Right panel: THE-94 canonical jet8-40 output with "
            "a smooth fit guide fitted to THE-94 itself. No point rescaling."
        ),
        "the94_point_scale_factor": 1.0,
        "the94_smooth_fit": {
            "type": "degree5_loglog_polynomial",
            "coefficients_high_to_low": [float(x) for x in the94_coeff],
            "fit_input": "all THE-94 current_kept_value bins, jet8 through jet40",
            "ratio_stats_all_bins": ratio_stats(right_rows, "the94_canonical_over_the94_smoothfit"),
        },
        "ppg12_left_panel": {
            "samples": ["jet12", "jet20", "jet30", "jet40"],
            "ratio_stats": ratio_stats(left_rows, "ppg12_ian_over_fit"),
            "jet8_omitted_reason": "saved PPG12 jet8 is documented as a historical source/normalization exception",
        },
        "boundary_summary": {
            "THE94_last_jet8_bin": [last_jet8["bin_low"], last_jet8["bin_high"]],
            "THE94_first_jet12_bin": [first_jet12["bin_low"], first_jet12["bin_high"]],
            "THE94_last_jet8_over_THE94_fit": last_jet8["the94_canonical_over_the94_smoothfit"],
            "THE94_first_jet12_over_THE94_fit": first_jet12["the94_canonical_over_the94_smoothfit"],
            "THE94_boundary_ratio_last_jet8_to_first_jet12": float(
                last_jet8["the94_canonical_over_the94_smoothfit"]
            )
            / float(first_jet12["the94_canonical_over_the94_smoothfit"]),
        },
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    OUT_NOTE.write_text(
        "\n".join(
            [
                "Split inclusive-jet Fig.6 comparison.",
                "",
                "Left: PPG12 IAN/source reference restricted to jet12-40 with the IAN PPG12 fit guide.",
                "Right: THE-94 canonical output from jet8-40 with a smooth fit guide fitted to THE-94 points.",
                "No THE-94 points are rescaled. The saved PPG12 jet8 object is not shown as a reference target here.",
                "",
                json.dumps(manifest["boundary_summary"], indent=2, sort_keys=True),
                "",
            ]
        )
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

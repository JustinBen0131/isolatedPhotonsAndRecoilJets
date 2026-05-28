#!/usr/bin/env python3
from __future__ import annotations

import json
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppPhotonMLPipeline/"
    "ppg12_basev3E_currentIAN_finalShuhang_20260526_2115/validation/fullsim_shuhang_overlay"
)
DEFAULT_SUMMARY = (
    ROOT
    / "pp_currentian_basev3e_bdt_score_overlay_vs_ppg12_pp_noCent_bdt_22_28_noNPB_eta0-pt3-cut0_summary.json"
)
DEFAULT_DATA_PROJECTION = ROOT / "ppg12_data_h2d_bdt_eta0_pt3_cut0_projection.json"
DEFAULT_OUT = ROOT / "pp_currentian_basev3e_bdt_fourcurve_tall_overlay_fig13_noNPB.png"


def load_json_with_optional_log_prefix(path: Path) -> dict:
    text = path.read_text()
    if not text.lstrip().startswith("{"):
        text = text[text.index("{") :]
    return json.loads(text)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Make the tall four-curve pp baseV3E BDT overlay used for the Fig. 13-style slide RHS."
    )
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--data-projection", type=Path, default=DEFAULT_DATA_PROJECTION)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    data = load_json_with_optional_log_prefix(args.summary)
    data_projection = load_json_with_optional_log_prefix(args.data_projection)
    bins = np.asarray(data["bins"], dtype=float)
    centers = 0.5 * (bins[:-1] + bins[1:])

    curves = {
        "This analysis signal": np.asarray(data["this_analysis_signal_hist"], dtype=float),
        "This analysis inclusive": np.asarray(data["this_analysis_inclusive_hist"], dtype=float),
        "PPG12 signal": np.asarray(data["shuhang_signal_hist"], dtype=float),
        "PPG12 inclusive": np.asarray(data["shuhang_inclusive_hist"], dtype=float),
    }
    ppg12_data = np.asarray(data_projection["hist"], dtype=float)

    fig = plt.figure(figsize=(8.2, 10.8), dpi=190)
    gs = fig.add_gridspec(
        2,
        1,
        height_ratios=[3.7, 1.05],
        hspace=0.05,
        left=0.16,
        right=0.96,
        top=0.95,
        bottom=0.11,
    )
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)

    styles = {
        "This analysis signal": dict(color="#e41a1c", lw=2.2, ls="-"),
        "This analysis inclusive": dict(color="#1f4fff", lw=2.2, ls="-"),
        "PPG12 signal": dict(color="#8b1a1a", lw=2.0, ls="--"),
        "PPG12 inclusive": dict(color="#173a9a", lw=2.0, ls="--"),
    }
    for label, hist in curves.items():
        ax.step(centers, hist, where="mid", label=label, **styles[label])

    peak_specs = [
        ("This analysis inclusive", "This analysis incl. peak", 0.27),
        ("PPG12 inclusive", "PPG12 incl. peak", 0.26),
    ]
    for curve_label, text, label_x in peak_specs:
        hist = curves[curve_label]
        peak_idx = int(np.nanargmax(hist))
        peak_xy = (float(centers[peak_idx]), float(hist[peak_idx]))
        color = styles[curve_label]["color"]
        ax.hlines(
            peak_xy[1],
            xmin=0.0,
            xmax=0.50 if curve_label == "This analysis inclusive" else 0.62,
            colors=color,
            linestyles=(0, (5, 3)),
            linewidth=1.8,
            alpha=0.75,
            zorder=4,
        )
        ax.plot(
            [peak_xy[0]],
            [peak_xy[1]],
            marker="o",
            ms=8.5,
            mfc="white",
            mec=color,
            mew=2.4,
            linestyle="none",
            zorder=8,
        )
        ax.text(
            label_x,
            peak_xy[1] - (0.018 if curve_label == "This analysis inclusive" else -0.008),
            text,
            ha="left",
            va="top" if curve_label == "This analysis inclusive" else "bottom",
            fontsize=12.7,
            color=color,
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec=color, lw=1.0, alpha=0.92),
            zorder=9,
        )

    ax.text(
        0.58,
        0.70,
        r"$\bf{\it{sPHENIX}}$ Internal" + "\n"
        r"$p$+$p$ $\sqrt{s}=200$ GeV" + "\n"
        r"$|\eta| < 0.7$" + "\n"
        r"$22 < E_T < 28$ GeV" + "\n"
        r"w/o NPB cut",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=17,
        linespacing=1.15,
    )
    ax.legend(loc="upper right", bbox_to_anchor=(0.985, 0.985), frameon=False, fontsize=14.2, handlelength=2.2)
    ax.set_ylabel("unit-normalized counts", fontsize=23)
    ax.set_ylim(0.0, 0.62)
    ax.set_xlim(0.0, 1.0)
    ax.tick_params(axis="both", which="major", direction="in", top=True, right=True, length=8, labelsize=19)
    ax.tick_params(axis="both", which="minor", direction="in", top=True, right=True, length=4)
    ax.minorticks_on()
    ax.tick_params(labelbottom=False)

    diff = ppg12_data - curves["This analysis inclusive"]
    rax.axhline(0.0, color="black", lw=1.1, ls="--")
    rax.plot(centers, diff, "o", ms=4.6, color="black", label="Data - this analysis Inc. M.C.")
    rax.set_ylabel("Data - Inc. M.C.", fontsize=20)
    rax.set_xlabel("BDT score", fontsize=24, ha="right", x=1.0)
    finite_diff = diff[np.isfinite(diff)]
    if finite_diff.size:
        margin = max(0.02, 0.12 * float(np.max(finite_diff) - np.min(finite_diff)))
        rax.set_ylim(float(np.min(finite_diff)) - margin, float(np.max(finite_diff)) + margin)
    else:
        rax.set_ylim(-0.18, 0.55)
    rax.tick_params(axis="both", which="major", direction="in", top=True, right=True, length=8, labelsize=18)
    rax.tick_params(axis="both", which="minor", direction="in", top=True, right=True, length=4)
    rax.minorticks_on()
    rax.legend(loc="upper right", frameon=False, fontsize=11)

    out = args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    print(out)


if __name__ == "__main__":
    main()

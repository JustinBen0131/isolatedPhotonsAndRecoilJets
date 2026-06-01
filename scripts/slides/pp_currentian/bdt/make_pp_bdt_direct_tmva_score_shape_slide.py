#!/usr/bin/env python3
"""Build a slide-ready pp BDT score-shape comparison from same-row direct scores."""

from __future__ import annotations
# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath
_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
_CODEX_SCRIPTS_DIR_STR = str(_CODEX_SCRIPTS_DIR)
if _CODEX_SCRIPTS_DIR_STR not in _codex_sys.path:
    _codex_sys.path.append(_CODEX_SCRIPTS_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_SCRIPTS_DIR_STR

import argparse
import json
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch


INK = "#111827"
MUTED = "#4B5563"
GRID = "#D1D5DB"
BLUE = "#2563EB"
RED = "#DC2626"
PANEL_EDGE = "#CBD5E1"
SOFT_BLUE = "#EFF6FF"
SOFT_GREEN = "#ECFDF5"
SOFT_AMBER = "#FFF7ED"


def set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": ["Times New Roman", "Times", "DejaVu Serif"],
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": INK,
            "axes.linewidth": 1.35,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "mathtext.fontset": "dejavuserif",
        }
    )


def rounded_box(fig: plt.Figure, x: float, y: float, w: float, h: float, fc: str, ec: str = PANEL_EDGE) -> None:
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            transform=fig.transFigure,
            boxstyle="round,pad=0.012,rounding_size=0.012",
            linewidth=1.1,
            edgecolor=ec,
            facecolor=fc,
            zorder=-1,
        )
    )


def fig_text(
    fig: plt.Figure,
    x: float,
    y: float,
    text: str,
    *,
    size: float,
    color: str = INK,
    weight: str = "normal",
    linespacing: float = 1.18,
) -> None:
    fig.text(x, y, text, ha="left", va="top", fontsize=size, color=color, fontweight=weight, linespacing=linespacing)


def wrapped_text(
    fig: plt.Figure,
    x: float,
    y: float,
    text: str,
    *,
    width: int,
    size: float,
    color: str = INK,
    linespacing: float = 1.18,
) -> None:
    wrapped = "\n".join(textwrap.fill(part, width=width) for part in text.split("\n"))
    fig_text(fig, x, y, wrapped, size=size, color=color, linespacing=linespacing)


def natural_sample_key(sample: str) -> tuple[str, int]:
    head = "".join(ch for ch in sample if not ch.isdigit())
    digits = "".join(ch for ch in sample if ch.isdigit())
    return head, int(digits or 0)


def unit_hist(values: np.ndarray, bins: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    counts, _ = np.histogram(values, bins=bins)
    total = counts.sum()
    if total <= 0:
        return counts.astype(float)
    return counts.astype(float) / float(total)


def dist_stats(values: np.ndarray, bins: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    hist = unit_hist(values, bins)
    peak_idx = int(np.argmax(hist)) if hist.size else 0
    center = 0.5 * (bins[peak_idx] + bins[peak_idx + 1]) if hist.size else float("nan")
    return {
        "rows": int(values.size),
        "mean": float(np.mean(values)) if values.size else float("nan"),
        "median": float(np.median(values)) if values.size else float("nan"),
        "peak_score_bin_center": float(center),
        "peak_bin_fraction": float(hist[peak_idx]) if hist.size else float("nan"),
        "frac_score_lt_0p1": float(np.mean(values < 0.1)) if values.size else float("nan"),
        "frac_score_gt_0p8": float(np.mean(values > 0.8)) if values.size else float("nan"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--summary-json", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    set_style()
    args.outdir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.csv)
    required = {"class", "sample", "cluster_Et", "our_score", "shuhang_split_score"}
    missing = required.difference(frame.columns)
    if missing:
        raise SystemExit(f"Missing required columns in {args.csv}: {sorted(missing)}")

    frame = frame[frame["class"].isin(["signal", "inclusive"])].copy()
    frame["our_score"] = np.clip(frame["our_score"].astype(float), 0.0, 1.0)
    frame["shuhang_split_score"] = np.clip(frame["shuhang_split_score"].astype(float), 0.0, 1.0)

    bins = np.linspace(0.0, 1.0, 51)
    signal = frame[frame["class"] == "signal"]
    inclusive = frame[frame["class"] == "inclusive"]

    hists = {
        "this_signal": unit_hist(signal["our_score"].to_numpy(), bins),
        "ppg12_signal": unit_hist(signal["shuhang_split_score"].to_numpy(), bins),
        "this_inclusive": unit_hist(inclusive["our_score"].to_numpy(), bins),
        "ppg12_inclusive": unit_hist(inclusive["shuhang_split_score"].to_numpy(), bins),
    }
    stats = {
        "this_analysis_signal": dist_stats(signal["our_score"].to_numpy(), bins),
        "ppg12_direct_tmva_signal": dist_stats(signal["shuhang_split_score"].to_numpy(), bins),
        "this_analysis_inclusive": dist_stats(inclusive["our_score"].to_numpy(), bins),
        "ppg12_direct_tmva_inclusive": dist_stats(inclusive["shuhang_split_score"].to_numpy(), bins),
    }

    signal_samples = sorted(signal["sample"].astype(str).unique(), key=natural_sample_key)
    inclusive_samples = sorted(inclusive["sample"].astype(str).unique(), key=natural_sample_key)

    out_png = args.outdir / "pp_basev3e_direct_tmva_score_shape_slide20_candidate.png"
    plot_png = args.outdir / "pp_basev3e_direct_tmva_score_shape_plotonly.png"
    out_json = args.outdir / "pp_basev3e_direct_tmva_score_shape_slide20_candidate.json"

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.subplots_adjust(0, 0, 1, 1)

    fig_text(fig, 0.045, 0.945, "pp baseV3E Score Shapes: Direct Model Application", size=29.0, weight="bold")
    fig_text(
        fig,
        0.045,
        0.902,
        "Same candidate rows scored with this-analysis XGBoost and Shuhang/PPG12 split-TMVA baseV3E model outputs.",
        size=15.3,
        color=MUTED,
    )

    ax = fig.add_axes([0.060, 0.145, 0.570, 0.695])
    ax.stairs(hists["this_signal"], bins, color=RED, lw=3.0, label="This analysis signal")
    ax.stairs(hists["ppg12_signal"], bins, color=RED, lw=3.0, ls="--", label="PPG12 direct-TMVA signal")
    ax.stairs(hists["this_inclusive"], bins, color=BLUE, lw=3.0, label="This analysis inclusive")
    ax.stairs(hists["ppg12_inclusive"], bins, color=BLUE, lw=3.0, ls="--", label="PPG12 direct-TMVA inclusive")
    ymax = max(np.max(v) for v in hists.values()) * 1.25
    this_signal_peak = stats["this_analysis_signal"]["peak_score_bin_center"]
    ppg12_signal_peak = stats["ppg12_direct_tmva_signal"]["peak_score_bin_center"]
    this_inclusive_peak = stats["this_analysis_inclusive"]["peak_score_bin_center"]
    ppg12_inclusive_peak = stats["ppg12_direct_tmva_inclusive"]["peak_score_bin_center"]
    this_inclusive_peak_height = stats["this_analysis_inclusive"]["peak_bin_fraction"]
    ppg12_inclusive_peak_height = stats["ppg12_direct_tmva_inclusive"]["peak_bin_fraction"]
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, max(0.06, ymax))
    ax.set_xlabel("BDT score", fontsize=21)
    ax.set_ylabel("Unit-normalized / 0.02 score bin", fontsize=18.3)
    ax.grid(color=GRID, lw=0.9, alpha=0.72)
    ax.tick_params(axis="both", labelsize=15, length=7)
    ax.legend(
        loc="upper right",
        bbox_to_anchor=(0.985, 0.988),
        ncol=1,
        fontsize=12.7,
        frameon=True,
        facecolor="white",
        edgecolor="#D1D5DB",
        framealpha=0.96,
    )
    x0, x1 = 0.035, 0.285
    ax.hlines(
        this_inclusive_peak_height,
        x0,
        x1,
        color=BLUE,
        lw=3.0,
        alpha=0.78,
        zorder=4,
    )
    ax.hlines(
        ppg12_inclusive_peak_height,
        x0,
        x1,
        color=BLUE,
        lw=3.0,
        ls=(0, (7, 4)),
        alpha=0.78,
        zorder=4,
    )
    label_x = x1 + 0.012
    ax.text(
        label_x,
        this_inclusive_peak_height,
        f"This analysis inclusive peak height: {this_inclusive_peak_height:.2f}",
        ha="left",
        va="center",
        fontsize=11.4,
        color=BLUE,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 2.2},
    )
    ax.text(
        label_x,
        ppg12_inclusive_peak_height,
        f"PPG12 inclusive peak height: {ppg12_inclusive_peak_height:.2f}",
        ha="left",
        va="center",
        fontsize=11.4,
        color=BLUE,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 2.2},
    )
    ax.text(
        0.965,
        0.742,
        r"$\bf{\it{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=16.5,
        color=INK,
    )
    ax.text(
        0.965,
        0.680,
        r"$p{+}p$ $\sqrt{s}=200$ GeV" "\n" r"$|\eta|<0.7$, $22<E_T<28$ GeV",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=13.4,
        color=INK,
        linespacing=1.15,
    )

    rhs_box_x = 0.665
    rhs_x = 0.686
    rhs_w = 0.295

    rounded_box(fig, rhs_box_x, 0.690, rhs_w, 0.165, SOFT_BLUE, "#BFDBFE")
    fig_text(fig, rhs_x, 0.825, "Direct-model comparison", size=19.2, color=BLUE, weight="bold")
    wrapped_text(
        fig,
        rhs_x,
        0.780,
        "PPG12 curves here are Shuhang's split-TMVA model applied directly to the same rows, not pre-made ROOT BDT histograms.",
        width=43,
        size=12.9,
        linespacing=1.22,
    )

    rounded_box(fig, rhs_box_x, 0.465, rhs_w, 0.175, SOFT_GREEN, "#A7F3D0")
    fig_text(fig, rhs_x, 0.612, "What this shows", size=19.2, color="#047857", weight="bold")
    this_inc_low = stats["this_analysis_inclusive"]["frac_score_lt_0p1"]
    ppg12_inc_low = stats["ppg12_direct_tmva_inclusive"]["frac_score_lt_0p1"]
    wrapped_text(
        fig,
        rhs_x,
        0.566,
        f"Both direct models put a large inclusive fraction below 0.1: this analysis {this_inc_low:.1%}, PPG12 direct {ppg12_inc_low:.1%}. The low-score spike follows the row sample.",
        width=43,
        size=12.9,
        linespacing=1.23,
    )

    rounded_box(fig, rhs_box_x, 0.250, rhs_w, 0.165, "white", PANEL_EDGE)
    fig_text(fig, rhs_x, 0.387, "Input rows", size=18.5, color=INK, weight="bold")
    input_text = (
        f"{len(frame):,} same rows\n"
        f"{len(signal):,} signal, {len(inclusive):,} inclusive\n"
        f"Signal: {', '.join(signal_samples)}\n"
        f"Inclusive: {', '.join(inclusive_samples)}"
    )
    wrapped_text(fig, rhs_x, 0.343, input_text, width=45, size=12.2, linespacing=1.23)

    rounded_box(fig, rhs_box_x, 0.070, rhs_w, 0.130, SOFT_AMBER, "#FDBA74")
    fig.text(rhs_x, 0.168, "Most likely signal-peak source:", ha="left", va="top", fontsize=13.0, fontweight="bold", color=INK)
    fig.text(
        rhs_x,
        0.134,
        "same rows + direct TMVA scoring\n"
        "isolate this to the learned model.\n"
        "The 0.93 to 0.95 shift is likely\n"
        "score calibration/training sample.",
        ha="left",
        va="top",
        fontsize=11.7,
        color=INK,
        linespacing=1.16,
    )

    fig.savefig(out_png, dpi=160)
    plt.close(fig)

    plot_fig, plot_ax = plt.subplots(figsize=(9.8, 7.0), dpi=180)
    plot_ax.stairs(hists["this_signal"], bins, color=RED, lw=2.8, label="This analysis signal")
    plot_ax.stairs(hists["ppg12_signal"], bins, color=RED, lw=2.8, ls="--", label="PPG12 direct-TMVA signal")
    plot_ax.stairs(hists["this_inclusive"], bins, color=BLUE, lw=2.8, label="This analysis inclusive")
    plot_ax.stairs(hists["ppg12_inclusive"], bins, color=BLUE, lw=2.8, ls="--", label="PPG12 direct-TMVA inclusive")
    plot_ax.hlines(this_inclusive_peak_height, 0.035, 0.285, color=BLUE, lw=2.6, alpha=0.78, zorder=4)
    plot_ax.hlines(ppg12_inclusive_peak_height, 0.035, 0.285, color=BLUE, lw=2.6, ls=(0, (7, 4)), alpha=0.78, zorder=4)
    plot_ax.set_xlim(0.0, 1.0)
    plot_ax.set_ylim(0.0, max(0.06, ymax))
    plot_ax.set_xlabel("BDT score", fontsize=18)
    plot_ax.set_ylabel("Unit-normalized / 0.02 score bin", fontsize=16)
    plot_ax.grid(color=GRID, lw=0.8, alpha=0.72)
    plot_ax.tick_params(axis="both", labelsize=13.5, length=6)
    plot_ax.legend(loc="upper right", fontsize=11.5, frameon=True, facecolor="white", edgecolor="#D1D5DB", framealpha=0.96)
    plot_ax.text(0.965, 0.690, r"$\bf{\it{sPHENIX}}$ Internal", transform=plot_ax.transAxes, ha="right", va="top", fontsize=14.8, color=INK)
    plot_ax.text(
        0.965,
        0.625,
        r"$p{+}p$ $\sqrt{s}=200$ GeV" "\n" r"$|\eta|<0.7$, $22<E_T<28$ GeV" "\n" "same-row direct model scores",
        transform=plot_ax.transAxes,
        ha="right",
        va="top",
        fontsize=11.8,
        color=INK,
        linespacing=1.15,
    )
    plot_fig.tight_layout(pad=1.2)
    plot_fig.savefig(plot_png, dpi=180)
    plt.close(plot_fig)

    payload = {
        "schema": "PP_BASEV3E_DIRECT_TMVA_SCORE_SHAPE_SLIDE20_V1",
        "output_png": str(out_png.resolve()),
        "plot_only_png": str(plot_png.resolve()),
        "canvas_px": [2560, 1440],
        "source_csv": str(args.csv.resolve()),
        "summary_json": str(args.summary_json.resolve()) if args.summary_json else None,
        "rows": int(len(frame)),
        "signal_rows": int(len(signal)),
        "inclusive_rows": int(len(inclusive)),
        "signal_samples": signal_samples,
        "inclusive_samples": inclusive_samples,
        "stats": stats,
        "note": "PPG12 curves are direct Shuhang split-TMVA scores on the same rows, not projected ROOT histograms.",
    }
    out_json.write_text(json.dumps(payload, indent=2) + "\n")
    print(out_png)
    print(plot_png)
    print(out_json)


if __name__ == "__main__":
    main()

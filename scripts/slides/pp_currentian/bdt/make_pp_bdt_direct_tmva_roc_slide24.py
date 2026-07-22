#!/usr/bin/env python3
"""Build a slide-ready pp ROC comparison from same-row direct model scores."""

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
GREEN = "#047857"
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
            linewidth=1.15,
            edgecolor=ec,
            facecolor=fc,
            zorder=-1,
        )
    )


def text_block(
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


def wrapped_block(fig: plt.Figure, x: float, y: float, text: str, *, width: int, size: float, color: str = INK, linespacing: float = 1.18) -> None:
    wrapped = "\n".join(textwrap.fill(part, width=width) for part in text.split("\n"))
    text_block(fig, x, y, wrapped, size=size, color=color, linespacing=linespacing)


def roc_for(frame: pd.DataFrame, score_col: str) -> tuple[np.ndarray, np.ndarray, float]:
    y = (frame["class"] == "signal").astype(int).to_numpy()
    score = frame[score_col].to_numpy(dtype=float)
    valid = np.isfinite(score)
    y = y[valid]
    score = score[valid]
    order = np.argsort(-score, kind="mergesort")
    y_sorted = y[order]
    score_sorted = score[order]
    positives = float(np.sum(y_sorted == 1))
    negatives = float(np.sum(y_sorted == 0))
    if positives <= 0 or negatives <= 0:
        raise ValueError(f"Cannot compute ROC for {score_col}: positives={positives}, negatives={negatives}")

    distinct = np.r_[np.where(np.diff(score_sorted))[0], y_sorted.size - 1]
    tps = np.cumsum(y_sorted == 1)[distinct]
    fps = np.cumsum(y_sorted == 0)[distinct]
    tpr = np.r_[0.0, tps / positives, 1.0]
    fpr = np.r_[0.0, fps / negatives, 1.0]
    return fpr, tpr, float(np.trapezoid(tpr, fpr))


def load_summary(path: Path | None) -> dict:
    if path is None or not path.exists():
        return {}
    return json.loads(path.read_text())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--summary-json", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    set_style()
    args.outdir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.csv)
    summary = load_summary(args.summary_json)
    required = {"class", "sample", "cluster_Et", "our_score", "shuhang_split_score"}
    missing = required.difference(frame.columns)
    if missing:
        raise SystemExit(f"Missing required columns in {args.csv}: {sorted(missing)}")

    frame = frame[frame["class"].isin(["signal", "inclusive"])].copy()
    fpr_our, tpr_our, auc_our = roc_for(frame, "our_score")
    fpr_ppg12, tpr_ppg12, auc_ppg12 = roc_for(frame, "shuhang_split_score")

    file_counts = summary.get("file_counts", {})
    signal_samples = sorted(k.split("/", 1)[1] for k in file_counts if k.startswith("signal/"))
    inclusive_samples = sorted(k.split("/", 1)[1] for k in file_counts if k.startswith("inclusive/"))
    rows = len(frame)
    n_signal = int((frame["class"] == "signal").sum())
    n_inclusive = int((frame["class"] == "inclusive").sum())

    out_png = args.outdir / "pp_basev3e_bdt_roc_direct_tmva_slide24_candidate.png"
    out_json = args.outdir / "pp_basev3e_bdt_roc_direct_tmva_slide24_candidate.json"

    def natural_sample_key(sample: str) -> tuple[str, int]:
        head = "".join(ch for ch in sample if not ch.isdigit())
        digits = "".join(ch for ch in sample if ch.isdigit())
        return head, int(digits or 0)

    signal_samples = sorted(signal_samples, key=natural_sample_key)
    inclusive_samples = sorted(inclusive_samples, key=natural_sample_key)

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.subplots_adjust(0, 0, 1, 1)

    text_block(fig, 0.045, 0.945, "pp baseV3E ROC: Direct Model Application", size=31, weight="bold")
    text_block(
        fig,
        0.045,
        0.902,
        "Same candidate rows scored with this analysis XGBoost model and Shuhang/PPG12 split TMVA baseV3E model.",
        size=15.4,
        color=MUTED,
    )

    ax = fig.add_axes([0.060, 0.150, 0.535, 0.685])
    ax.plot(fpr_our, tpr_our, color=BLUE, lw=3.2, label=f"This analysis XGBoost, AUC = {auc_our:.3f}")
    ax.plot(fpr_ppg12, tpr_ppg12, color=RED, lw=3.2, ls="--", label=f"PPG12 split TMVA, AUC = {auc_ppg12:.3f}")
    ax.plot([0, 1], [0, 1], color="#9CA3AF", lw=1.6, ls=":", label="random classifier")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Inclusive MC efficiency", fontsize=20)
    ax.set_ylabel("Signal MC efficiency", fontsize=20)
    ax.set_title("Same-row ROC", loc="left", fontsize=22, fontweight="bold", color=INK, pad=10)
    ax.grid(color=GRID, lw=0.9, alpha=0.70)
    ax.tick_params(axis="both", labelsize=15, length=7)
    ax.legend(loc="lower right", frameon=True, facecolor="white", edgecolor="#D1D5DB", framealpha=0.95, fontsize=14.5)

    ax.text(
        0.965,
        0.940,
        r"$\bf{\it{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=18.8,
        color=INK,
    )
    ax.text(
        0.965,
        0.878,
        r"Pythia $\sqrt{s}=200$ GeV" "\n" r"$|\eta|<0.7$, $22<E_T<28$ GeV",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=15.5,
        color=INK,
        linespacing=1.18,
    )

    rhs_box_x = 0.625
    rhs_x = 0.646
    rhs_w = 0.345
    rounded_box(fig, rhs_box_x, 0.695, rhs_w, 0.155, SOFT_BLUE, "#BFDBFE")
    text_block(fig, rhs_x, 0.823, "What is plotted", size=20.2, color=BLUE, weight="bold")
    wrapped_block(
        fig,
        rhs_x,
        0.780,
        "Same 22-28 GeV candidates scored row-by-row with each model.\n\n"
        "Red dashed: direct PPG12 split-TMVA baseV3E score, not projected ROOT histograms.",
        width=62,
        size=12.8,
        linespacing=1.18,
    )

    rounded_box(fig, rhs_box_x, 0.475, rhs_w, 0.170, SOFT_GREEN, "#A7F3D0")
    text_block(fig, rhs_x, 0.617, "Direct-model result", size=20.2, color=GREEN, weight="bold")
    wrapped_block(
        fig,
        rhs_x,
        0.572,
        f"This analysis AUC = {auc_our:.3f}\n"
        f"PPG12 split TMVA AUC = {auc_ppg12:.3f}\n\n"
        "Direct scoring puts both models in the expected high-AUC baseV3E range.",
        width=58,
        size=13.2,
        linespacing=1.22,
    )

    rounded_box(fig, rhs_box_x, 0.235, rhs_w, 0.190, "white", PANEL_EDGE)
    text_block(fig, rhs_x, 0.394, "Same-row audit inputs", size=19.0, color=INK, weight="bold")
    signal_line = ", ".join(signal_samples) if signal_samples else "photonjet5, photonjet10, photonjet20"
    inclusive_line = ", ".join(inclusive_samples) if inclusive_samples else "jet8, jet12, jet20, jet30, jet40"
    input_text = (
        f"{rows:,} candidates: {n_signal:,} signal and {n_inclusive:,} inclusive\n"
        f"Signal MC: {signal_line}\n"
        f"Inclusive MC: {inclusive_line}\n"
        "Scores are computed from model files, not pre-made BDT histograms."
    )
    wrapped_block(
        fig,
        rhs_x,
        0.352,
        input_text,
        width=59,
        size=12.1,
        linespacing=1.27,
    )

    rounded_box(fig, rhs_box_x, 0.070, rhs_w, 0.115, SOFT_AMBER, "#FDBA74")
    fig.text(rhs_x, 0.150, "Takeaway:", ha="left", va="top", fontsize=13.8, fontweight="bold", color=INK)
    fig.text(
        0.735,
        0.150,
        "similar ROC curve and AUC, slight\n"
        "variations between them need\n"
        "further study.",
        ha="left",
        va="top",
        fontsize=12.5,
        color=INK,
        linespacing=1.15,
    )

    fig.savefig(out_png, dpi=160, bbox_inches=None, pad_inches=0)
    plt.close(fig)

    out_json.write_text(
        json.dumps(
            {
                "schema": "PP_BASEV3E_DIRECT_TMVA_ROC_SLIDE24_V1",
                "output_png": str(out_png),
                "canvas_px": [2560, 1440],
                "source_csv": str(args.csv),
                "summary_json": str(args.summary_json) if args.summary_json else None,
                "rows": rows,
                "signal_rows": n_signal,
                "inclusive_rows": n_inclusive,
                "auc_this_analysis_xgboost": auc_our,
                "auc_ppg12_split_tmva_direct": auc_ppg12,
                "signal_samples": signal_samples,
                "inclusive_samples": inclusive_samples,
                "note": "ROC curves are computed from same-row direct model scores, not projected ROOT histograms.",
            },
            indent=2,
        )
        + "\n"
    )
    print(out_png)
    print(out_json)


if __name__ == "__main__":
    main()

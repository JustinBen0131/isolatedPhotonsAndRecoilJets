#!/usr/bin/env python3
"""Render the direct, same-row PPG12 baseV3E ROC comparison as a slide-use PNG.

The input CSV must be the direct score-equivalence sample: identical accepted
MC candidates scored by the analysis retraining and the out-of-box PPG12
baseV3E TMVA.  This helper deliberately computes the ROC from the raw rows,
rather than from a pre-rendered diagnostic or a slide screenshot.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def roc_curve(labels: np.ndarray, scores: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Return FPR, TPR and AUC for binary labels, with no sklearn dependency."""
    order = np.argsort(-scores, kind="mergesort")
    labels = labels[order].astype(bool)
    positives = int(labels.sum())
    negatives = int((~labels).sum())
    if positives == 0 or negatives == 0:
        raise ValueError("ROC requires both signal and inclusive-background rows.")

    true_positive = np.cumsum(labels)
    false_positive = np.cumsum(~labels)
    distinct = np.r_[np.where(np.diff(scores[order]))[0], len(scores) - 1]
    tpr = np.r_[0.0, true_positive[distinct] / positives]
    fpr = np.r_[0.0, false_positive[distinct] / negatives]
    return fpr, tpr, float(np.trapezoid(tpr, fpr))


def value_at_signal_efficiency(
    fpr: np.ndarray, tpr: np.ndarray, target: float
) -> float:
    """Interpolate background acceptance at a fixed signal efficiency."""
    unique_tpr, first_indices = np.unique(tpr, return_index=True)
    return float(np.interp(target, unique_tpr, fpr[first_indices]))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, required=True, help="Direct same-row score sample CSV")
    parser.add_argument("--audit", type=Path, required=True, help="Score-contract audit JSON")
    parser.add_argument("--outdir", type=Path, required=True, help="Destination for PNG and manifest")
    args = parser.parse_args()

    frame = pd.read_csv(args.csv)
    required = {"class", "our_score", "shuhang_split_score"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    labels = frame["class"].eq("signal").to_numpy(dtype=bool)
    fpr_ours, tpr_ours, auc_ours = roc_curve(labels, frame["our_score"].to_numpy(float))
    fpr_ppg12, tpr_ppg12, auc_ppg12 = roc_curve(
        labels, frame["shuhang_split_score"].to_numpy(float)
    )
    wp_signal_efficiency = 0.80
    wp_ours = value_at_signal_efficiency(fpr_ours, tpr_ours, wp_signal_efficiency)
    wp_ppg12 = value_at_signal_efficiency(fpr_ppg12, tpr_ppg12, wp_signal_efficiency)

    audit = json.loads(args.audit.read_text())
    args.outdir.mkdir(parents=True, exist_ok=True)
    png_path = args.outdir / "pp_basev3e_same_row_roc_overlay.png"
    manifest_path = args.outdir / "pp_basev3e_same_row_roc_overlay.json"

    plt.rcParams.update(
        {
            "font.family": ["Helvetica Neue", "Arial", "DejaVu Sans"],
            "font.size": 16,
            "axes.labelsize": 18,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 13,
            "axes.linewidth": 1.25,
            "mathtext.fontset": "dejavusans",
        }
    )
    fig, ax = plt.subplots(figsize=(8.6, 7.0), constrained_layout=True)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    blue = "#1976B9"
    red = "#D62828"
    ax.plot([0, 1], [0, 1], color="#9AA0A6", linestyle=":", linewidth=1.2, zorder=1)
    ax.plot(
        fpr_ours,
        tpr_ours,
        color=blue,
        linewidth=3.0,
        label=f"This-analysis retraining, AUC = {auc_ours:.3f}",
        zorder=3,
    )
    ax.plot(
        fpr_ppg12,
        tpr_ppg12,
        color=red,
        linewidth=2.8,
        linestyle="--",
        label=f"PPG12 baseV3E (out-of-box TMVA), AUC = {auc_ppg12:.3f}",
        zorder=3,
    )
    ax.scatter([wp_ours], [wp_signal_efficiency], s=68, marker="s", color=blue, edgecolor="white", linewidth=0.8, zorder=4)
    ax.scatter([wp_ppg12], [wp_signal_efficiency], s=68, marker="s", color=red, edgecolor="white", linewidth=0.8, zorder=4)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Background acceptance")
    ax.set_ylabel("Signal efficiency")
    ax.set_xticks(np.arange(0, 1.01, 0.2))
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.tick_params(direction="in", top=True, right=True)
    ax.grid(True, color="#D6D6D6", linewidth=0.55, alpha=0.65)

    # The upper-middle region is intentionally used: it avoids the steep ROC
    # turn-on at upper left while retaining the standard sPHENIX in-plot label.
    ax.text(0.53, 0.935, r"$\bf{sPHENIX}$ Internal", transform=ax.transAxes, va="top", fontsize=14)
    ax.text(
        0.53,
        0.875,
        r"$p{+}p$, $\sqrt{s}=200$ GeV" "\n" r"$22 < E_T < 28$ GeV, $|\eta| < 0.7$",
        transform=ax.transAxes,
        va="top",
        fontsize=12.5,
    )
    ax.text(
        0.955,
        0.455,
        "80% signal-efficiency points\n"
        f"This analysis: $\\epsilon_{{\\mathrm{{bkg}}}}$ = {wp_ours:.3f}\n"
        f"PPG12 TMVA: $\\epsilon_{{\\mathrm{{bkg}}}}$ = {wp_ppg12:.3f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=11.3,
        bbox={"boxstyle": "round,pad=0.32", "facecolor": "white", "edgecolor": "#B8B8B8", "alpha": 0.96},
    )
    ax.legend(loc="lower right", frameon=False, handlelength=2.5, borderaxespad=0.6)
    fig.savefig(png_path, dpi=300, facecolor="white")
    plt.close(fig)

    manifest = {
        "schema": "JSTG_PPG12_SAME_ROW_ROC_V1",
        "plot": str(png_path),
        "input_csv": str(args.csv),
        "input_csv_sha256": sha256(args.csv),
        "audit_json": str(args.audit),
        "selection": {"pt_range_GeV": audit["pt_range"], "abs_eta_max": audit["eta_max"]},
        "sample": audit.get("sample"),
        "comparison": {
            "analysis_model": "This-analysis baseV3E retraining (score_ours)",
            "ppg12_model": "PPG12 baseV3E out-of-box split-TMVA (score_shuhang_split)",
            "same_row_candidate_count": int(len(frame)),
            "signal_count": int(labels.sum()),
            "inclusive_count": int((~labels).sum()),
        },
        "metrics_recomputed_from_raw_scores": {
            "analysis_auc": auc_ours,
            "ppg12_auc": auc_ppg12,
            "signal_efficiency_working_point": wp_signal_efficiency,
            "analysis_background_acceptance_at_working_point": wp_ours,
            "ppg12_background_acceptance_at_working_point": wp_ppg12,
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {png_path}")
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Restyle the BDT/MLP/LogReg ROC overlay for slide use."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPO = Path(__file__).resolve().parents[1]
SOURCE = (
    REPO
    / "dataOutput/auauMLPDiagnosticPlots/"
    / "roc_overlay_bdt_mlp_logreg_fullstat_15to35_20260513"
)
OUT_DIR = (
    REPO
    / "dataOutput/auauMLPDiagnosticPlots/"
    / "roc_overlay_bdt_mlp_logreg_slide43_20260515"
)

CURVE_POINTS = SOURCE / "roc_overlay_bdt_mlp_logreg_fullstat_curve_points.csv"
AUC_SUMMARY = SOURCE / "roc_overlay_bdt_mlp_logreg_fullstat_auc_summary.csv"
METADATA = SOURCE / "roc_overlay_bdt_mlp_logreg_fullstat_metadata.json"

CENT_ORDER = ["0-20%", "20-50%", "50-80%"]
MODEL_ORDER = ["BDT", "MLP", "LogReg"]
MODEL_LABEL = {
    "BDT": "BDT",
    "MLP": "MLP",
    "LogReg": "Logistic regression",
}
COLORS = {
    "BDT": "#0072B2",
    "MLP": "#D55E00",
    "LogReg": "#009E73",
}


def read_curves() -> dict[tuple[str, str], tuple[list[float], list[float]]]:
    curves: dict[tuple[str, str], tuple[list[float], list[float]]] = defaultdict(lambda: ([], []))
    with CURVE_POINTS.open() as handle:
        for row in csv.DictReader(handle):
            fpr, tpr = curves[(row["centrality"], row["model"])]
            fpr.append(float(row["fpr"]))
            tpr.append(float(row["tpr"]))
    return curves


def read_auc() -> dict[tuple[str, str], dict]:
    out = {}
    with AUC_SUMMARY.open() as handle:
        for row in csv.DictReader(handle):
            out[(row["centrality"], row["model"])] = row
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    curves = read_curves()
    auc = read_auc()
    metadata = json.loads(METADATA.read_text()) if METADATA.exists() else {}

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.linewidth": 1.15,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.95), sharex=True, sharey=True)

    for ax, cent in zip(axes, CENT_ORDER, strict=True):
        for model in MODEL_ORDER:
            key = (cent, model)
            fpr, tpr = curves[key]
            info = auc[key]
            ax.plot(
                fpr,
                tpr,
                lw=2.7,
                color=COLORS[model],
                label=f"{MODEL_LABEL[model]}  AUC={float(info['auc']):.3f}",
            )
        ax.plot([0, 1], [0, 1], color="0.70", lw=1.15, ls="--", zorder=0)
        info0 = auc[(cent, "BDT")]
        ax.text(
            0.045,
            0.065,
            f"N={int(info0['entries']):,}\nS={int(info0['signal_entries']):,}, B={int(info0['background_entries']):,}",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=8.8,
            color="0.25",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.8},
        )
        ax.set_title(f"{cent} centrality", fontsize=13.5, fontweight="bold", pad=10)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.grid(True, color="0.87", lw=0.65)
        ax.set_xlabel("False positive rate / background efficiency", fontsize=11)
        ax.tick_params(labelsize=10.5, length=6)
        ax.minorticks_on()
        ax.legend(
            loc="lower right",
            fontsize=9.0,
            frameon=True,
            framealpha=0.94,
            edgecolor="0.75",
            borderpad=0.55,
        )

    axes[0].set_ylabel("True positive rate / signal efficiency", fontsize=11)

    fig.text(0.045, 0.965, r"$\bf{\it{sPHENIX}}$ Internal", ha="left", va="top", fontsize=14.5)
    fig.text(0.045, 0.915, "Au+Au embedded validation", ha="left", va="top", fontsize=10.8)
    fig.text(0.50, 0.965, "ROC curves by centrality", ha="center", va="top", fontsize=18, fontweight="bold")
    fig.text(
        0.50,
        0.912,
        r"$15 < E_T < 35$ GeV",
        ha="center",
        va="top",
        fontsize=12.5,
        color="0.28",
    )

    fig.tight_layout(rect=[0.035, 0.070, 0.995, 0.835], w_pad=1.9)

    out_png = OUT_DIR / "roc_overlay_bdt_mlp_logreg_by_centrality_pt15to35.png"
    out_meta = OUT_DIR / "roc_overlay_bdt_mlp_logreg_by_centrality_pt15to35_metadata.json"
    out_auc = OUT_DIR / "roc_overlay_bdt_mlp_logreg_by_centrality_pt15to35_auc_summary.csv"
    fig.savefig(out_png, dpi=220)
    plt.close(fig)

    payload = {
        "schema": "AUAU_ROC_OVERLAY_BDT_MLP_LOGREG_SLIDE43_V1",
        "source_curve_points": str(CURVE_POINTS.relative_to(REPO)),
        "source_auc_summary": str(AUC_SUMMARY.relative_to(REPO)),
        "source_metadata": str(METADATA.relative_to(REPO)),
        "selection": metadata.get("selection", "15 < cluster_Et < 35 GeV"),
        "rows_seen": metadata.get("rows_seen"),
        "rows_selected_pt15to35": metadata.get("rows_selected_pt15to35"),
        "note": "Plot is a visual restyle of the existing full-stat BDT/MLP/LogReg ROC curves; curve points and AUC values are unchanged.",
    }
    out_meta.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    out_auc.write_text(AUC_SUMMARY.read_text())
    print(out_png)
    print(out_auc)
    print(out_meta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

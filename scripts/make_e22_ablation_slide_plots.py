#!/usr/bin/env python3
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auauMLDiagnosticRuns/"
    "global_etcent_inclusive3_sixpack_20260516_135439"
)
PULL_ROOT = ROOT / "validation" / "basev3e_w33_e22ratio_20260518_192630"
MODEL_DIR = PULL_ROOT / "bdt_models"
VALIDATION_DIR = PULL_ROOT / "validation"
OUTDIR = ROOT / "slideReady" / "basev3e_w33_e22ratio_diagnostics"

LABELS = {
    "cluster_Et": r"cluster $E_T$",
    "cluster_Eta": r"cluster $\eta$",
    "vertexz": r"$z_{vtx}$",
    "centrality": "centrality",
    "cluster_weta_cogx": r"$w_{\eta}$",
    "cluster_wphi_cogx": r"$w_{\phi}$",
    "cluster_weta33_cogx": r"$w_{\eta}^{3x3}$",
    "cluster_wphi33_cogx": r"$w_{\phi}^{3x3}$",
    "cluster_et1": "cluster et1",
    "cluster_et2": "cluster et2",
    "cluster_et3": "cluster et3",
    "cluster_et4": "cluster et4",
    "e11_over_e33": r"$E_{11}/E_{33}$",
    "e32_over_e35": r"$E_{32}/E_{35}$",
    "e22_over_e37": r"$E_{22}/E_{37}$",
    "e22_over_e53": r"$E_{22}/E_{53}$",
}


def read_json_with_preamble(path: Path) -> dict:
    text = path.read_text()
    return json.loads(text[text.find("{") :])


PRODUCT_CONFIG = {
    "baseBDT_v3E_withCentrality_w33_E22E37": {
        "suffix": "E22E37",
        "title": r"$E_{22}/E_{37}$ ablation BDT",
        "input_line": r"Inputs: base v3E + centrality + 3x3 widths + $E_{22}/E_{37}$",
        "overlay_subtitle": r"Rows ordered by split-gain usage in the $E_{22}/E_{37}$ ablation BDT",
    },
    "baseBDT_v3E_withCentrality_w33_E22E37_E22E53": {
        "suffix": "E22E37_E22E53",
        "title": "two-ratio ablation BDT",
        "input_line": r"Inputs: base v3E + centrality + 3x3 widths + $E_{22}/E_{37}$ + $E_{22}/E_{53}$",
        "overlay_subtitle": r"Rows ordered by split-gain usage in the two-ratio ablation BDT",
    },
}


def split_gain(product: str) -> pd.DataFrame:
    meta = json.loads((MODEL_DIR / f"auau_tight_bdt_{product}_tmva.metadata.json").read_text())
    features = list(meta["features"])
    model = read_json_with_preamble(MODEL_DIR / f"auau_tight_bdt_{product}_tmva.xgb.json")
    gain = defaultdict(float)
    count = defaultdict(int)
    for tree in model["learner"]["gradient_booster"]["model"]["trees"]:
        for left, right, idx, loss in zip(
            tree["left_children"],
            tree["right_children"],
            tree["split_indices"],
            tree["loss_changes"],
        ):
            if left < 0 and right < 0:
                continue
            idx = int(idx)
            if 0 <= idx < len(features):
                feature = features[idx]
                gain[feature] += max(float(loss), 0.0)
                count[feature] += 1
    total = sum(gain.values())
    rows = [
        {
            "feature": feature,
            "split_gain": gain[feature],
            "split_count": count[feature],
            "split_gain_fraction": gain[feature] / total if total else 0.0,
        }
        for feature in features
    ]
    return pd.DataFrame(rows).sort_values("split_gain_fraction", ascending=False)


def style_axis(ax) -> None:
    ax.tick_params(axis="both", direction="in", top=True, right=True, length=5.5, width=1.15)
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)
        spine.set_color("0.12")
    ax.grid(axis="x", color="0.88", linewidth=0.9)
    ax.set_axisbelow(True)


def add_sphenix(fig) -> None:
    fig.text(0.865, 0.946, "sPHENIX", ha="right", va="bottom", fontsize=22, fontstyle="italic", fontweight="bold")
    fig.text(0.965, 0.946, "Internal", ha="right", va="bottom", fontsize=22)


def make_split_gain_plot(product: str, df: pd.DataFrame) -> Path:
    cfg = PRODUCT_CONFIG[product]
    out = OUTDIR / f"basev3e_w33_{cfg['suffix']}_split_gain_usage.png"
    plot_df = df.sort_values("split_gain_fraction", ascending=True).copy()
    plot_df["fraction_percent"] = 100.0 * plot_df["split_gain_fraction"]
    fig, ax = plt.subplots(figsize=(14.4, 9.6), dpi=190)
    fig.subplots_adjust(left=0.17, right=0.975, top=0.76, bottom=0.12)
    y = np.arange(len(plot_df))
    colors = plt.cm.Blues(np.linspace(0.34, 0.90, len(plot_df)))
    bars = ax.barh(y, plot_df["fraction_percent"], color=colors, edgecolor="none", height=0.66)
    ax.set_yticks(y)
    ax.set_yticklabels([LABELS.get(v, v) for v in plot_df["feature"]], fontsize=16)
    ax.set_xlabel("Fraction of total split gain [%]", fontsize=20, labelpad=10)
    ax.tick_params(axis="x", labelsize=15)
    xmax = max(10.0, float(plot_df["fraction_percent"].max()) * 1.22)
    ax.set_xlim(0.0, xmax)
    style_axis(ax)
    for bar, (_, row) in zip(bars, plot_df.iterrows()):
        value = float(row["fraction_percent"])
        label = f"{value:.1f}% ({int(row['split_count'])})"
        x = min(value + 0.35, xmax - 4.8)
        ha = "left"
        if value > xmax - 5.8:
            x = xmax - 0.35
            ha = "right"
        ax.text(x, bar.get_y() + bar.get_height() / 2, label, va="center", ha=ha, fontsize=14.8, fontweight="bold", clip_on=True)
    fig.text(0.145, 0.885, f"XGBoost split-gain usage: {cfg['title']}", ha="left", fontsize=21.0, fontweight="bold")
    fig.text(0.145, 0.852, cfg["input_line"], ha="left", fontsize=15.0, color="0.25")
    fig.text(0.145, 0.820, "gain summed over all tree splits; labels give gain fraction and split count", ha="left", fontsize=14.8, color="0.30")
    fig.text(0.145, 0.790, r"Embedded Photon12+20 vs embedded inclusive Jet12+20+30, 15 < cluster $E_T$ < 35 GeV", ha="left", fontsize=13.2, color="0.22")
    add_sphenix(fig)
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    return out


def make_energy_ratio_overlay(product: str, df_gain: pd.DataFrame) -> Path:
    cfg = PRODUCT_CONFIG[product]
    hist_path = VALIDATION_DIR / f"energy_ratio_histograms_{cfg['suffix']}_fullstat.csv"
    hist = pd.read_csv(hist_path)
    ratio_order = [
        feature
        for feature in df_gain["feature"].tolist()
        if feature.startswith("e") and "_over_" in feature and feature in set(hist["feature"])
    ]
    out = OUTDIR / f"basev3e_w33_{cfg['suffix']}_energy_ratio_overlays_by_cent.png"
    cent_order = ["0-20%", "20-50%", "50-80%"]
    fig, axes = plt.subplots(len(ratio_order), len(cent_order), figsize=(15.6, 11.2), dpi=190, sharex=True)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.690, bottom=0.085, wspace=0.14, hspace=0.25)
    colors = {"signal": "#1f77b4", "background": "#d95f02"}
    labels = {"signal": "signal", "background": "background"}
    for r, feature in enumerate(ratio_order):
        row_data = hist[hist["feature"].eq(feature)]
        ymax = max(float(row_data["density"].max()) * 1.18, 0.1)
        for c, cent in enumerate(cent_order):
            ax = axes[r, c] if len(ratio_order) > 1 else axes[c]
            for cls in ["signal", "background"]:
                sub = row_data[row_data["cent_bin"].eq(cent) & row_data["class"].eq(cls)].sort_values("bin_low")
                x = 0.5 * (sub["bin_low"].to_numpy() + sub["bin_high"].to_numpy())
                y = sub["density"].to_numpy()
                ax.step(x, y, where="mid", color=colors[cls], linewidth=2.15, label=labels[cls])
                ax.fill_between(x, y, step="mid", color=colors[cls], alpha=0.12)
            ax.set_ylim(0, ymax)
            ax.set_xlim(0.0, 1.2)
            ax.tick_params(axis="both", labelsize=11.5)
            ax.grid(color="0.90", linewidth=0.8)
            for spine in ax.spines.values():
                spine.set_linewidth(1.0)
                spine.set_color("0.18")
            if r == 0:
                ax.set_title(cent, fontsize=14.5, fontweight="bold", pad=7)
            if c == 0:
                gain = 100.0 * float(df_gain.set_index("feature").loc[feature, "split_gain_fraction"])
                ax.set_ylabel(f"{LABELS.get(feature, feature)}\nDensity", fontsize=12.5)
                ax.text(0.03, 0.90, f"gain {gain:.1f}%", transform=ax.transAxes, fontsize=10.8, color="0.25")
            if r == len(ratio_order) - 1:
                ax.set_xlabel("Feature value", fontsize=12.8)
    handles, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, lab, loc="upper center", bbox_to_anchor=(0.58, 0.755), ncol=2, frameon=False, fontsize=14.0)
    fig.text(0.075, 0.905, "Signal/background input-shape overlays for ratio variables", ha="left", fontsize=20.5, fontweight="bold")
    fig.text(0.075, 0.870, cfg["overlay_subtitle"], ha="left", fontsize=14.5, color="0.30")
    fig.text(0.075, 0.840, r"Normalized full-stat validation distributions, 15 < cluster $E_T$ < 35 GeV", ha="left", fontsize=13.2, color="0.25")
    add_sphenix(fig)
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    return out


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for product, cfg in PRODUCT_CONFIG.items():
        df = split_gain(product)
        df.to_csv(OUTDIR / f"basev3e_w33_{cfg['suffix']}_split_gain.csv", index=False)
        print(make_split_gain_plot(product, df))
        hist_path = VALIDATION_DIR / f"energy_ratio_histograms_{cfg['suffix']}_fullstat.csv"
        if hist_path.exists():
            print(make_energy_ratio_overlay(product, df))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auauMLDiagnosticRuns/"
    "global_etcent_inclusive3_sixpack_20260516_135439/slideReady/"
    "basev3e_w33_feature_diagnostics"
)
CSV_PATH = ROOT / "basev3e_w33_split_gain.csv"
OUT_PATH = ROOT / "basev3e_w33_split_gain_usage_slide40_style.png"


LABELS = {
    "e32_over_e35": r"$E_{32}/E_{35}$",
    "cluster_et1": "cluster et1",
    "cluster_wphi33_cogx": r"$w_{\phi}^{3x3}$",
    "cluster_weta_cogx": r"$w_{\eta}$",
    "cluster_weta33_cogx": r"$w_{\eta}^{3x3}$",
    "cluster_wphi_cogx": r"$w_{\phi}$",
    "centrality": "centrality",
    "cluster_et3": "cluster et3",
    "e11_over_e33": r"$E_{11}/E_{33}$",
    "cluster_et2": "cluster et2",
    "cluster_Et": r"cluster $E_T$",
    "cluster_Eta": r"cluster $\eta$",
    "cluster_et4": "cluster et4",
    "vertexz": r"$z_{vtx}$",
}


def main() -> None:
    df = pd.read_csv(CSV_PATH).copy()
    df["fraction_percent"] = 100.0 * df["split_gain_fraction"]
    df = df.sort_values("fraction_percent", ascending=True)

    fig, ax = plt.subplots(figsize=(14.4, 8.1), dpi=180)
    fig.subplots_adjust(left=0.145, right=0.968, top=0.695, bottom=0.14)

    y = range(len(df))
    colors = plt.cm.Blues([0.35 + 0.55 * i / max(len(df) - 1, 1) for i in range(len(df))])
    bars = ax.barh(y, df["fraction_percent"], color=colors, edgecolor="none", height=0.72)

    ax.set_yticks(list(y))
    ax.set_yticklabels([LABELS.get(v, v) for v in df["feature"]], fontsize=16)
    ax.set_xlabel("Fraction of total split gain [%]", fontsize=20, labelpad=10)
    ax.set_xlim(0.0, 39.0)
    ax.set_xticks(range(0, 40, 5))
    ax.tick_params(axis="x", labelsize=16, direction="in", top=True, length=7, width=1.3)
    ax.tick_params(axis="y", labelsize=16, direction="in", right=True, length=7, width=1.3)
    ax.grid(axis="x", color="0.86", linewidth=1.2)
    ax.set_axisbelow(True)

    for spine in ax.spines.values():
        spine.set_linewidth(1.4)
        spine.set_color("0.12")

    xmax = ax.get_xlim()[1]
    for bar, (_, row) in zip(bars, df.iterrows()):
        value = float(row["fraction_percent"])
        count = int(row["split_count"])
        label = f"{value:.1f}%  ({count} splits)"
        x = min(value + 0.55, xmax - 6.2)
        ha = "left"
        if value > xmax - 7.5:
            x = xmax - 0.55
            ha = "right"
        ax.text(
            x,
            bar.get_y() + bar.get_height() / 2.0,
            label,
            va="center",
            ha=ha,
            fontsize=15.5,
            fontweight="bold",
            color="0.05",
            clip_on=True,
        )

    fig.text(
        0.155,
        0.835,
        "XGBoost split-gain usage: base v3E + centrality + 3x3 widths",
        ha="left",
        va="bottom",
        fontsize=22,
        fontweight="bold",
    )
    fig.text(
        0.155,
        0.802,
        "gain summed over all tree splits in the trained global BDT",
        ha="left",
        va="bottom",
        fontsize=16,
        color="0.30",
    )
    fig.text(
        0.155,
        0.768,
        r"Photon12+20 & Jet12+20+30, 15 < cluster $E_T$ < 35 GeV",
        ha="left",
        va="bottom",
        fontsize=15,
        color="0.22",
    )
    fig.text(
        0.865,
        0.942,
        "sPHENIX",
        ha="right",
        va="bottom",
        fontsize=23,
        fontstyle="italic",
        fontweight="bold",
    )
    fig.text(0.965, 0.942, "Internal", ha="right", va="bottom", fontsize=23)
    fig.savefig(OUT_PATH, facecolor="white")
    plt.close(fig)
    print(OUT_PATH)


if __name__ == "__main__":
    main()

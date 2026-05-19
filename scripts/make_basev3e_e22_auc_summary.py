#!/usr/bin/env python3
"""Make a slide-ready AUC summary for the base-v3E E22-ratio ablation."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


OUT_DIR = Path(
    "dataOutput/auauMLDiagnosticRuns/"
    "global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/basev3e_w33_e22ratio_diagnostics"
)

CENTRALITY_LABEL = "0-20%"
CENTRALITY_KEY = "0_20"
FINITE_ENTRIES = 1_397_721
SIGNAL_ENTRIES = 1_197_217
BACKGROUND_ENTRIES = 200_504


ROWS = [
    {
        "model": "Base v3E + centrality + 3x3 widths",
        "short": "base",
        "inputs": 14,
        "auc": 0.7947099269024895,
        "centrality_bin": CENTRALITY_KEY,
        "provenance": "0-20% AUC from basev3e_controls_20260518_1110 centrality-resolved validation diagnostics.",
    },
    {
        "model": "+ E22/E37",
        "short": "E22/E37",
        "inputs": 15,
        "auc": 0.8039423313452887,
        "centrality_bin": CENTRALITY_KEY,
        "provenance": "0-20% AUC from basev3e_w33_e22ratio_20260518_192630 validation metrics.",
    },
    {
        "model": "+ E22/E53",
        "short": "E22/E53",
        "inputs": 15,
        "auc": 0.8014658462098146,
        "centrality_bin": CENTRALITY_KEY,
        "provenance": "0-20% AUC from basev3e_w33_e22ratio_20260518_192630 validation metrics.",
    },
    {
        "model": "+ E22/E37 and E22/E53",
        "short": "both ratios",
        "inputs": 16,
        "auc": 0.807601576784641,
        "centrality_bin": CENTRALITY_KEY,
        "provenance": "0-20% AUC from basev3e_w33_e22ratio_20260518_192630 validation metrics.",
    },
]


def add_sphenix_label(ax: plt.Axes) -> None:
    ax.text(
        0.865,
        1.105,
        "sPHENIX",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=20,
        fontstyle="italic",
        fontweight="bold",
        clip_on=False,
    )
    ax.text(
        0.868,
        1.105,
        " Internal",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=20,
        clip_on=False,
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(ROWS)
    base_auc = float(df.loc[0, "auc"])
    df["delta_auc"] = df["auc"] - base_auc
    df["delta_auc_x1e3"] = 1000.0 * df["delta_auc"]
    df["delta_auc_percent"] = 100.0 * df["delta_auc"] / base_auc
    df["rank"] = df["auc"].rank(method="min", ascending=False).astype(int)

    csv_path = OUT_DIR / "basev3e_w33_e22ratio_auc_summary.csv"
    df.to_csv(csv_path, index=False)

    colors = {
        "base": "#6B7280",
        "E22/E37": "#2563EB",
        "E22/E53": "#F97316",
        "both ratios": "#16A34A",
    }
    marker_edges = {
        "base": "#374151",
        "E22/E37": "#1E40AF",
        "E22/E53": "#C2410C",
        "both ratios": "#166534",
    }

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.labelsize": 20,
            "xtick.labelsize": 17,
            "ytick.labelsize": 18,
        }
    )

    fig, ax = plt.subplots(figsize=(16.0, 9.0), dpi=180)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Preserve the conceptual build-up order from top to bottom.
    plot_df = df.iloc[::-1].reset_index(drop=True)
    y_positions = list(range(len(plot_df)))

    xmin = min(df["auc"]) - 0.0032
    xmax = max(df["auc"]) + 0.0052
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-0.65, len(plot_df) - 0.35)

    for i, (_, row) in enumerate(plot_df.iterrows()):
        y = y_positions[i]
        if i % 2 == 0:
            ax.axhspan(y - 0.46, y + 0.46, color="#F8FAFC", zorder=0)

        auc = float(row["auc"])
        short = str(row["short"])
        color = colors[short]
        edge = marker_edges[short]

        ax.hlines(
            y,
            min(base_auc, auc),
            max(base_auc, auc),
            color=color,
            linewidth=9 if short != "base" else 0,
            alpha=0.32,
            zorder=1,
        )
        ax.scatter(
            [auc],
            [y],
            s=340,
            color=color,
            edgecolor=edge,
            linewidth=2.1,
            zorder=3,
        )

        delta = float(row["delta_auc_x1e3"])
        delta_percent = float(row["delta_auc_percent"])
        if short == "base":
            delta_text = "reference"
        else:
            delta_text = f"+{delta_percent:.2f}% gain in AUC"

        ax.text(
            auc + 0.00048,
            y + 0.13,
            f"AUC {auc:.4f}",
            ha="left",
            va="center",
            fontsize=20,
            fontweight="bold",
            color="#111827",
            zorder=4,
        )
        ax.text(
            auc + 0.00048,
            y - 0.20,
            delta_text,
            ha="left",
            va="center",
            fontsize=17,
            color=color if short != "base" else "#4B5563",
            zorder=4,
        )

        ax.text(
            xmin + 0.00045,
            y,
            f"{int(row['inputs'])} inputs",
            ha="left",
            va="center",
            fontsize=15,
            color="#4B5563",
            bbox={
                "facecolor": "white",
                "edgecolor": "#CBD5E1",
                "boxstyle": "round,pad=0.25",
                "linewidth": 0.8,
            },
        )

    ax.axvline(
        base_auc,
        color="#374151",
        linestyle=(0, (4, 4)),
        linewidth=2.0,
        alpha=0.9,
        zorder=2,
    )

    display_labels = {
        "Base v3E + centrality + 3x3 widths": "Base v3E + centrality\n+ 3x3 widths",
        "+ E22/E37": "+ E22/E37",
        "+ E22/E53": "+ E22/E53",
        "+ E22/E37 and E22/E53": "+ E22/E37\n+ E22/E53",
    }
    ax.set_yticks(y_positions)
    ax.set_yticklabels([display_labels[m] for m in plot_df["model"].tolist()])
    ax.tick_params(axis="y", length=0, pad=14)
    ax.set_xlabel(f"{CENTRALITY_LABEL} centrality validation AUC")
    ax.grid(axis="x", color="#D1D5DB", linewidth=1.0, alpha=0.85)
    ax.grid(axis="y", visible=False)
    ax.set_axisbelow(True)

    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#111827")
    ax.spines["bottom"].set_linewidth(1.2)

    fig.text(
        0.265,
        0.93,
        f"AUC gain from E22 energy-ratio inputs ({CENTRALITY_LABEL} centrality)",
        ha="left",
        va="bottom",
        fontsize=25,
        fontweight="bold",
    )
    fig.text(
        0.265,
        0.902,
        (
            "Photon12+20 signal vs Jet12+20+30 inclusive background; same XGBoost settings"
        ),
        ha="left",
        va="bottom",
        fontsize=15,
        color="#4B5563",
    )
    fig.text(
        0.265,
        0.875,
        f"AUCs are computed in the {CENTRALITY_LABEL} centrality bin",
        ha="left",
        va="bottom",
        fontsize=15,
        color="#4B5563",
    )
    ax.text(
        0.0,
        -0.18,
        (
            f"{CENTRALITY_LABEL} finite-score entries: {FINITE_ENTRIES:,} "
            f"({SIGNAL_ENTRIES:,} signal, {BACKGROUND_ENTRIES:,} background).\n"
            "Base row comes from base-controls centrality diagnostics; E22 rows come from the READY E22 sidecar validation summary."
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12.8,
        color="#6B7280",
        clip_on=False,
    )

    add_sphenix_label(ax)

    fig.subplots_adjust(left=0.265, right=0.945, top=0.80, bottom=0.23)
    png_path = OUT_DIR / "basev3e_w33_e22ratio_auc_summary.png"
    fig.savefig(png_path)
    plt.close(fig)

    print(png_path)
    print(csv_path)


if __name__ == "__main__":
    main()

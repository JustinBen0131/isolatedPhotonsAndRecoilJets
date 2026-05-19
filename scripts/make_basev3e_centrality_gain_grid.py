#!/usr/bin/env python3
"""Plot base-v3E AUC comparison with and without centrality in a 3x8 grid."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: E402


DEFAULT_CENT_ORDER = ["0-20", "20-50", "50-80"]
DEFAULT_ET_ORDER = ["15-17", "17-19", "19-21", "21-23", "23-25", "25-27", "27-30", "30-35"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--prefix", default="basev3e_centrality_auc_gain_3x8")
    parser.add_argument("--title", default="Base v3E BDT: AUC with centrality vs without centrality")
    return parser.parse_args()


def label_low(label: str) -> float:
    return float(str(label).split("-")[0])


def ordered_labels(values: pd.Series, fallback: list[str]) -> list[str]:
    labels = sorted({str(v) for v in values}, key=label_low)
    if set(labels) == set(fallback):
        return fallback
    return labels


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input)
    args.outdir.mkdir(parents=True, exist_ok=True)

    cent_order = ordered_labels(df["centrality_bin"], DEFAULT_CENT_ORDER)
    et_order = ordered_labels(df["et_bin"], DEFAULT_ET_ORDER)

    delta = np.full((len(cent_order), len(et_order)), np.nan)
    auc_with = np.full_like(delta, np.nan, dtype=float)
    auc_without = np.full_like(delta, np.nan, dtype=float)
    entries = np.zeros_like(delta, dtype=int)

    for _, row in df.iterrows():
        cent = str(row["centrality_bin"])
        et = str(row["et_bin"])
        if cent not in cent_order or et not in et_order:
            continue
        i = cent_order.index(cent)
        j = et_order.index(et)
        auc_with[i, j] = float(row["auc_with_centrality"])
        auc_without[i, j] = float(row["auc_without_centrality"])
        delta[i, j] = float(row["delta_with_minus_without"])
        entries[i, j] = int(row["entries"])

    # Slide labels show three decimals, so ties are defined at that precision.
    rounded_delta = np.round(auc_with, 3) - np.round(auc_without, 3)
    winner_grid = np.where(rounded_delta > 0.0, 2.0, np.where(rounded_delta < 0.0, 0.0, 1.0))
    cmap = ListedColormap(["#d95f5f", "#3c9f5f", "#2f78bd"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

    fig_w = 10.4
    fig_h = 2.55 + 0.62 * len(cent_order)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=260)
    ax.imshow(winner_grid, cmap=cmap, norm=norm, aspect="auto")

    ax.set_xticks(np.arange(len(et_order)))
    ax.set_xticklabels([x.replace("-", "-") for x in et_order], fontsize=10)
    ax.set_yticks(np.arange(len(cent_order)))
    ax.set_yticklabels([f"{x}%" for x in cent_order], fontsize=11, fontweight="bold")
    ax.set_xlabel(r"Cluster $E_T$ bin [GeV]", fontsize=11)
    ax.set_ylabel("Centrality", fontsize=11)
    ax.tick_params(which="both", length=0)

    for x in np.arange(-0.5, len(et_order), 1):
        ax.axvline(x, color="white", lw=1.4)
    for y in np.arange(-0.5, len(cent_order), 1):
        ax.axhline(y, color="white", lw=1.4)

    fontsize = 7.8 if len(cent_order) <= 3 else 6.25
    for i in range(len(cent_order)):
        for j in range(len(et_order)):
            d = delta[i, j]
            if not np.isfinite(d) or entries[i, j] == 0:
                text = "no entries"
                color = "0.35"
            else:
                display_delta = round(float(auc_with[i, j]), 3) - round(float(auc_without[i, j]), 3)
                if display_delta > 0.0:
                    winner = "cent winner"
                elif display_delta < 0.0:
                    winner = "cent loser"
                else:
                    winner = "tie"
                text = (
                    f"C {auc_with[i, j]:.3f}\n"
                    f"N {auc_without[i, j]:.3f}\n"
                    f"{display_delta:+.3f}\n{winner}"
                )
                color = "white"
            ax.text(j, i, text, ha="center", va="center", fontsize=fontsize, fontweight="bold", color=color, linespacing=0.96)

    ax.set_title(args.title, fontsize=14, fontweight="bold", pad=18)
    fig.text(
        0.985,
        0.965,
        r"Photon12+20 vs Jet12+20+30, $15<E_T<35$ GeV",
        ha="right",
        va="top",
        fontsize=9.6,
        color="0.32",
    )

    cbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=ax,
        pad=0.015,
        fraction=0.045,
        ticks=[0, 1, 2],
    )
    cbar.ax.set_yticklabels(["centrality lower", "tie", "centrality higher"])
    cbar.set_label("AUC comparison", fontsize=9.2)
    cbar.ax.tick_params(labelsize=8.5)

    fig.tight_layout(rect=[0.035, 0.04, 0.98, 0.90])
    png = args.outdir / f"{args.prefix}.png"
    csv_out = args.outdir / f"{args.prefix}.csv"
    fig.savefig(png)
    plt.close(fig)
    df.to_csv(csv_out, index=False)
    print(png)
    print(csv_out)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Make THE-8 Branch A ladder score-separation slide from compact validation histograms."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CENTRALITY = [("0_20", "0-20%"), ("20_50", "20-50%"), ("50_80", "50-80%")]


def binned_auc(signal_counts: np.ndarray, background_counts: np.ndarray) -> float:
    sig_total = float(np.sum(signal_counts))
    bkg_total = float(np.sum(background_counts))
    if sig_total <= 0 or bkg_total <= 0:
        return float("nan")
    bkg_below = np.cumsum(background_counts) - background_counts
    wins = float(np.sum(signal_counts * bkg_below))
    ties = float(np.sum(signal_counts * background_counts))
    return (wins + 0.5 * ties) / (sig_total * bkg_total)


def step(ax: plt.Axes, edges: np.ndarray, density: np.ndarray, color: str, label: str, linestyle: str) -> None:
    y = np.r_[density, density[-1]]
    ax.step(edges, y, where="post", color=color, lw=1.9, linestyle=linestyle, label=label)


def draw_density(ax: plt.Axes, edges: np.ndarray, density: np.ndarray, color: str, label: str, *, linestyle: str) -> None:
    y = np.r_[density, density[-1]]
    ax.fill_between(edges, y, step="post", color=color, alpha=0.12, linewidth=0)
    ax.step(edges, y, where="post", color=color, lw=2.25, linestyle=linestyle, label=label)


def fmt_millions(value: int) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f}M"
    if value >= 1_000:
        return f"{value / 1_000:.0f}k"
    return str(value)


def auc_color(auc: float) -> str:
    if not np.isfinite(auc):
        return "#64748B"
    if auc >= 0.88:
        return "#047857"
    if auc >= 0.82:
        return "#2563EB"
    return "#B45309"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--tag", default="the8_branchA_ladder_score_separation_slide23_style")
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    branches = payload["branches"]
    args.outdir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.linewidth": 1.05,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )

    fig, axes = plt.subplots(3, 3, figsize=(16, 9), sharex=True, sharey=True, dpi=160)
    fig.patch.set_facecolor("white")
    sig_color = "#1D4ED8"
    bkg_color = "#DC2626"
    ink = "#111827"
    muted = "#475569"
    light_grid = "#E5E7EB"
    row_colors = ["#E0F2FE", "#ECFDF5", "#FFF7ED"]
    summary_rows: list[dict[str, object]] = []

    ymax = 0.0
    cache: list[dict[str, object]] = []
    for branch in branches:
        edges = np.asarray(branch["bin_edges"], dtype=float)
        for key, label in CENTRALITY:
            cent = branch["by_centrality"][key]
            sig_density = np.asarray(cent["signal"]["density"], dtype=float)
            bkg_density = np.asarray(cent["background"]["density"], dtype=float)
            sig_counts = np.asarray(cent["signal"]["counts"], dtype=float)
            bkg_counts = np.asarray(cent["background"]["counts"], dtype=float)
            item = {
                "branch": branch,
                "cent_key": key,
                "cent_label": label,
                "edges": edges,
                "sig_density": sig_density,
                "bkg_density": bkg_density,
                "sig_counts": sig_counts,
                "bkg_counts": bkg_counts,
                "auc": binned_auc(sig_counts, bkg_counts),
            }
            cache.append(item)
            ymax = max(ymax, float(np.nanmax(sig_density)), float(np.nanmax(bkg_density)))
    ymax = 1.18 * ymax if ymax > 0 else 1.0

    for irow, branch in enumerate(branches):
        for icol, (_, cent_label) in enumerate(CENTRALITY):
            ax = axes[irow, icol]
            item = next(x for x in cache if x["branch"] is branch and x["cent_label"] == cent_label)
            draw_density(ax, item["edges"], item["sig_density"], sig_color, "Signal", linestyle="-")
            draw_density(ax, item["edges"], item["bkg_density"], bkg_color, "Background", linestyle="--")
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(0.0, ymax)
            ax.grid(True, color=light_grid, lw=0.55, alpha=0.85)
            ax.tick_params(labelsize=9.2, pad=2)
            ax.set_facecolor("#FFFFFF")
            for spine in ax.spines.values():
                spine.set_color("#111827")
            if irow == 0:
                ax.set_title(f"{cent_label} centrality", fontsize=14.6, fontweight="bold", pad=8, color=ink)
            if icol == 0:
                ax.set_ylabel("Area-normalized density", fontsize=10.4, color=ink)
                ax.text(
                    -0.255,
                    0.5,
                    branch["label"],
                    transform=ax.transAxes,
                    rotation=90,
                    ha="center",
                    va="center",
                    fontsize=16.0,
                    fontweight="bold",
                    color=ink,
                    bbox=dict(boxstyle="round,pad=0.26", fc=row_colors[irow], ec="#CBD5E1", lw=0.8),
                )
            if irow == 2:
                ax.set_xlabel("BDT score", fontsize=10.8, color=ink)
            signal_entries = int(np.sum(item["sig_counts"]))
            background_entries = int(np.sum(item["bkg_counts"]))
            auc = float(item["auc"])
            ax.text(
                0.965,
                0.93,
                f"AUC {auc:.3f}",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=11.0,
                fontweight="bold",
                color=ink,
                bbox=dict(boxstyle="round,pad=0.24", fc="white", ec="#CBD5E1", alpha=0.94),
            )
            ax.text(
                0.965,
                0.775,
                f"S {fmt_millions(signal_entries)}   B {fmt_millions(background_entries)}",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=9.0,
                color=muted,
            )
            summary_rows.append(
                {
                    "branch": branch["label"],
                    "centrality": cent_label,
                    "auc_binned": item["auc"],
                    "signal_entries": signal_entries,
                    "background_entries": background_entries,
                    "inclusive_auc": branch["summary"].get("globalEtCent1535_bdt_noIso_auc", ""),
                    "finite_score_fraction": branch["summary"].get("finite_score_fraction", ""),
                    "total_entries": branch["summary"].get("total_entries", ""),
                    "scored_entries": branch["summary"].get("scored_entries", ""),
                }
            )

    jet12_auc = float(branches[0]["summary"].get("globalEtCent1535_bdt_noIso_auc", "nan"))
    jet123_auc = float(branches[1]["summary"].get("globalEtCent1535_bdt_noIso_auc", "nan"))
    jet1234_auc = float(branches[2]["summary"].get("globalEtCent1535_bdt_noIso_auc", "nan"))
    fig.text(
        0.035,
        0.980,
        "Jet30 and Jet40 strengthen score separation",
        ha="left",
        va="top",
        fontsize=28.0,
        fontweight="bold",
        color=ink,
    )
    fig.text(
        0.035,
        0.922,
        r"Embedded background ladder comparison; full-stat validation, $15 < E_T < 35$ GeV",
        ha="left",
        va="top",
        fontsize=15.0,
        color=muted,
    )
    bullet_font = 13.8
    fig.text(
        0.035,
        0.885,
        "• Each row scores the same validation slice while adding harder embedded-jet background samples.",
        ha="left",
        va="top",
        fontsize=bullet_font,
        color="#334155",
    )
    fig.text(
        0.035,
        0.858,
        "• AUC gain alone is a diagnostic; fake rate, composition, and ABCD closure still decide the default.",
        ha="left",
        va="top",
        fontsize=bullet_font,
        color="#334155",
    )
    fig.text(
        0.738,
        0.976,
        "Inclusive AUC",
        ha="left",
        va="top",
        fontsize=10.6,
        color=muted,
    )
    auc_text = f"{jet12_auc:.3f}  ->  {jet123_auc:.3f}  ->  {jet1234_auc:.3f}"
    fig.text(
        0.738,
        0.944,
        auc_text,
        ha="left",
        va="top",
        fontsize=22.0,
        fontweight="bold",
        color="#047857",
    )
    sphenix_x = 0.060
    sphenix_y = 0.795
    fig.text(
        sphenix_x,
        sphenix_y,
        "sPHENIX",
        ha="left",
        va="bottom",
        fontsize=14.8,
        fontstyle="italic",
        fontweight="bold",
        color=ink,
    )
    fig.text(
        sphenix_x + 0.078,
        sphenix_y,
        "Internal  Au+Au embedded validation",
        ha="left",
        va="bottom",
        fontsize=14.2,
        color=ink,
    )
    handles, labels = axes[0, 1].get_legend_handles_labels()
    leg = fig.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.738, 0.905),
        ncol=2,
        frameon=True,
        fontsize=12.8,
        handlelength=3.1,
        columnspacing=1.35,
        borderpad=0.48,
        labelspacing=0.45,
    )
    leg.get_frame().set_facecolor("white")
    leg.get_frame().set_edgecolor("#CBD5E1")
    leg.get_frame().set_linewidth(1.2)
    leg.get_frame().set_alpha(0.98)
    for text in leg.get_texts():
        text.set_color(ink)
    fig.text(
        0.035,
        0.030,
        "Panel AUCs are computed from compact binned score histograms; final default choice still needs fake-rate, composition, and ABCD-closure checks.",
        ha="left",
        va="bottom",
        fontsize=12.2,
        color="#334155",
    )
    fig.tight_layout(rect=[0.060, 0.095, 0.985, 0.785], h_pad=0.72, w_pad=1.42)

    png = args.outdir / f"{args.tag}.png"
    csv_path = args.outdir / f"{args.tag}.csv"
    json_path = args.outdir / f"{args.tag}.json"
    fig.savefig(png)
    plt.close(fig)

    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    json_path.write_text(
        json.dumps(
            {
                "schema": "THE8_BRANCH_A_LADDER_SCORE_SEPARATION_PANEL_V1",
                "input": str(args.input),
                "png": str(png),
                "summary_csv": str(csv_path),
                "rows": summary_rows,
                "note": "AUC values in the panel are computed from compact binned histograms; inclusive AUC values are copied from validation_summary.txt.",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(png)
    print(csv_path)
    print(json_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

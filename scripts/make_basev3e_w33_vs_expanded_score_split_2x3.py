#!/usr/bin/env python3
"""Make 2x3 score separation plot for compact base-v3E+w33 vs expanded BDT."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


TOP_PRODUCT = "baseBDT_v3E_withCentrality_w33"
EXPANDED_PRODUCT = "globalEtCent1535_bdt_noIso"
CENT_BINS = [
    ("0-20", "0_20", ["0-10", "10-20"]),
    ("20-50", "20_50", ["20-30", "30-40", "40-50"]),
    ("50-80", "50_80", ["50-60", "60-80"]),
]
ROWS = [
    ("compact_w33", "Base PPG12 baseV3E\n+ centrality + weta33/wphi33"),
    ("expanded", "Expanded feature-list\nglobal BDT"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compact-hist-csv", type=Path, required=True)
    parser.add_argument("--expanded-fine-hist-csv", type=Path, required=True)
    parser.add_argument("--expanded-metrics-json", type=Path, required=True)
    parser.add_argument("--auc-reference-csv", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--prefix", default="basev3e_w33_vs_expanded_score_split_2x3")
    return parser.parse_args()


def step(ax, edges: np.ndarray, density: np.ndarray, color: str, label: str) -> None:
    y = np.r_[density, density[-1]]
    ax.step(edges, y, where="post", lw=2.0, color=color, label=label)
    ax.fill_between(edges, y, step="post", color=color, alpha=0.08, linewidth=0)


def densities_from_counts(edges: np.ndarray, sig: np.ndarray, bkg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    widths = np.diff(edges)
    sig_total = float(sig.sum())
    bkg_total = float(bkg.sum())
    sig_density = np.divide(sig, sig_total * widths, out=np.zeros_like(sig, dtype=float), where=(sig_total * widths) > 0)
    bkg_density = np.divide(bkg, bkg_total * widths, out=np.zeros_like(bkg, dtype=float), where=(bkg_total * widths) > 0)
    return sig_density, bkg_density


def load_auc_reference(path: Path) -> dict[tuple[str, str], float]:
    df = pd.read_csv(path)
    wanted = {
        "compact_w33": "base v3E + w33 + centrality",
        "expanded": "32-feature global",
    }
    out: dict[tuple[str, str], float] = {}
    for row_key, model in wanted.items():
        for cent_label, cent_key, _ in CENT_BINS:
            part = df[(df["model"] == model) & (df["centrality_bin"] == cent_key)]
            if part.empty:
                raise SystemExit(f"Missing AUC reference row for {model} {cent_key}")
            out[(row_key, cent_label)] = float(part["auc_entries_weighted"].iloc[0])
    return out


def load_compact(path: Path, auc_reference: dict[tuple[str, str], float]) -> dict[tuple[str, str], dict[str, object]]:
    df = pd.read_csv(path)
    out: dict[tuple[str, str], dict[str, object]] = {}
    for cent_label, _, _ in CENT_BINS:
        part = df[(df["product"] == TOP_PRODUCT) & (df["centrality_bin"] == cent_label)].sort_values("score_low")
        if part.empty:
            raise SystemExit(f"Missing compact rows for {TOP_PRODUCT} {cent_label}")
        edges = np.r_[part["score_low"].to_numpy(float), float(part["score_high"].iloc[-1])]
        sig_density, bkg_density = densities_from_counts(
            edges,
            part["signal_count"].to_numpy(float),
            part["background_count"].to_numpy(float),
        )
        out[("compact_w33", cent_label)] = {
            "edges": edges,
            "signal_density": sig_density,
            "background_density": bkg_density,
            "auc": float(auc_reference[("compact_w33", cent_label)]),
            "signal_count": part["signal_count"].to_numpy(float),
            "background_count": part["background_count"].to_numpy(float),
        }
    return out


def load_expanded(hist_path: Path, metrics_path: Path, auc_reference: dict[tuple[str, str], float]) -> dict[tuple[str, str], dict[str, object]]:
    df = pd.read_csv(hist_path)
    df = df[(df["product"] == EXPANDED_PRODUCT) & (df["et_low"] >= 15) & (df["et_high"] <= 35)].copy()
    json.loads(metrics_path.read_text())

    out: dict[tuple[str, str], dict[str, object]] = {}
    for cent_label, cent_key, fine_cent_bins in CENT_BINS:
        part = df[df["centrality_bin"].isin(fine_cent_bins)].copy()
        if part.empty:
            raise SystemExit(f"Missing expanded rows for {EXPANDED_PRODUCT} {cent_label}")
        grouped = (
            part.groupby(["score_low", "score_high"], as_index=False)[["signal_count", "background_count"]]
            .sum()
            .sort_values("score_low")
        )
        edges = np.r_[grouped["score_low"].to_numpy(float), float(grouped["score_high"].iloc[-1])]
        sig_density, bkg_density = densities_from_counts(
            edges,
            grouped["signal_count"].to_numpy(float),
            grouped["background_count"].to_numpy(float),
        )
        out[("expanded", cent_label)] = {
            "edges": edges,
            "signal_density": sig_density,
            "background_density": bkg_density,
            "auc": float(auc_reference[("expanded", cent_label)]),
            "signal_count": grouped["signal_count"].to_numpy(float),
            "background_count": grouped["background_count"].to_numpy(float),
        }
    return out


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    auc_reference = load_auc_reference(args.auc_reference_csv)
    prepared = {}
    prepared.update(load_compact(args.compact_hist_csv, auc_reference))
    prepared.update(load_expanded(args.expanded_fine_hist_csv, args.expanded_metrics_json, auc_reference))

    max_y = max(
        float(np.nanmax(cell["signal_density"])) if len(cell["signal_density"]) else 0.0
        for cell in prepared.values()
    )
    max_y = max(
        max_y,
        max(float(np.nanmax(cell["background_density"])) if len(cell["background_density"]) else 0.0 for cell in prepared.values()),
    )

    sig_color = "#1f77b4"
    bkg_color = "#ff7f0e"
    fig, axes = plt.subplots(2, 3, figsize=(13.8, 6.75), dpi=240, sharex=True, sharey=True)

    for row_i, (row_key, row_label) in enumerate(ROWS):
        for col_i, (cent_label, _, _) in enumerate(CENT_BINS):
            ax = axes[row_i, col_i]
            cell = prepared[(row_key, cent_label)]
            step(ax, cell["edges"], cell["signal_density"], sig_color, "Signal")
            step(ax, cell["edges"], cell["background_density"], bkg_color, "Background")
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(0.0, max_y * 1.16)
            ax.grid(True, color="0.90", linewidth=0.55)
            ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=10.5)
            for spine in ax.spines.values():
                spine.set_linewidth(1.0)
            if row_i == 0:
                ax.set_title(f"{cent_label}%", fontsize=15, fontweight="bold", pad=8)
            ax.text(
                0.055,
                0.90,
                f"$E_T$-wtd AUC {cell['auc']:.3f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=10.4,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="0.82", alpha=0.92),
            )
            if col_i == 0:
                ax.text(
                    0.055,
                    0.60,
                    row_label,
                    transform=ax.transAxes,
                    ha="left",
                    va="center",
                    fontsize=10.4,
                    fontweight="bold",
                    color="0.10",
                    bbox=dict(boxstyle="round,pad=0.26", fc="white", ec="0.78", alpha=0.94),
                )

    axes[0, 2].legend(loc="upper right", frameon=False, fontsize=12)
    for ax in axes[-1, :]:
        ax.set_xlabel("BDT score", fontsize=13)
    fig.text(0.025, 0.49, "Area-normalized density", rotation=90, ha="center", va="center", fontsize=12.8)
    fig.suptitle("BDT score separation: compact 3x3-width inputs vs expanded feature list", fontsize=16.4, fontweight="bold", y=0.978)
    fig.text(
        0.52,
        0.931,
        r"Photon12+20 signal vs Jet12+20+30 background, $15<E_T<35$ GeV; AUC labels use entries-weighted $E_T$ bins",
        ha="center",
        va="top",
        fontsize=10.8,
        color="0.32",
    )
    fig.tight_layout(rect=[0.055, 0.06, 0.995, 0.895], w_pad=1.0, h_pad=0.95)

    png = args.outdir / f"{args.prefix}.png"
    summary_rows = []
    for row_key, row_label in ROWS:
        for cent_label, _, _ in CENT_BINS:
            cell = prepared[(row_key, cent_label)]
            summary_rows.append(
                {
                    "row": row_key,
                    "label": row_label.replace("\n", " "),
                    "centrality_bin": cent_label,
                    "auc": cell["auc"],
                    "signal_entries": int(np.sum(cell["signal_count"])),
                    "background_entries": int(np.sum(cell["background_count"])),
                }
            )
    csv_out = args.outdir / f"{args.prefix}.csv"
    fig.savefig(png)
    plt.close(fig)
    pd.DataFrame(summary_rows).to_csv(csv_out, index=False)
    print(png)
    print(csv_out)


if __name__ == "__main__":
    main()

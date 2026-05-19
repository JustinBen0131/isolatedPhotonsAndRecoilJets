#!/usr/bin/env python3
"""Make 1-row centrality tables of baseline BDT score separation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


PRODUCT = "globalEtCent1535_bdt_noIso"
FINE_BINS = [
    ("0-10", 0, 10),
    ("10-20", 10, 20),
    ("20-30", 20, 30),
    ("30-40", 30, 40),
    ("40-50", 40, 50),
    ("50-60", 50, 60),
    ("60-80", 60, 80),
]
COARSE_BINS = [
    ("0-20", 0, 20),
    ("20-50", 20, 50),
    ("50-80", 50, 80),
]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hist-csv", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--product", default=PRODUCT)
    ap.add_argument("--model-label", default="Global base-v3E + centrality model")
    ap.add_argument("--title-prefix", default="Baseline BDT score separation")
    ap.add_argument("--auc-csv", type=Path)
    ap.add_argument("--auc-model")
    return ap.parse_args()


def auc_from_binned_counts(signal: np.ndarray, background: np.ndarray) -> float:
    signal = np.asarray(signal, dtype=np.float64)
    background = np.asarray(background, dtype=np.float64)
    n_sig = float(signal.sum())
    n_bkg = float(background.sum())
    if n_sig <= 0 or n_bkg <= 0:
        return float("nan")
    bkg_below = np.r_[0.0, np.cumsum(background)[:-1]]
    favorable = float(np.sum(signal * bkg_below))
    ties = 0.5 * float(np.sum(signal * background))
    return (favorable + ties) / (n_sig * n_bkg)


def aggregate(df: pd.DataFrame, product: str, bins: list[tuple[str, float, float]]) -> list[dict[str, object]]:
    rows = []
    base = df[
        (df["product"] == product)
        & (df["et_low"] >= 15.0)
        & (df["et_high"] <= 35.0)
    ].copy()
    if base.empty:
        raise SystemExit(f"No rows found for {product} in {df}")

    for label, lo, hi in bins:
        exact = base[(base["centrality_low"] == lo) & (base["centrality_high"] == hi)].copy()
        part = exact
        if part.empty:
            part = base[(base["centrality_low"] >= lo) & (base["centrality_high"] <= hi)].copy()
        if part.empty:
            raise SystemExit(f"No histogram rows for centrality {label}")
        grouped = (
            part.groupby(["score_low", "score_high"], as_index=False)[["signal_count", "background_count"]]
            .sum()
            .sort_values("score_low")
        )
        widths = grouped["score_high"].to_numpy(dtype=float) - grouped["score_low"].to_numpy(dtype=float)
        sig = grouped["signal_count"].to_numpy(dtype=float)
        bkg = grouped["background_count"].to_numpy(dtype=float)
        sig_total = float(sig.sum())
        bkg_total = float(bkg.sum())
        sig_density = np.divide(sig, sig_total * widths, out=np.zeros_like(sig), where=(sig_total * widths) > 0)
        bkg_density = np.divide(bkg, bkg_total * widths, out=np.zeros_like(bkg), where=(bkg_total * widths) > 0)
        rows.append(
            {
                "label": label,
                "lo": lo,
                "hi": hi,
                "edges": np.r_[grouped["score_low"].to_numpy(dtype=float), grouped["score_high"].iloc[-1]],
                "signal_density": sig_density,
                "background_density": bkg_density,
                "signal_count": sig,
                "background_count": bkg,
                "signal_entries": int(round(sig_total)),
                "background_entries": int(round(bkg_total)),
                "auc": auc_from_binned_counts(sig, bkg),
            }
        )
    return rows


def load_auc_overrides(path: Path | None, model: str | None) -> dict[str, float]:
    if path is None or model is None:
        return {}
    table = pd.read_csv(path)
    if not {"model", "centrality_bin", "auc_entries_weighted"}.issubset(table.columns):
        return {}
    out: dict[str, float] = {}
    for row in table[table["model"] == model].to_dict("records"):
        key = str(row["centrality_bin"]).replace("_", "-")
        out[key] = float(row["auc_entries_weighted"])
    return out


def step(ax, edges: np.ndarray, density: np.ndarray, color: str, label: str) -> None:
    y = np.r_[density, density[-1]]
    ax.step(edges, y, where="post", lw=1.8, color=color, label=label)
    ax.fill_between(edges, y, step="post", color=color, alpha=0.08, linewidth=0)


def draw(rows: list[dict[str, object]], output: Path, title: str, subtitle: str, auc_overrides: dict[str, float]) -> None:
    n = len(rows)
    figsize = (21.0, 3.9) if n == 7 else (12.8, 3.9)
    title_size = 14.8 if n == 7 else 16.0
    subtitle_size = 9.8 if n == 7 else 10.8
    ylabel_size = 10.8 if n == 7 else 12.0
    fig, axes = plt.subplots(1, n, figsize=figsize, dpi=220, sharex=True)
    if n == 1:
        axes = [axes]

    sig_color = "#1f77b4"
    bkg_color = "#ff7f0e"
    max_y = 0.0
    for row in rows:
        max_y = max(max_y, float(np.nanmax(row["signal_density"])), float(np.nanmax(row["background_density"])))

    for ax, row in zip(axes, rows, strict=True):
        edges = row["edges"]
        step(ax, edges, row["signal_density"], sig_color, "Signal")
        step(ax, edges, row["background_density"], bkg_color, "Background")
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, max_y * 1.18)
        ax.set_title(f"{row['label']}%", fontsize=13.0, fontweight="bold", pad=7)
        auc = auc_overrides.get(str(row["label"]), row["auc"])
        ax.text(
            0.05,
            0.90,
            f"AUC {auc:.3f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=10.5,
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.20", fc="white", ec="0.82", alpha=0.92),
        )
        ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=9.5)
        ax.grid(True, color="0.90", linewidth=0.55)
        for spine in ax.spines.values():
            spine.set_linewidth(1.0)

    axes[0].set_ylabel("Area-normalized density", fontsize=ylabel_size)
    for ax in axes:
        ax.set_xlabel("BDT score", fontsize=11)
    axes[-1].legend(loc="upper right", frameon=False, fontsize=10.5)

    fig.text(0.022, 0.965, r"$\bf{\it{sPHENIX}}$ Internal", ha="left", va="top", fontsize=13.5)
    fig.text(0.50, 0.965, title, ha="center", va="top", fontsize=title_size, fontweight="bold")
    fig.text(0.50, 0.915, subtitle, ha="center", va="top", fontsize=subtitle_size, color="0.30")
    fig.tight_layout(rect=[0.018, 0.030, 0.995, 0.850], w_pad=1.0)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def write_summary(rows: list[dict[str, object]], output: Path, auc_overrides: dict[str, float]) -> None:
    fields = ["centrality_bin", "auc_from_display", "auc_from_pooled_histogram", "signal_entries", "background_entries"]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "centrality_bin": row["label"],
                    "auc_from_display": f"{auc_overrides.get(str(row['label']), row['auc']):.8g}",
                    "auc_from_pooled_histogram": f"{row['auc']:.8g}",
                    "signal_entries": row["signal_entries"],
                    "background_entries": row["background_entries"],
                }
            )


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.hist_csv)
    fine = aggregate(df, args.product, FINE_BINS)
    coarse = aggregate(df, args.product, COARSE_BINS)
    auc_overrides = load_auc_overrides(args.auc_csv, args.auc_model)
    args.outdir.mkdir(parents=True, exist_ok=True)

    safe_product = args.product.replace("/", "_")
    fine_png = args.outdir / f"{safe_product}_score_split_fine7_centrality_1x7.png"
    coarse_png = args.outdir / f"{safe_product}_score_split_coarse3_centrality_1x3.png"
    draw(
        fine,
        fine_png,
        f"{args.title_prefix} by fine centrality",
        rf"{args.model_label}; Photon12+20 signal vs Jet12+20+30 embedded inclusive background, $15 < E_T < 35$ GeV.",
        auc_overrides,
    )
    draw(
        coarse,
        coarse_png,
        f"{args.title_prefix} by coarse centrality",
        r"Same validated BDT, integrated over $15 < E_T < 35$ GeV; curves are area-normalized.",
        auc_overrides,
    )
    write_summary(fine, args.outdir / f"{safe_product}_score_split_fine7_centrality_1x7.csv", auc_overrides)
    write_summary(coarse, args.outdir / f"{safe_product}_score_split_coarse3_centrality_1x3.csv", auc_overrides)
    print(fine_png)
    print(coarse_png)


if __name__ == "__main__":
    main()

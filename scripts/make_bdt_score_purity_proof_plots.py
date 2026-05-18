#!/usr/bin/env python3
"""Make slide-ready plots proving how BDT score tracks truth purity.

The inputs are the validated sixpack BDT score histograms.  The plots are
deliberately diagnostic: they show that the score is a signal-likeness ranker
whose thresholds map to measured truth purity in simulation, while also showing
that the raw score is not itself a calibrated probability.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


PRODUCT_LABELS = {
    "globalEtCent1535_bdt_noIso": "Baseline BDT",
    "globalEtCent1535_bdt_iso": "Isolation-input BDT",
}

PRODUCT_COLORS = {
    "globalEtCent1535_bdt_noIso": "#2f6fbb",
    "globalEtCent1535_bdt_iso": "#b23a8f",
}

CENT_ORDER = ["0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-80"]
COARSE_CENT_ORDER = ["0-20", "20-50", "50-80"]
ET_ORDER = ["15-17", "17-19", "19-21", "21-23", "23-25", "25-27", "27-30", "30-35"]


def _coarse_cent(label: str) -> str:
    low = int(str(label).split("-")[0])
    if low < 20:
        return "0-20"
    if low < 50:
        return "20-50"
    return "50-80"


def _safe_div(num: np.ndarray | float, den: np.ndarray | float) -> np.ndarray | float:
    return np.divide(num, den, out=np.full_like(np.asarray(num, dtype=float), np.nan), where=np.asarray(den) != 0)


def _wilson_err(k: np.ndarray, n: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return approximate 68% Wilson lower/upper errors for binomial fractions."""
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)
    z = 1.0
    p = _safe_div(k, n)
    good = n > 0
    denom = np.full_like(n, np.nan, dtype=float)
    center = np.full_like(n, np.nan, dtype=float)
    half = np.full_like(n, np.nan, dtype=float)
    denom[good] = 1.0 + z * z / n[good]
    center[good] = (p[good] + z * z / (2.0 * n[good])) / denom[good]
    half[good] = z * np.sqrt((p[good] * (1.0 - p[good]) + z * z / (4.0 * n[good])) / n[good]) / denom[good]
    lo = np.clip(center - half, 0.0, 1.0)
    hi = np.clip(center + half, 0.0, 1.0)
    return p - lo, hi - p


def _weighted_corr(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    if finite.sum() < 3:
        return np.nan
    x = x[finite]
    y = y[finite]
    w = w[finite]
    wx = np.average(x, weights=w)
    wy = np.average(y, weights=w)
    cov = np.average((x - wx) * (y - wy), weights=w)
    vx = np.average((x - wx) ** 2, weights=w)
    vy = np.average((y - wy) ** 2, weights=w)
    if vx <= 0 or vy <= 0:
        return np.nan
    return float(cov / np.sqrt(vx * vy))


def _aggregate_score_bins(hist: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    cols = group_cols + ["score_low", "score_high"]
    out = (
        hist.groupby(cols, as_index=False)[["signal_count", "background_count"]]
        .sum()
        .sort_values(cols)
    )
    out["score_center"] = 0.5 * (out["score_low"] + out["score_high"])
    out["total_count"] = out["signal_count"] + out["background_count"]
    out["truth_purity_bin"] = out["signal_count"] / out["total_count"].replace(0, np.nan)
    return out


def _add_cumulative_above(score_bins: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    pieces = []
    for _, group in score_bins.groupby(group_cols, dropna=False):
        g = group.sort_values("score_low").copy()
        sig = g["signal_count"].to_numpy(dtype=float)
        bkg = g["background_count"].to_numpy(dtype=float)
        g["signal_above"] = sig[::-1].cumsum()[::-1]
        g["background_above"] = bkg[::-1].cumsum()[::-1]
        g["total_above"] = g["signal_above"] + g["background_above"]
        g["truth_purity_above"] = g["signal_above"] / g["total_above"].replace(0, np.nan)
        g["signal_eff_above"] = g["signal_above"] / sig.sum() if sig.sum() > 0 else np.nan
        g["background_fake_above"] = g["background_above"] / bkg.sum() if bkg.sum() > 0 else np.nan
        pieces.append(g)
    return pd.concat(pieces, ignore_index=True)


def _style_axis(ax):
    ax.grid(True, color="#e5e5e5", linewidth=0.8)
    ax.tick_params(direction="in", top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)


def plot_global_calibration(score_global: pd.DataFrame, out: Path, thresholds: dict[str, float]):
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2), dpi=180)
    ax0, ax1 = axes
    xdiag = np.linspace(0, 1, 200)
    ax0.plot(xdiag, xdiag, color="#7a7a7a", linestyle="--", linewidth=1.2, label="score = purity")

    summary = {}
    ax1b = ax1.twinx()
    # Draw the twin axis only once, then add both products to it.
    ax1b.set_ylabel("Signal efficiency kept", fontsize=11)
    ax1b.set_ylim(0, 1.03)
    ax1b.tick_params(direction="in", right=True)

    for product, label in PRODUCT_LABELS.items():
        color = PRODUCT_COLORS[product]
        g = score_global[score_global["product"] == product].sort_values("score_center")
        x = g["score_center"].to_numpy()
        purity = g["truth_purity_bin"].to_numpy()
        total = g["total_count"].to_numpy(dtype=float)
        sig = g["signal_count"].to_numpy(dtype=float)
        yerr_lo, yerr_hi = _wilson_err(sig, total)
        ax0.plot(x, purity, color=color, linewidth=2.4, label=label)
        ax0.fill_between(x, purity - yerr_lo, purity + yerr_hi, color=color, alpha=0.12, linewidth=0)
        ax1.plot(g["score_low"], g["truth_purity_above"], color=color, linewidth=2.6, label=label)
        ax1b.plot(g["score_low"], g["signal_eff_above"], color=color, linewidth=1.8, linestyle=":", alpha=0.85)
        threshold = thresholds.get(product, np.nan)
        if np.isfinite(threshold):
            ax1.axvline(threshold, color=color, linewidth=1.4, linestyle="--", alpha=0.8)
            near = g.iloc[(g["score_low"] - threshold).abs().argsort()[:1]]
            if len(near):
                summary[product] = {
                    "label": label,
                    "wp80_threshold": float(threshold),
                    "wp80_truth_purity": float(near["truth_purity_above"].iloc[0]),
                    "wp80_signal_efficiency": float(near["signal_eff_above"].iloc[0]),
                    "global_raw_truth_purity": float(g["signal_count"].sum() / g["total_count"].sum()),
                    "score_purity_weighted_corr": _weighted_corr(
                        g["score_center"].to_numpy(),
                        g["truth_purity_bin"].to_numpy(),
                        g["total_count"].to_numpy(),
                    ),
                }

    ax0.set_title("Measured truth purity in narrow score slices", fontsize=14, fontweight="bold")
    ax0.set_xlabel("BDT score")
    ax0.set_ylabel("Truth purity in score slice: S/(S+B)")
    ax0.set_xlim(0, 1)
    ax0.set_ylim(0, 1.02)
    ax0.legend(loc="lower right", frameon=False, fontsize=10)
    _style_axis(ax0)

    ax1.set_title("Purity of candidates kept above a score threshold", fontsize=14, fontweight="bold")
    ax1.set_xlabel("Minimum BDT score threshold")
    ax1.set_ylabel("Truth purity after threshold")
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1.02)
    ax1.legend(loc="lower right", frameon=False, fontsize=10)
    _style_axis(ax1)

    fig.suptitle("BDT score is a signal-likeness ranker, not a calibrated probability", fontsize=16, fontweight="bold", y=0.98)
    fig.text(
        0.5,
        0.01,
        "Dotted lines on the right show signal efficiency kept; dashed vertical lines mark each model's WP80 threshold.",
        ha="center",
        fontsize=10,
        color="#444444",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    fig.savefig(out / "bdt_score_truth_purity_global_calibration.png")
    plt.close(fig)
    return summary


def plot_centrality_ladder(score_coarse: pd.DataFrame, out: Path, thresholds: dict[str, float]):
    fig, axes = plt.subplots(2, 3, figsize=(14.2, 7.4), dpi=180, sharex=True, sharey=True)
    for row, (product, label) in enumerate(PRODUCT_LABELS.items()):
        color = PRODUCT_COLORS[product]
        for col, cent in enumerate(COARSE_CENT_ORDER):
            ax = axes[row, col]
            g = score_coarse[(score_coarse["product"] == product) & (score_coarse["coarse_centrality"] == cent)]
            g = g.sort_values("score_low")
            raw = g["signal_count"].sum() / max(g["total_count"].sum(), 1)
            ax.plot(g["score_low"], g["truth_purity_above"], color=color, linewidth=2.7)
            ax.axhline(raw, color="#666666", linewidth=1.2, linestyle=":", label="raw purity" if (row, col) == (0, 0) else None)
            threshold = thresholds.get(product, np.nan)
            if np.isfinite(threshold):
                ax.axvline(threshold, color=color, linewidth=1.4, linestyle="--", alpha=0.85)
            ax.text(
                0.04,
                0.92,
                f"{label}\nCentrality {cent}%\nraw={raw:.3f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=9.5,
                bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#d4d4d4", alpha=0.92),
            )
            if row == 1:
                ax.set_xlabel("Minimum BDT score")
            if col == 0:
                ax.set_ylabel("Truth purity after threshold")
            ax.set_xlim(0, 1)
            ax.set_ylim(0.45, 1.02)
            _style_axis(ax)
    fig.suptitle("Thresholding the BDT score raises measured truth purity in every coarse centrality region", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out / "bdt_score_truth_purity_threshold_ladder_by_centrality.png")
    plt.close(fig)


def _pivot_fine(summary: pd.DataFrame, product: str, column: str) -> pd.DataFrame:
    d = summary[summary["product"] == product].copy()
    d["centrality_bin"] = pd.Categorical(d["centrality_bin"], CENT_ORDER, ordered=True)
    d["et_bin"] = pd.Categorical(d["et_bin"], ET_ORDER, ordered=True)
    return d.pivot_table(index="centrality_bin", columns="et_bin", values=column, observed=False).reindex(index=CENT_ORDER, columns=ET_ORDER)


def _annotate_heatmap(ax, data: pd.DataFrame, fmt: str, cmap_name: str, vmin=None, vmax=None, center=None):
    arr = data.to_numpy(dtype=float)
    if center is not None:
        norm = TwoSlopeNorm(vmin=vmin, vcenter=center, vmax=vmax)
        im = ax.imshow(arr, cmap=cmap_name, norm=norm, aspect="auto")
    else:
        im = ax.imshow(arr, cmap=cmap_name, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(ET_ORDER)), ET_ORDER, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(len(CENT_ORDER)), [f"{x}%" for x in CENT_ORDER], fontsize=8)
    finite = arr[np.isfinite(arr)]
    midpoint = np.nanmedian(finite) if finite.size else 0.5
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            val = arr[i, j]
            if np.isfinite(val):
                color = "white" if val > midpoint else "black"
                ax.text(j, i, format(val, fmt), ha="center", va="center", fontsize=7.2, color=color)
    return im


def plot_fine_bin_purity_proof(score_fine: pd.DataFrame, out: Path, thresholds: dict[str, float]):
    rows = []
    for (product, cent, et), g in score_fine.groupby(["product", "centrality_bin", "et_bin"]):
        g = g.sort_values("score_low")
        raw_purity = g["signal_count"].sum() / max(g["total_count"].sum(), 1)
        threshold = thresholds.get(product, np.nan)
        if np.isfinite(threshold):
            gthr = g.iloc[(g["score_low"] - threshold).abs().argsort()[:1]]
            wp80_purity = float(gthr["truth_purity_above"].iloc[0])
            wp80_eff = float(gthr["signal_eff_above"].iloc[0])
            wp80_fake = float(gthr["background_fake_above"].iloc[0])
        else:
            wp80_purity = np.nan
            wp80_eff = np.nan
            wp80_fake = np.nan
        high = g[g["score_low"] >= 0.75]
        low = g[g["score_high"] <= 0.25]
        high_purity = high["signal_count"].sum() / max(high["total_count"].sum(), 1)
        low_purity = low["signal_count"].sum() / max(low["total_count"].sum(), 1)
        corr = _weighted_corr(
            g["score_center"].to_numpy(),
            g["truth_purity_bin"].to_numpy(),
            g["total_count"].to_numpy(),
        )
        rows.append(
            {
                "product": product,
                "model_label": PRODUCT_LABELS.get(product, product),
                "centrality_bin": cent,
                "et_bin": et,
                "raw_truth_purity": raw_purity,
                "wp80_truth_purity": wp80_purity,
                "wp80_signal_efficiency": wp80_eff,
                "wp80_background_fake_rate": wp80_fake,
                "purity_gain_abs": wp80_purity - raw_purity,
                "low_score_truth_purity_score_lt_025": low_purity,
                "high_score_truth_purity_score_ge_075": high_purity,
                "high_minus_low_score_purity": high_purity - low_purity,
                "score_vs_truth_purity_weighted_corr": corr,
                "signal_entries": float(g["signal_count"].sum()),
                "background_entries": float(g["background_count"].sum()),
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "bdt_score_truth_purity_fine7x8_proof_summary.csv", index=False)

    fig, axes = plt.subplots(2, 3, figsize=(16.0, 8.1), dpi=180)
    specs = [
        ("wp80_truth_purity", "WP80 truth purity", "viridis", 0.75, 1.0, None, ".2f"),
        ("purity_gain_abs", "WP80 gain over raw", "magma", 0.0, 0.35, None, ".2f"),
        ("high_minus_low_score_purity", "High-score minus low-score purity", "coolwarm", -0.2, 1.0, 0.0, ".2f"),
    ]
    for row, (product, label) in enumerate(PRODUCT_LABELS.items()):
        for col, (column, title, cmap, vmin, vmax, center, fmt) in enumerate(specs):
            ax = axes[row, col]
            data = _pivot_fine(summary, product, column)
            im = _annotate_heatmap(ax, data, fmt, cmap, vmin, vmax, center)
            ax.set_title(f"{label}: {title}", fontsize=10.8, fontweight="bold")
            if col == 0:
                ax.set_ylabel("Centrality bin")
            if row == 1:
                ax.set_xlabel("$E_T$ bin [GeV]")
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
            cb.ax.tick_params(labelsize=8)
    fig.suptitle("Fine-bin proof that higher BDT score selects higher-purity simulated photons", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out / "bdt_score_truth_purity_fine7x8_proof_maps.png")
    plt.close(fig)
    return summary


def plot_score_to_purity_lookup(score_global: pd.DataFrame, out: Path):
    fig, ax = plt.subplots(figsize=(8.8, 5.7), dpi=180)
    lookup_rows = []
    for product, label in PRODUCT_LABELS.items():
        color = PRODUCT_COLORS[product]
        g = score_global[score_global["product"] == product].sort_values("score_center")
        ax.plot(g["score_low"], g["truth_purity_above"], color=color, linewidth=2.8, label=label)
        for target in [0.90, 0.95, 0.97]:
            gg = g[g["truth_purity_above"] >= target]
            if len(gg):
                x = float(gg["score_low"].iloc[0])
                ax.scatter([x], [target], color=color, s=34, zorder=4)
                lookup_rows.append((target, label, x))
    for y in [0.90, 0.95, 0.97]:
        ax.axhline(y, color="#c7c7c7", linewidth=1.0, linestyle="--", zorder=0)
    ax.set_title("Empirical score-to-purity lookup from truth-labeled simulation", fontsize=15, fontweight="bold")
    ax.set_xlabel("Minimum BDT score threshold")
    ax.set_ylabel("Measured truth purity after threshold")
    ax.set_xlim(0, 1)
    ax.set_ylim(0.72, 1.01)
    ax.legend(loc="lower right", frameon=False)
    _style_axis(ax)

    by_target = {target: {} for target in [0.90, 0.95, 0.97]}
    for target, label, x in lookup_rows:
        by_target[target][label] = x
    table_lines = ["threshold needed for target purity"]
    for target in [0.90, 0.95, 0.97]:
        base = by_target[target].get("Baseline BDT", np.nan)
        iso = by_target[target].get("Isolation-input BDT", np.nan)
        table_lines.append(f"{target:.0%}: baseline {base:.2f}, isolation-input {iso:.2f}")
    ax.text(
        0.03,
        0.06,
        "\n".join(table_lines),
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9.5,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#d5d5d5", alpha=0.95),
    )
    fig.tight_layout()
    fig.savefig(out / "bdt_score_to_truth_purity_lookup.png")
    plt.close(fig)


def plot_score_bin_composition(score_global: pd.DataFrame, out: Path):
    fig, axes = plt.subplots(2, 1, figsize=(10.8, 7.2), dpi=180, sharex=True, sharey=True)
    for ax, (product, label) in zip(axes, PRODUCT_LABELS.items()):
        g = score_global[score_global["product"] == product].sort_values("score_center")
        x = g["score_center"].to_numpy()
        width = 0.022
        sig_frac = g["truth_purity_bin"].to_numpy()
        bkg_frac = 1.0 - sig_frac
        ax.bar(x, bkg_frac, width=width, color="#d84c4c", alpha=0.82, label="background fraction")
        ax.bar(x, sig_frac, bottom=bkg_frac, width=width, color="#2f80c9", alpha=0.88, label="signal fraction")
        ax.plot(x, sig_frac, color="#0d3f75", linewidth=2.0)
        ax.text(
            0.02,
            0.90,
            label,
            transform=ax.transAxes,
            fontsize=12,
            fontweight="bold",
            ha="left",
            va="top",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#d4d4d4", alpha=0.94),
        )
        ax.set_ylabel("Fraction within score bin")
        ax.set_ylim(0, 1.02)
        _style_axis(ax)
    axes[0].legend(loc="center right", frameon=False)
    axes[-1].set_xlabel("BDT score bin")
    fig.suptitle("Truth composition by BDT score bin: high score is high purity, not perfect purity", fontsize=15.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out / "bdt_score_truth_composition_by_score_bin.png")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--validation-dir",
        type=Path,
        default=Path("dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/validation/bdt_finished_only_20260516_180817"),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/slideReady/bdt_score_purity_proof"),
    )
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    hist_path = args.validation_dir / "fine7x8_score_histograms_globalEtCent1535_bdt_iso_noIso.csv"
    wp_path = args.validation_dir / "bdt_working_points_target80.csv"
    hist = pd.read_csv(hist_path)
    hist = hist[hist["product"].isin(PRODUCT_LABELS)].copy()
    hist["coarse_centrality"] = hist["centrality_bin"].map(_coarse_cent)

    wp = pd.read_csv(wp_path)
    threshold_col = "threshold" if "threshold" in wp.columns else "intercept"
    thresholds = dict(zip(wp["product"], wp[threshold_col]))

    score_global = _add_cumulative_above(_aggregate_score_bins(hist, ["product"]), ["product"])
    score_coarse = _add_cumulative_above(_aggregate_score_bins(hist, ["product", "coarse_centrality"]), ["product", "coarse_centrality"])
    score_fine = _add_cumulative_above(_aggregate_score_bins(hist, ["product", "centrality_bin", "et_bin"]), ["product", "centrality_bin", "et_bin"])

    global_summary = plot_global_calibration(score_global, args.outdir, thresholds)
    plot_centrality_ladder(score_coarse, args.outdir, thresholds)
    fine_summary = plot_fine_bin_purity_proof(score_fine, args.outdir, thresholds)
    plot_score_to_purity_lookup(score_global, args.outdir)
    plot_score_bin_composition(score_global, args.outdir)

    with (args.outdir / "bdt_score_truth_purity_global_summary.json").open("w") as f:
        json.dump(global_summary, f, indent=2, sort_keys=True)

    print(f"[purityProof] wrote {args.outdir}")
    print(f"[purityProof] fine bins: {len(fine_summary)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

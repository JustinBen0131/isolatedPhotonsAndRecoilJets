#!/usr/bin/env python3
"""Diagnose AUC trends in the Au+Au BDT train/test split stress test."""

from __future__ import annotations
# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath
_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
_CODEX_SCRIPTS_DIR_STR = str(_CODEX_SCRIPTS_DIR)
if _CODEX_SCRIPTS_DIR_STR not in _codex_sys.path:
    _codex_sys.path.append(_CODEX_SCRIPTS_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_SCRIPTS_DIR_STR

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


INK = "#111827"
MUTED = "#4B5563"
GRID = "#E5E7EB"
BLUE = "#0072B2"
RED = "#CC334E"
ORANGE = "#E69F00"
GREEN = "#009E73"
PURPLE = "#7E57C2"
COLORS = {"90/10": BLUE, "50/50": ORANGE, "10/90": GREEN}
SPLIT_TRAIN_FRACTION = {"90/10": 0.90, "50/50": 0.50, "10/90": 0.10}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, required=True, help="split-study compact JSON")
    ap.add_argument("--outdir", type=Path, required=True)
    return ap.parse_args()


def density_moments(edges: np.ndarray, density: np.ndarray) -> dict[str, float]:
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)
    mass = float(np.sum(density * widths))
    mean = float(np.sum(centers * density * widths) / mass)
    var = float(np.sum(((centers - mean) ** 2) * density * widths) / mass)
    return {"mass": mass, "mean": mean, "std": float(np.sqrt(max(var, 0.0)))}


def density_quantile(edges: np.ndarray, density: np.ndarray, q: float) -> float:
    widths = np.diff(edges)
    probs = density * widths
    total = float(np.sum(probs))
    if total <= 0:
        return float("nan")
    cdf = np.cumsum(probs) / total
    idx = int(np.searchsorted(cdf, q, side="left"))
    idx = max(0, min(idx, len(widths) - 1))
    prev = 0.0 if idx == 0 else float(cdf[idx - 1])
    frac = 0.0 if cdf[idx] <= prev else (q - prev) / (float(cdf[idx]) - prev)
    return float(edges[idx] + frac * widths[idx])


def survival_at_threshold(edges: np.ndarray, density: np.ndarray, threshold: float) -> float:
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)
    probs = density * widths
    total = float(np.sum(probs))
    if total <= 0:
        return float("nan")
    return float(np.sum(probs[centers >= threshold]) / total)


def wp_signal_threshold(edges: np.ndarray, signal_density: np.ndarray, target_eff: float = 0.80) -> float:
    return density_quantile(edges, signal_density, 1.0 - target_eff)


def hanley_mcneil_se(auc: float, n_signal: int, n_background: int) -> float:
    q1 = auc / (2.0 - auc)
    q2 = 2.0 * auc * auc / (1.0 + auc)
    var = (
        auc * (1.0 - auc)
        + (n_signal - 1) * (q1 - auc * auc)
        + (n_background - 1) * (q2 - auc * auc)
    ) / (n_signal * n_background)
    return float(np.sqrt(max(var, 0.0)))


def approx_auc_from_hist(edges: np.ndarray, sig_density: np.ndarray, bkg_density: np.ndarray) -> float:
    widths = np.diff(edges)
    sig = sig_density * widths
    bkg = bkg_density * widths
    sig = sig / np.sum(sig)
    bkg = bkg / np.sum(bkg)
    bkg_cdf_less = np.cumsum(bkg) - bkg
    ties = 0.5 * bkg
    return float(np.sum(sig * (bkg_cdf_less + ties)))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    payload = json.loads(args.input.read_text())
    args.outdir.mkdir(parents=True, exist_ok=True)

    hist_map = {(h["split_label"], h["cent_label"]): h for h in payload["histograms"]}
    splits = [r["label"] for r in payload["rows"]]
    cents = [c["label"] for c in payload["centrality_bins"]]
    baseline = "90/10"

    rows: list[dict[str, object]] = []
    for cent in cents:
        base_auc = hist_map[(baseline, cent)]["auc"]
        base_wp_fake = None
        for split in splits:
            h = hist_map[(split, cent)]
            edges = np.asarray(h["bin_edges"], dtype=float)
            sig = np.asarray(h["signal_density"], dtype=float)
            bkg = np.asarray(h["background_density"], dtype=float)
            sig_m = density_moments(edges, sig)
            bkg_m = density_moments(edges, bkg)
            wp80_thr = wp_signal_threshold(edges, sig, 0.80)
            wp80_fake = survival_at_threshold(edges, bkg, wp80_thr)
            if split == baseline:
                base_wp_fake = wp80_fake
            overlap = float(np.sum(np.minimum(sig, bkg) * np.diff(edges)))
            sep = (sig_m["mean"] - bkg_m["mean"]) / np.sqrt(0.5 * (sig_m["std"] ** 2 + bkg_m["std"] ** 2))
            auc = float(h["auc"])
            se = hanley_mcneil_se(auc, int(h["signal_entries"]), int(h["background_entries"]))
            rows.append(
                {
                    "centrality": cent,
                    "split": split,
                    "train_fraction": SPLIT_TRAIN_FRACTION[split],
                    "auc": auc,
                    "auc_delta_vs_90_10": auc - base_auc,
                    "auc_se_hanley_mcneil": se,
                    "auc_delta_vs_90_10_over_independent_se": (auc - base_auc) / np.sqrt(2 * se * se),
                    "hist_approx_auc": approx_auc_from_hist(edges, sig, bkg),
                    "signal_entries": int(h["signal_entries"]),
                    "background_entries": int(h["background_entries"]),
                    "signal_mean_score": sig_m["mean"],
                    "background_mean_score": bkg_m["mean"],
                    "mean_score_gap": sig_m["mean"] - bkg_m["mean"],
                    "signal_q20": density_quantile(edges, sig, 0.20),
                    "signal_q50": density_quantile(edges, sig, 0.50),
                    "signal_q80": density_quantile(edges, sig, 0.80),
                    "background_q20": density_quantile(edges, bkg, 0.20),
                    "background_q50": density_quantile(edges, bkg, 0.50),
                    "background_q80": density_quantile(edges, bkg, 0.80),
                    "density_overlap_integral": overlap,
                    "standardized_mean_separation": sep,
                    "wp80_threshold_hist": wp80_thr,
                    "wp80_background_fake_hist": wp80_fake,
                    "wp80_fake_delta_vs_90_10": 0.0 if base_wp_fake is None else wp80_fake - base_wp_fake,
                    "background_survival_score_ge_0p5": survival_at_threshold(edges, bkg, 0.5),
                    "signal_survival_score_ge_0p5": survival_at_threshold(edges, sig, 0.5),
                }
            )

    write_csv(args.outdir / "splitstudy_deep_dive_metrics.csv", rows)

    focus = [r for r in rows if r["centrality"] == "0-20%"]
    summary = {
        "input": str(args.input),
        "caveat": "This first-pass diagnostic uses compact histogram summaries, not paired per-event score caches.",
        "focus_0_20": focus,
        "main_findings": [
            "The 0-20% per-centrality AUC rise is small: +0.00068 for 50/50 and +0.00173 for 10/90 relative to 90/10.",
            "Approximate independent-sample AUC significance is at most order 1 sigma from compact counts; paired score-cache bootstrap is needed for a decisive claim.",
            "The histogram-estimated WP80 background fake rate does not improve monotonically with the AUC rise in 0-20%.",
        ],
    }
    (args.outdir / "splitstudy_deep_dive_summary.json").write_text(json.dumps(summary, indent=2))

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )

    fig, axs = plt.subplots(2, 2, figsize=(11.5, 8.0), dpi=180)
    fig.suptitle("Split-study first-pass diagnostics from compact score histograms", fontsize=17, fontweight="bold", x=0.04, ha="left")
    ax_auc, ax_delta, ax_wp, ax_sep = axs.ravel()
    xs = [SPLIT_TRAIN_FRACTION[s] for s in splits]
    for cent in cents:
        yr = [next(r for r in rows if r["centrality"] == cent and r["split"] == s) for s in splits]
        ax_auc.plot(xs, [r["auc"] for r in yr], marker="o", label=cent, lw=2)
        ax_delta.plot(xs, [r["auc_delta_vs_90_10"] for r in yr], marker="o", label=cent, lw=2)
        ax_wp.plot(xs, [r["wp80_background_fake_hist"] for r in yr], marker="o", label=cent, lw=2)
        ax_sep.plot(xs, [r["standardized_mean_separation"] for r in yr], marker="o", label=cent, lw=2)
    for ax in axs.ravel():
        ax.grid(True, color=GRID)
        ax.set_xlabel("training fraction")
        ax.invert_xaxis()
    ax_auc.set_ylabel("AUC")
    ax_delta.axhline(0, color=MUTED, lw=1)
    ax_delta.set_ylabel("AUC delta vs 90/10")
    ax_wp.set_ylabel("hist-estimated WP80 bkg fake")
    ax_sep.set_ylabel("standardized mean score gap")
    ax_auc.legend(frameon=False, fontsize=9)
    fig.text(
        0.04,
        0.025,
        "Caveat: compact histograms support trend triage; paired event-level score caches are needed for a final overfitting claim.",
        fontsize=9.5,
        color=MUTED,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.94])
    fig.savefig(args.outdir / "splitstudy_deep_dive_auc_wp80_summary.png")
    plt.close(fig)

    fig, axs = plt.subplots(2, 2, figsize=(11.5, 7.6), dpi=180)
    fig.suptitle("0-20% centrality: what actually changes with train fraction?", fontsize=17, fontweight="bold", x=0.04, ha="left")
    ax_s, ax_b, ax_ds, ax_db = axs.ravel()
    base = hist_map[(baseline, "0-20%")]
    base_edges = np.asarray(base["bin_edges"], dtype=float)
    centers = 0.5 * (base_edges[:-1] + base_edges[1:])
    base_sig = np.asarray(base["signal_density"], dtype=float)
    base_bkg = np.asarray(base["background_density"], dtype=float)
    for split in splits:
        h = hist_map[(split, "0-20%")]
        edges = np.asarray(h["bin_edges"], dtype=float)
        ctr = 0.5 * (edges[:-1] + edges[1:])
        sig = np.asarray(h["signal_density"], dtype=float)
        bkg = np.asarray(h["background_density"], dtype=float)
        ax_s.step(ctr, sig, where="mid", label=split, color=COLORS[split], lw=2)
        ax_b.step(ctr, bkg, where="mid", label=split, color=COLORS[split], lw=2)
        if split != baseline:
            ax_ds.step(ctr, sig - base_sig, where="mid", label=f"{split} - 90/10", color=COLORS[split], lw=2)
            ax_db.step(ctr, bkg - base_bkg, where="mid", label=f"{split} - 90/10", color=COLORS[split], lw=2)
    ax_s.set_title("Signal score density")
    ax_b.set_title("Background score density")
    ax_ds.set_title("Signal density difference")
    ax_db.set_title("Background density difference")
    for ax in axs.ravel():
        ax.grid(True, color=GRID)
        ax.set_xlabel("BDT score")
    ax_s.set_ylabel("density")
    ax_b.set_ylabel("density")
    ax_ds.set_ylabel("delta density")
    ax_db.set_ylabel("delta density")
    ax_s.legend(frameon=False, fontsize=9)
    ax_ds.axhline(0, color=MUTED, lw=1)
    ax_db.axhline(0, color=MUTED, lw=1)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(args.outdir / "splitstudy_0to20_score_shape_deltas.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.5, 3.8), dpi=180)
    focus_rows = [r for r in rows if r["centrality"] == "0-20%"]
    labels = [r["split"] for r in focus_rows]
    auc_delta = [r["auc_delta_vs_90_10"] for r in focus_rows]
    wp_delta = [r["wp80_fake_delta_vs_90_10"] for r in focus_rows]
    x = np.arange(len(labels))
    ax.bar(x - 0.18, auc_delta, width=0.36, color=PURPLE, label="AUC delta")
    ax.bar(x + 0.18, wp_delta, width=0.36, color=RED, label="WP80 bkg fake delta")
    ax.axhline(0, color=INK, lw=1)
    ax.set_xticks(x, labels)
    ax.set_ylabel("delta vs 90/10")
    ax.set_title("0-20%: AUC rise is tiny and not mirrored by a clean WP80 improvement", fontweight="bold")
    ax.grid(True, axis="y", color=GRID)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(args.outdir / "splitstudy_0to20_auc_vs_wp80_delta.png")
    plt.close(fig)

    print(args.outdir)


if __name__ == "__main__":
    main()

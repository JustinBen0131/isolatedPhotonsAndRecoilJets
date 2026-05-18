#!/usr/bin/env python3
"""Make a 2x3 BDT-vs-MLP score-separation comparison for the baseline AuAu models."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


CENT_BINS = [("0_20", "0-20%"), ("20_50", "20-50%"), ("50_80", "50-80%")]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bdt-validation", type=Path, required=True)
    ap.add_argument("--mlp-validation", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    return ap.parse_args()


def auc_from_scores(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=np.int8)
    score = np.asarray(score, dtype=np.float64)
    mask = np.isfinite(score) & np.isin(y, [0, 1])
    y = y[mask]
    score = score[mask]
    n_sig = int(np.sum(y == 1))
    n_bkg = int(np.sum(y == 0))
    if n_sig == 0 or n_bkg == 0:
        return float("nan")
    order = np.argsort(score, kind="mergesort")
    sorted_score = score[order]
    ranks = np.empty(score.size, dtype=np.float64)
    i = 0
    while i < score.size:
        j = i + 1
        while j < score.size and sorted_score[j] == sorted_score[i]:
            j += 1
        avg_rank = 0.5 * (i + 1 + j)
        ranks[order[i:j]] = avg_rank
        i = j
    rank_sum_sig = float(np.sum(ranks[y == 1]))
    return (rank_sum_sig - n_sig * (n_sig + 1) / 2.0) / (n_sig * n_bkg)


def load_mlp_hist_from_csv(path: Path) -> dict[str, dict[str, np.ndarray | float | int]]:
    rows: dict[str, list[dict[str, str]]] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle):
            rows.setdefault(row["cent_bin"], []).append(row)
    out = {}
    for cent, cent_rows in rows.items():
        cent_rows = sorted(cent_rows, key=lambda r: float(r["bin_lo"]))
        out[cent] = {
            "edges": np.array([float(r["bin_lo"]) for r in cent_rows] + [float(cent_rows[-1]["bin_hi"])]),
            "signal_density": np.array([float(r["signal_density"]) for r in cent_rows]),
            "background_density": np.array([float(r["background_density"]) for r in cent_rows]),
            "auc": float(cent_rows[0]["auc"]),
            "signal_entries": int(cent_rows[0]["signal_entries"]),
            "background_entries": int(cent_rows[0]["background_entries"]),
        }
    return out


def load_mlp_hist_from_score_caches(validation_dir: Path, edges: np.ndarray, score_key: str) -> dict[str, dict[str, np.ndarray | float | int]]:
    manifest = validation_dir / "score_caches.list"
    if not manifest.exists():
        raise FileNotFoundError(f"Missing {manifest}")
    paths = [Path(line.strip()) for line in manifest.read_text().splitlines() if line.strip()]
    if not paths:
        raise FileNotFoundError(f"No score caches listed in {manifest}")

    payload = {}
    for cent_key, cent_label in CENT_BINS:
        lo, hi = [float(x) for x in cent_label.rstrip("%").split("-")]
        sig_scores: list[np.ndarray] = []
        bkg_scores: list[np.ndarray] = []
        for cache in paths:
            with np.load(cache, allow_pickle=True) as z:
                missing = [k for k in ("is_signal", "centrality", "cluster_Et", score_key) if k not in z.files]
                if missing:
                    raise KeyError(f"{cache} missing keys {missing}")
                y = np.asarray(z["is_signal"], dtype=np.int8)
                cent = np.asarray(z["centrality"], dtype=np.float32)
                et = np.asarray(z["cluster_Et"], dtype=np.float32)
                score = np.asarray(z[score_key], dtype=np.float32)
                mask = (
                    np.isfinite(cent)
                    & np.isfinite(et)
                    & np.isfinite(score)
                    & (cent >= lo)
                    & (cent < hi)
                    & (et >= 15.0)
                    & (et < 35.0)
                    & np.isin(y, [0, 1])
                )
                if not np.any(mask):
                    continue
                sig_scores.append(score[mask & (y == 1)])
                bkg_scores.append(score[mask & (y == 0)])
        sig = np.concatenate(sig_scores) if sig_scores else np.array([], dtype=np.float32)
        bkg = np.concatenate(bkg_scores) if bkg_scores else np.array([], dtype=np.float32)
        sig_density = np.histogram(sig, bins=edges, density=True)[0] if sig.size else np.zeros(len(edges) - 1)
        bkg_density = np.histogram(bkg, bins=edges, density=True)[0] if bkg.size else np.zeros(len(edges) - 1)
        y_all = np.r_[np.ones(sig.size, dtype=np.int8), np.zeros(bkg.size, dtype=np.int8)]
        score_all = np.r_[sig, bkg]
        payload[cent_key] = {
            "edges": edges,
            "signal_density": sig_density,
            "background_density": bkg_density,
            "auc": auc_from_scores(y_all, score_all),
            "signal_entries": int(sig.size),
            "background_entries": int(bkg.size),
        }
    return payload


def load_mlp_hist(validation_dir: Path, edges: np.ndarray, score_key: str) -> dict[str, dict[str, np.ndarray | float | int]]:
    csv_path = validation_dir / f"coarse_centrality_score_histograms_{score_key.removeprefix('score_')}.csv"
    if csv_path.exists():
        return load_mlp_hist_from_csv(csv_path)
    return load_mlp_hist_from_score_caches(validation_dir, edges, score_key)


def step_plot(ax, edges, density, *, color, label):
    y = np.r_[density, density[-1]]
    ax.step(edges, y, where="post", color=color, lw=2.0, label=label)
    ax.fill_between(edges, y, step="post", color=color, alpha=0.08)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    bdt_hist = json.loads((args.bdt_validation / "score_histograms.json").read_text())
    edges = np.array(bdt_hist["bin_edges"], dtype=float)
    bdt_metrics = json.loads((args.bdt_validation / "validation_metrics.json").read_text())
    mlp_hist = load_mlp_hist(args.mlp_validation, edges, "score_globalEtCent1535_mlp_noIso")

    bdt_product = bdt_hist["products"]["globalEtCent1535_bdt_noIso"]["by_centrality"]
    bdt_auc = bdt_metrics["products"]["globalEtCent1535_bdt_noIso"]["auc_by_centrality"]

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.linewidth": 1.05,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "axes.grid": True,
            "grid.color": "0.90",
            "grid.linewidth": 0.55,
        }
    )

    fig, axes = plt.subplots(2, 3, figsize=(16.0, 8.6), sharex=True)
    summary = []
    row_defs = [
        ("Global BDT", "Classifier score", bdt_product),
        ("Global MLP", "Classifier score", mlp_hist),
    ]
    sig_color = "#1f77b4"
    bkg_color = "#d62728"

    for irow, (model_label, x_label, payload) in enumerate(row_defs):
        for icol, (cent_key, cent_label) in enumerate(CENT_BINS):
            ax = axes[irow, icol]
            if irow == 0:
                cent_payload = payload[cent_key]
                sig_density = np.array(cent_payload["signal"]["density"], dtype=float)
                bkg_density = np.array(cent_payload["background"]["density"], dtype=float)
                auc = float(bdt_auc[cent_key])
                sig_entries = int(cent_payload["signal"]["entries"])
                bkg_entries = int(cent_payload["background"]["entries"])
                local_edges = edges
            else:
                cent_payload = payload[cent_key]
                sig_density = np.array(cent_payload["signal_density"], dtype=float)
                bkg_density = np.array(cent_payload["background_density"], dtype=float)
                auc = float(cent_payload["auc"])
                sig_entries = int(cent_payload["signal_entries"])
                bkg_entries = int(cent_payload["background_entries"])
                local_edges = np.array(cent_payload["edges"], dtype=float)

            step_plot(ax, local_edges, sig_density, color=sig_color, label="Signal")
            step_plot(ax, local_edges, bkg_density, color=bkg_color, label="Background")
            ymax = max(float(np.nanmax(sig_density)), float(np.nanmax(bkg_density)))
            ax.set_ylim(0.0, ymax * 1.22 if ymax > 0 else 1.0)
            ax.set_xlim(0.0, 1.0)
            ax.text(
                0.045,
                0.88,
                f"AUC = {auc:.3f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=13,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.80", alpha=0.92),
            )
            if irow == 0:
                ax.set_title(f"Centrality {cent_label}", fontsize=15, fontweight="bold", pad=10)
            if icol == 0:
                ax.set_ylabel("Area-normalized density", fontsize=13)
                ax.text(
                    -0.18,
                    0.50,
                    model_label,
                    transform=ax.transAxes,
                    rotation=90,
                    ha="center",
                    va="center",
                    fontsize=15,
                    fontweight="bold",
                )
            if irow == 1:
                ax.set_xlabel(x_label, fontsize=13)
            if irow == 0 and icol == 2:
                ax.legend(loc="upper right", bbox_to_anchor=(1.0, 1.30), ncol=2, frameon=False, fontsize=12)
            summary.append(
                {
                    "model": model_label,
                    "cent_bin": cent_key,
                    "auc": auc,
                    "signal_entries": sig_entries,
                    "background_entries": bkg_entries,
                }
            )

    fig.text(0.040, 0.975, r"$\bf{\it{sPHENIX}}$ Internal", ha="left", va="top", fontsize=15.5)
    fig.text(0.040, 0.945, r"Au+Au embedded photon-ID validation", ha="left", va="top", fontsize=11.5)
    fig.text(0.040, 0.920, r"$15 < E_T < 35$ GeV", ha="left", va="top", fontsize=11.5)
    fig.text(
        0.57,
        0.975,
        "Signal/background score separation: BDT vs MLP",
        ha="center",
        va="top",
        fontsize=18.5,
        fontweight="bold",
    )
    fig.text(
        0.57,
        0.945,
        "Curves are area-normalized; panel labels show AUC.",
        ha="center",
        va="top",
        fontsize=11.5,
        color="0.30",
    )
    fig.text(
        0.57,
        0.020,
        "Both rows use the baseline feature set.",
        ha="center",
        va="bottom",
        fontsize=10.5,
        color="0.35",
    )
    fig.tight_layout(rect=[0.110, 0.065, 0.995, 0.875], h_pad=1.25, w_pad=1.75)

    png = args.outdir / "global_noiso_bdt_vs_mlp_score_separation_by_centrality.png"
    csv_path = args.outdir / "global_noiso_bdt_vs_mlp_score_separation_by_centrality.csv"
    fig.savefig(png, dpi=220)
    plt.close(fig)
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)
    print(png)
    print(csv_path)


if __name__ == "__main__":
    main()

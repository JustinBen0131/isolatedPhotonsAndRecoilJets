#!/usr/bin/env python3
"""Compare signal/background BDT-score separation with and without isolation inputs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CENT_BINS = [(0.0, 20.0, "0-20%"), (20.0, 50.0, "20-50%"), (50.0, 80.0, "50-80%")]
MODEL_ROWS = [
    {
        "label": "Baseline inputs",
        "score_key": "score_globalEtCent1535_bdt_noIso_ptCent7",
        "slug": "baseline_inputs",
    },
    {
        "label": "Baseline + isolation inputs",
        "score_key": "score_globalEtCent1535_bdt_iso_ptCent7",
        "slug": "baseline_plus_isolation_inputs",
    },
]


def score_cache_paths(path: Path) -> list[Path]:
    if path.is_dir():
        manifest = path / "score_caches.local.list"
        if not manifest.is_file():
            manifest = path / "score_caches.list"
    else:
        manifest = path
    if manifest.suffix == ".npz" and manifest.is_file():
        return [manifest]
    if not manifest.is_file():
        raise SystemExit(f"Missing score-cache manifest: {manifest}")
    paths = [Path(line.strip()) for line in manifest.read_text().splitlines() if line.strip()]
    if not paths:
        raise SystemExit(f"No score-cache paths in {manifest}")
    return paths


def auc_from_scores(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=np.int8)
    score = np.asarray(score, dtype=np.float64)
    mask = np.isin(y, [0, 1]) & np.isfinite(score)
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


def load_model_split(
    paths: list[Path], score_key: str, bins: np.ndarray
) -> tuple[dict[str, dict[str, object]], list[dict[str, object]]]:
    payload: dict[str, dict[str, object]] = {}
    rows: list[dict[str, object]] = []
    for clo, chi, cent_label in CENT_BINS:
        sig_parts: list[np.ndarray] = []
        bkg_parts: list[np.ndarray] = []
        for path in paths:
            with np.load(path, allow_pickle=True) as z:
                missing = [key for key in ("is_signal", "cluster_Et", "centrality", score_key) if key not in z.files]
                if missing:
                    raise SystemExit(f"{path} missing {missing}; available={z.files}")
                y = np.asarray(z["is_signal"], dtype=np.int8)
                et = np.asarray(z["cluster_Et"], dtype=np.float32)
                cent = np.asarray(z["centrality"], dtype=np.float32)
                score = np.asarray(z[score_key], dtype=np.float32)
                mask = (
                    np.isin(y, [0, 1])
                    & np.isfinite(et)
                    & np.isfinite(cent)
                    & np.isfinite(score)
                    & (et >= 15.0)
                    & (et < 35.0)
                    & (cent >= clo)
                    & (cent < chi)
                )
                if np.any(mask & (y == 1)):
                    sig_parts.append(score[mask & (y == 1)])
                if np.any(mask & (y == 0)):
                    bkg_parts.append(score[mask & (y == 0)])
        sig = np.concatenate(sig_parts) if sig_parts else np.array([], dtype=np.float32)
        bkg = np.concatenate(bkg_parts) if bkg_parts else np.array([], dtype=np.float32)
        sig_density = np.histogram(sig, bins=bins, density=True)[0] if sig.size else np.zeros(len(bins) - 1)
        bkg_density = np.histogram(bkg, bins=bins, density=True)[0] if bkg.size else np.zeros(len(bins) - 1)
        all_y = np.r_[np.ones(sig.size, dtype=np.int8), np.zeros(bkg.size, dtype=np.int8)]
        all_score = np.r_[sig, bkg]
        auc = auc_from_scores(all_y, all_score)
        record = {
            "centrality_bin": cent_label,
            "score_key": score_key,
            "signal_entries": int(sig.size),
            "background_entries": int(bkg.size),
            "auc": float(auc),
            "signal_mean_score": float(np.mean(sig)) if sig.size else float("nan"),
            "background_mean_score": float(np.mean(bkg)) if bkg.size else float("nan"),
            "signal_median_score": float(np.median(sig)) if sig.size else float("nan"),
            "background_median_score": float(np.median(bkg)) if bkg.size else float("nan"),
        }
        payload[cent_label] = {
            "signal_density": sig_density,
            "background_density": bkg_density,
            "edges": bins,
            **record,
        }
        rows.append(record)
    return payload, rows


def add_sphenix_label(fig: plt.Figure) -> None:
    fig.text(0.047, 0.872, "sPHENIX", ha="left", va="center", fontsize=18, fontstyle="italic", fontweight="bold")
    fig.text(0.124, 0.872, " Internal", ha="left", va="center", fontsize=18)


def step_density(ax: plt.Axes, edges: np.ndarray, density: np.ndarray, *, color: str, label: str, linestyle: str) -> None:
    y = np.r_[density, density[-1]]
    ax.step(edges, y, where="post", color=color, linewidth=2.3, linestyle=linestyle, label=label)
    ax.fill_between(edges, y, step="post", color=color, alpha=0.075)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--tag", default="binned_bdt_ptcent7_signal_background_score_split_iso_input_comparison")
    args = ap.parse_args()

    paths = score_cache_paths(args.manifest)
    bins = np.linspace(0.0, 1.0, 51)
    all_payload: dict[str, dict[str, dict[str, object]]] = {}
    all_rows: list[dict[str, object]] = []
    for row in MODEL_ROWS:
        payload, rows = load_model_split(paths, row["score_key"], bins)
        all_payload[row["slug"]] = payload
        for item in rows:
            item = dict(item)
            item["model_row"] = row["label"]
            item["model_slug"] = row["slug"]
            all_rows.append(item)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.linewidth": 1.05,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig, axes = plt.subplots(2, 3, figsize=(17.5, 9.0), sharex=True, sharey=False, dpi=185)
    fig.patch.set_facecolor("white")
    sig_color = "#0072B2"
    bkg_color = "#D55E00"

    for irow, model in enumerate(MODEL_ROWS):
        for icol, (_, _, cent_label) in enumerate(CENT_BINS):
            ax = axes[irow, icol]
            item = all_payload[model["slug"]][cent_label]
            step_density(ax, item["edges"], item["signal_density"], color=sig_color, label="Signal: Photon12+20", linestyle="-")
            step_density(
                ax,
                item["edges"],
                item["background_density"],
                color=bkg_color,
                label="Background: Jet12+20+30",
                linestyle="--",
            )
            ymax = max(float(np.nanmax(item["signal_density"])), float(np.nanmax(item["background_density"])))
            ax.set_ylim(0.0, max(2.0, 1.24 * ymax))
            ax.set_xlim(0.0, 1.0)
            ax.grid(True, color="#D1D5DB", linewidth=0.95, alpha=0.70)
            ax.minorticks_on()
            ax.tick_params(labelsize=10.2)
            if irow == 0:
                ax.set_title(f"{cent_label} centrality", fontsize=15.2, fontweight="bold", pad=7)
            if irow == 1:
                ax.set_xlabel("BDT score", fontsize=13.0)
            if icol == 0:
                ax.set_ylabel("Unit-normalized candidates", fontsize=12.2)
            ax.text(
                0.045,
                0.91,
                f"AUC = {item['auc']:.3f}\nS = {item['signal_entries']:,}\nB = {item['background_entries']:,}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=9.6,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.86, "pad": 2.7},
            )

    for irow, model in enumerate(MODEL_ROWS):
        pos = axes[irow, 0].get_position()
        y = 0.5 * (pos.y0 + pos.y1)
        fig.text(
            0.019,
            y,
            model["label"],
            ha="center",
            va="center",
            rotation=90,
            fontsize=14.0,
            fontweight="bold",
        )

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.52, 0.021),
        ncol=2,
        frameon=False,
        fontsize=12,
        handlelength=2.7,
        columnspacing=2.5,
    )

    fig.text(0.045, 0.965, "Signal/background BDT-score separation", ha="left", va="top", fontsize=24, fontweight="bold")
    fig.text(
        0.045,
        0.928,
        "8 $E_T$ x 7 centrality binned BDTs; top: baseline inputs, bottom: baseline + isolation inputs; 15 < cluster $E_T$ < 35 GeV",
        ha="left",
        va="top",
        fontsize=13.2,
        color="#374151",
    )
    add_sphenix_label(fig)
    fig.tight_layout(rect=[0.052, 0.105, 0.995, 0.832], h_pad=1.9, w_pad=1.35)

    args.outdir.mkdir(parents=True, exist_ok=True)
    png = args.outdir / f"{args.tag}.png"
    csv_path = args.outdir / f"{args.tag}.csv"
    json_path = args.outdir / f"{args.tag}.json"
    fig.savefig(png)
    plt.close(fig)

    fieldnames = [
        "model_row",
        "model_slug",
        "centrality_bin",
        "score_key",
        "signal_entries",
        "background_entries",
        "auc",
        "signal_mean_score",
        "background_mean_score",
        "signal_median_score",
        "background_median_score",
    ]
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)
    json_path.write_text(
        json.dumps(
            {
                "schema": "AUAU_BDT_ISO_INPUT_COMPARISON_SIGNAL_BACKGROUND_SPLIT_V1",
                "manifest": str(args.manifest),
                "models": MODEL_ROWS,
                "signal": "embedded Photon12+20",
                "background": "embedded inclusive Jet12+20+30",
                "pt_range_GeV": [15.0, 35.0],
                "rows": all_rows,
                "outputs": {"png": str(png), "csv": str(csv_path), "json": str(json_path)},
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

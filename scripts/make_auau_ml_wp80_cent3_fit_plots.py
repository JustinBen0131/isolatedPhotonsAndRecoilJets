#!/usr/bin/env python3
"""Plot target-80 score thresholds vs E_T for AuAu ML score-cache products."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


PT_EDGES = [15, 17, 19, 21, 23, 25, 27, 30, 35]
CENT3_EDGES = [0, 20, 50, 80]


def parse_edges(text: str, default: list[float]) -> list[float]:
    if not text:
        return [float(x) for x in default]
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def score_cache_paths(path: Path) -> list[Path]:
    if path.is_dir():
        local = path / "score_caches.local.list"
        manifest = local if local.is_file() else path / "score_caches.list"
    else:
        manifest = path
    if manifest.is_file() and manifest.suffix == ".list":
        paths = [Path(line.strip()) for line in manifest.read_text().splitlines() if line.strip()]
    elif manifest.is_file() and manifest.suffix == ".npz":
        paths = [manifest]
    else:
        raise SystemExit(f"Cannot find score caches from {path}")
    if not paths:
        raise SystemExit(f"No score cache entries in {manifest}")
    return paths


def load_score_arrays(manifest: Path, score_key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    paths = score_cache_paths(manifest)
    y_parts: list[np.ndarray] = []
    et_parts: list[np.ndarray] = []
    cent_parts: list[np.ndarray] = []
    score_parts: list[np.ndarray] = []
    for path in paths:
        with np.load(path, allow_pickle=True) as z:
            missing = [k for k in ("is_signal", "cluster_Et", "centrality", score_key) if k not in z.files]
            if missing:
                raise SystemExit(f"{path} missing keys: {missing}; available={z.files}")
            y_parts.append(np.asarray(z["is_signal"], dtype=np.int8))
            et_parts.append(np.asarray(z["cluster_Et"], dtype=np.float32))
            cent_parts.append(np.asarray(z["centrality"], dtype=np.float32))
            score_parts.append(np.asarray(z[score_key], dtype=np.float32))
    return (
        np.concatenate(y_parts),
        np.concatenate(et_parts),
        np.concatenate(cent_parts),
        np.concatenate(score_parts),
    )


def threshold_for_signal_efficiency(y: np.ndarray, score: np.ndarray, target: float) -> dict | None:
    valid = np.isin(y, [0, 1]) & np.isfinite(score)
    sig = score[valid & (y == 1)]
    bkg = score[valid & (y == 0)]
    if sig.size <= 0 or bkg.size <= 0:
        return None
    threshold = float(np.quantile(sig, max(0.0, min(1.0, 1.0 - target))))
    return {
        "threshold": threshold,
        "signal_efficiency": float(np.mean(sig > threshold)),
        "background_fake_rate": float(np.mean(bkg > threshold)),
        "signal_entries": int(sig.size),
        "background_entries": int(bkg.size),
    }


def derive_cells(
    y: np.ndarray,
    et: np.ndarray,
    cent: np.ndarray,
    score: np.ndarray,
    pt_edges: list[float],
    cent_edges: list[float],
    target: float,
) -> tuple[list[dict], dict]:
    finite = np.isfinite(et) & np.isfinite(cent) & np.isfinite(score) & np.isin(y, [0, 1])
    finite &= (et >= pt_edges[0]) & (et < pt_edges[-1])
    finite &= (cent >= cent_edges[0]) & (cent < cent_edges[-1])
    inclusive = threshold_for_signal_efficiency(y[finite], score[finite], target)
    cells: list[dict] = []
    for clo, chi in zip(cent_edges[:-1], cent_edges[1:]):
        cmask = finite & (cent >= clo) & (cent < chi)
        cent_fallback = threshold_for_signal_efficiency(y[cmask], score[cmask], target)
        for plo, phi in zip(pt_edges[:-1], pt_edges[1:]):
            mask = cmask & (et >= plo) & (et < phi)
            item = threshold_for_signal_efficiency(y[mask], score[mask], target)
            source = "cell"
            if item is None:
                item = cent_fallback
                source = "centrality_fallback"
            if item is None:
                item = inclusive
                source = "inclusive_fallback"
            if item is None:
                continue
            cells.append(
                {
                    "centrality_min": float(clo),
                    "centrality_max": float(chi),
                    "pt_min": float(plo),
                    "pt_max": float(phi),
                    "pt_center": float(0.5 * (plo + phi)),
                    "source": source,
                    **item,
                }
            )
    return cells, inclusive or {}


def fit_by_centrality(cells: list[dict]) -> dict[str, dict]:
    fits: dict[str, dict] = {}
    for clo, chi in zip(CENT3_EDGES[:-1], CENT3_EDGES[1:]):
        key = f"{clo:g}_{chi:g}"
        rows = [
            r
            for r in cells
            if r["centrality_min"] == float(clo)
            and r["centrality_max"] == float(chi)
            and r["source"] == "cell"
            and math.isfinite(r["threshold"])
        ]
        if len(rows) < 2:
            fits[key] = {
                "slope": math.nan,
                "intercept": math.nan,
                "max_abs_residual": math.nan,
                "n_fit_points": len(rows),
            }
            continue
        x = np.array([r["pt_center"] for r in rows], dtype=float)
        t = np.array([r["threshold"] for r in rows], dtype=float)
        w = np.sqrt(np.array([max(1, r["signal_entries"]) for r in rows], dtype=float))
        coeff = np.polyfit(x, t, 1, w=w)
        pred = coeff[0] * x + coeff[1]
        fits[key] = {
            "slope": float(coeff[0]),
            "intercept": float(coeff[1]),
            "max_abs_residual": float(np.max(np.abs(t - pred))),
            "n_fit_points": int(len(rows)),
        }
    return fits


def save_cells_csv(path: Path, cells: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "centrality_min",
        "centrality_max",
        "pt_min",
        "pt_max",
        "pt_center",
        "threshold",
        "signal_efficiency",
        "background_fake_rate",
        "signal_entries",
        "background_entries",
        "source",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(cells)


def add_sphenix_label(fig: plt.Figure) -> None:
    fig.text(0.057, 0.855, "sPHENIX", ha="left", va="center", fontsize=18, fontstyle="italic", fontweight="bold")
    fig.text(0.135, 0.855, " Internal", ha="left", va="center", fontsize=18)


def plot_cells(cells: list[dict], fits: dict[str, dict], out_png: Path, model_label: str, target: float) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(17.4, 5.9), sharex=True, sharey=True, dpi=185)
    fig.patch.set_facecolor("white")
    colors = ["#1F77B4", "#009E73", "#D55E00"]
    all_thresholds = [r["threshold"] for r in cells if math.isfinite(r["threshold"])]
    ylo = max(0.0, min(all_thresholds) - 0.055)
    yhi = min(1.0, max(all_thresholds) + 0.075)
    xfit = np.linspace(PT_EDGES[0], PT_EDGES[-1], 160)

    for idx, (ax, clo, chi) in enumerate(zip(axes, CENT3_EDGES[:-1], CENT3_EDGES[1:])):
        rows = [r for r in cells if r["centrality_min"] == float(clo) and r["centrality_max"] == float(chi)]
        x = np.array([r["pt_center"] for r in rows], dtype=float)
        y = np.array([r["threshold"] for r in rows], dtype=float)
        eff = np.array([r["signal_efficiency"] for r in rows], dtype=float)
        sig_n = np.array([max(1, r["signal_entries"]) for r in rows], dtype=float)
        score_err = np.sqrt(np.maximum(eff * (1.0 - eff), 0.0) / sig_n)
        ax.errorbar(
            x,
            y,
            yerr=score_err,
            fmt="o",
            color=colors[idx],
            markeredgecolor="black",
            markeredgewidth=0.5,
            markersize=6.8,
            capsize=2.7,
            linestyle="none",
            label="cell target-80 threshold",
        )
        fit = fits.get(f"{clo:g}_{chi:g}", {})
        if math.isfinite(fit.get("slope", math.nan)) and math.isfinite(fit.get("intercept", math.nan)):
            ax.plot(
                xfit,
                fit["slope"] * xfit + fit["intercept"],
                color="#111111",
                linewidth=2.0,
                linestyle="--",
                label="linear fit",
            )
            fit_text = f"threshold = {fit['intercept']:.3f} {fit['slope']:+.4f} E_T\nmax residual = {fit['max_abs_residual']:.3f}"
        else:
            fit_text = "fit unavailable"
        ax.text(
            0.045,
            0.935,
            fit_text,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=10.6,
            color="#111111",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.86, "pad": 3},
        )
        ax.set_title(f"{clo:g}-{chi:g}% centrality", fontsize=15, fontweight="bold")
        ax.set_xlim(14.5, 35.5)
        ax.set_ylim(ylo, yhi)
        ax.grid(True, color="#D1D5DB", linewidth=1.0, alpha=0.8)
        ax.tick_params(direction="in", which="both", top=True, right=True, labelsize=12)
        ax.minorticks_on()
        ax.set_xlabel(r"cluster $E_T$ bin center [GeV]", fontsize=14)
        if idx == 0:
            ax.set_ylabel(f"score threshold for {100 * target:.0f}% signal efficiency", fontsize=14)
        if idx == 2:
            ax.legend(loc="lower right", frameon=False, fontsize=11)

    fig.text(0.055, 0.94, f"{model_label}: target-80 score cuts", ha="left", va="top", fontsize=24, fontweight="bold")
    fig.text(
        0.055,
        0.897,
        "Derived from validation score caches: Photon12+20 signal vs Jet12+20+30 inclusive background",
        ha="left",
        va="top",
        fontsize=14,
        color="#4B5563",
    )
    add_sphenix_label(fig)
    fig.tight_layout(rect=[0.045, 0.085, 0.995, 0.80])
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--score-key", required=True)
    ap.add_argument("--model-label", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--target", type=float, default=0.80)
    ap.add_argument("--pt-edges", default="")
    args = ap.parse_args()

    pt_edges = parse_edges(args.pt_edges, PT_EDGES)
    y, et, cent, score = load_score_arrays(args.manifest, args.score_key)
    cells, inclusive = derive_cells(y, et, cent, score, pt_edges, [float(x) for x in CENT3_EDGES], args.target)
    fits = fit_by_centrality(cells)

    png = args.outdir / f"{args.tag}_wp80_cent3_threshold_fit_vs_et.png"
    csv_path = args.outdir / f"{args.tag}_wp80_cent3_threshold_cells.csv"
    json_path = args.outdir / f"{args.tag}_wp80_cent3_threshold_summary.json"
    save_cells_csv(csv_path, cells)
    plot_cells(cells, fits, png, args.model_label, args.target)
    json_path.write_text(
        json.dumps(
            {
                "schema": "AUAU_ML_WP80_CENT3_FIT_V1",
                "manifest": str(args.manifest),
                "score_key": args.score_key,
                "model_label": args.model_label,
                "tag": args.tag,
                "target_signal_efficiency": args.target,
                "pt_edges": pt_edges,
                "cent_edges": CENT3_EDGES,
                "inclusive": inclusive,
                "fits": fits,
                "cells": cells,
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

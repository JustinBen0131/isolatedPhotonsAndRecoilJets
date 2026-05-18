#!/usr/bin/env python3
"""Plot E_T-dependent WP80 BDT thresholds in centrality bins.

This is intentionally validation-side only: it reads score-cache npz files,
derives score thresholds that retain the requested signal efficiency in each
E_T x centrality cell, and fits threshold vs E_T independently per centrality
bin. It does not inspect or modify RecoilJets production ROOT outputs.
"""

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
CENT7_EDGES = [0, 10, 20, 30, 40, 50, 60, 80]


def parse_edges(text: str, default: list[float]) -> list[float]:
    if not text:
        return [float(x) for x in default]
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def threshold_for_signal_efficiency(y: np.ndarray, score: np.ndarray, target: float) -> dict | None:
    mask = np.isin(y, [0, 1]) & np.isfinite(score)
    sig = score[mask & (y == 1)]
    bkg = score[mask & (y == 0)]
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


def load_score_arrays(score_dir: Path, score_key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    paths = sorted(score_dir.glob("score_cache_*.npz"))
    if not paths:
        raise SystemExit(f"No score_cache_*.npz files found under {score_dir}")
    y_parts: list[np.ndarray] = []
    et_parts: list[np.ndarray] = []
    cent_parts: list[np.ndarray] = []
    score_parts: list[np.ndarray] = []
    for path in paths:
        with np.load(path, allow_pickle=True) as z:
            missing = [k for k in ("is_signal", "cluster_Et", "centrality", score_key) if k not in z.files]
            if missing:
                raise SystemExit(f"{path} missing keys: {missing}")
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


def derive_cells(
    y: np.ndarray,
    et: np.ndarray,
    cent: np.ndarray,
    score: np.ndarray,
    pt_edges: list[float],
    cent_edges: list[float],
    target: float,
) -> tuple[list[dict], dict]:
    cells: list[dict] = []
    finite = np.isfinite(et) & np.isfinite(cent) & np.isfinite(score) & np.isin(y, [0, 1])
    inclusive = threshold_for_signal_efficiency(y[finite], score[finite], target)
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
            row = {
                "centrality_min": float(clo),
                "centrality_max": float(chi),
                "pt_min": float(plo),
                "pt_max": float(phi),
                "pt_center": float(0.5 * (plo + phi)),
                "source": source,
                **item,
            }
            cells.append(row)
    return cells, inclusive or {}


def fit_by_centrality(cells: list[dict]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    by_cent: dict[str, list[dict]] = {}
    for row in cells:
        key = f"{row['centrality_min']:g}_{row['centrality_max']:g}"
        by_cent.setdefault(key, []).append(row)
    for key, rows in by_cent.items():
        good = [r for r in rows if r["source"] == "cell" and math.isfinite(r["threshold"])]
        if len(good) >= 2:
            x = np.array([r["pt_center"] for r in good], dtype=float)
            t = np.array([r["threshold"] for r in good], dtype=float)
            w = np.sqrt(np.array([max(1, r["signal_entries"]) for r in good], dtype=float))
            coeff = np.polyfit(x, t, 1, w=w)
            pred = coeff[0] * x + coeff[1]
            out[key] = {
                "slope": float(coeff[0]),
                "intercept": float(coeff[1]),
                "max_abs_residual": float(np.max(np.abs(t - pred))),
                "n_fit_points": len(good),
            }
        else:
            out[key] = {"slope": math.nan, "intercept": math.nan, "max_abs_residual": math.nan, "n_fit_points": len(good)}
    return out


def save_csv(path: Path, rows: list[dict]) -> None:
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
        writer.writerows(rows)


def style_axis(ax, xlabel: str | None, ylabel: str | None) -> None:
    ax.grid(True, color="#d7d7d7", alpha=0.7, linewidth=0.8)
    ax.tick_params(direction="in", which="both", top=True, right=True, labelsize=10)
    ax.minorticks_on()
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=12)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=12)


def plot_thresholds(
    cells: list[dict],
    fits: dict[str, dict],
    cent_edges: list[float],
    out_png: Path,
    target: float,
    model_label: str,
    tag: str,
) -> None:
    n_cent = len(cent_edges) - 1
    if n_cent <= 3:
        nrows, ncols = 1, n_cent
        figsize = (5.3 * ncols, 4.6)
    else:
        nrows, ncols = 2, 4
        figsize = (18.0, 8.2)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)
    axes_arr = np.atleast_1d(axes).ravel()
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#6F3FB5", "#E69F00", "#56B4E9"]
    xfit = np.linspace(PT_EDGES[0], PT_EDGES[-1], 120)
    all_thr = [r["threshold"] for r in cells if math.isfinite(r["threshold"])]
    ylo = max(0.0, min(all_thr) - 0.055) if all_thr else 0.0
    yhi = min(1.0, max(all_thr) + 0.065) if all_thr else 1.0
    for idx, (clo, chi) in enumerate(zip(cent_edges[:-1], cent_edges[1:])):
        ax = axes_arr[idx]
        rows = [r for r in cells if r["centrality_min"] == float(clo) and r["centrality_max"] == float(chi)]
        color = colors[idx % len(colors)]
        x = np.array([r["pt_center"] for r in rows], dtype=float)
        y = np.array([r["threshold"] for r in rows], dtype=float)
        eff = np.array([r["signal_efficiency"] for r in rows], dtype=float)
        sig_n = np.array([r["signal_entries"] for r in rows], dtype=float)
        yerr = np.sqrt(np.maximum(eff * (1.0 - eff), 0.0) / np.maximum(sig_n, 1.0))
        ax.errorbar(x, y, yerr=yerr, fmt="o", color=color, markersize=5.5, capsize=2.2, label="cell threshold")
        f = fits.get(f"{clo:g}_{chi:g}", {})
        if math.isfinite(f.get("slope", math.nan)) and math.isfinite(f.get("intercept", math.nan)):
            ax.plot(xfit, f["slope"] * xfit + f["intercept"], color="#111111", linewidth=1.8, linestyle="--", label="linear fit")
            txt = f"slope={f['slope']:+.4f}/GeV\nmax resid={f['max_abs_residual']:.3f}"
        else:
            txt = "fit unavailable"
        ax.text(0.04, 0.94, txt, transform=ax.transAxes, ha="left", va="top", fontsize=9.0)
        ax.set_title(f"Centrality {clo:g}-{chi:g}%", fontsize=13, fontweight="bold")
        ax.set_ylim(ylo, yhi)
        ax.set_xlim(14.6, 35.4)
        is_bottom_row = nrows == 1 or (idx // ncols) == (nrows - 1)
        style_axis(ax, r"$E_T^\gamma$ bin center [GeV]" if is_bottom_row else None, None)
        if idx == min(2, n_cent - 1):
            ax.legend(loc="best", frameon=False, fontsize=9)
    for ax in axes_arr[n_cent:]:
        ax.axis("off")
    fig.suptitle(
        f"{model_label}: $E_T$-dependent BDT WP{int(round(100*target))} cuts in {n_cent} centrality bins",
        fontsize=15,
        fontweight="bold",
        y=0.985,
    )
    fig.text(0.014, 0.52, "BDT score cut for 80% signal efficiency", va="center", rotation=90, fontsize=13)
    fig.text(
        0.5,
        0.018,
        f"Derived from validation score caches for {tag}; tight selection is score > threshold. These are proposed runtime cuts, not the previous flat-threshold production.",
        ha="center",
        fontsize=9.8,
        color="#444444",
    )
    fig.tight_layout(rect=[0.045, 0.055, 0.995, 0.94])
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=190)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--score-dir", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--score-key", default="score_globalEtCent1535_bdt_noIso")
    ap.add_argument("--model-label", default="Global BDT")
    ap.add_argument("--tag", default="globalEtCent1535_bdt_noIso")
    ap.add_argument("--target", type=float, default=0.80)
    ap.add_argument("--pt-edges", default="")
    args = ap.parse_args()

    pt_edges = parse_edges(args.pt_edges, PT_EDGES)
    y, et, cent, score = load_score_arrays(args.score_dir, args.score_key)
    summary = {
        "score_dir": str(args.score_dir),
        "score_key": args.score_key,
        "target_signal_efficiency": args.target,
        "pt_edges": pt_edges,
        "rows": int(y.size),
    }
    for label, cent_edges in (("cent3", CENT3_EDGES), ("cent7", CENT7_EDGES)):
        cells, inclusive = derive_cells(y, et, cent, score, pt_edges, [float(x) for x in cent_edges], args.target)
        fits = fit_by_centrality(cells)
        save_csv(args.outdir / f"{args.tag}_wp80_{label}_threshold_cells.csv", cells)
        plot_thresholds(
            cells,
            fits,
            [float(x) for x in cent_edges],
            args.outdir / f"{args.tag}_wp80_{label}_threshold_fit_vs_et.png",
            args.target,
            args.model_label,
            args.tag,
        )
        summary[label] = {"cent_edges": cent_edges, "inclusive": inclusive, "fits": fits}
    (args.outdir / f"{args.tag}_wp80_centbin_fit_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"[wp80CentFit] wrote {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Plot interim AuAu b002 shower-shape stage flows from a compact JSON cache."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CACHE = (
    REPO_ROOT
    / "dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_auauDataNewV008_b002_20260612"
    / "interim_complete_runs"
    / "the42_b002_interim_complete_run_shower_shape_hists.json"
)
DEFAULT_OUTDIR = DEFAULT_CACHE.parent

VARS = [
    ("e11_to_e33", r"$E_{11}/E_{33}$", (0.0, 1.0)),
    ("weta_cogx", r"$w_{\eta}^{\mathrm{cogX}}$", (0.0, 2.0)),
    ("e32_to_e35", r"$E_{32}/E_{35}$", (0.4, 1.0)),
]
STAGES = [
    ("cut0", "Before preselection"),
    ("cut1", "After preselection"),
    ("cut2", "After tight ID"),
]
CENTRALITIES = [
    ("cent0_20", "0-20%", "#111111"),
    ("cent20_50", "20-50%", "#1f77b4"),
    ("cent50_80", "50-80%", "#d95f02"),
]


def load_cache(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("schema") != "THE42_B002_INTERIM_COMPLETE_RUN_SHOWER_SHAPES_V1":
        raise RuntimeError(f"Unexpected cache schema in {path}")
    return data


def arrays(payload: dict, var: str, cent: str, stage: str) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    item = payload.get("hists", {}).get(var, {}).get(cent, {}).get(stage)
    if not item:
        return None
    x = np.asarray(item["x"], dtype=float)
    y = np.asarray(item["y"], dtype=float)
    e = np.asarray(item["e"], dtype=float)
    if x.size == 0 or y.size != x.size or e.size != x.size:
        return None
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(e)
    x, y, e = x[mask], y[mask], e[mask]
    total = float(np.sum(y))
    if total <= 0:
        return None
    return x, y / total, e / total


def rebin(x: np.ndarray, y: np.ndarray, e: np.ndarray, factor: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if factor <= 1 or len(x) < factor:
        return x, y, e
    n = len(x) // factor
    trim = n * factor
    xr = x[:trim].reshape(n, factor)
    yr = y[:trim].reshape(n, factor)
    er = e[:trim].reshape(n, factor)
    return xr.mean(axis=1), yr.sum(axis=1), np.sqrt(np.sum(er * er, axis=1))


def select_range(
    arr: tuple[np.ndarray, np.ndarray, np.ndarray],
    xlim: tuple[float, float],
    factor: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x, y, e = rebin(*arr, factor)
    mask = (x >= xlim[0]) & (x <= xlim[1])
    x, y, e = x[mask], y[mask], e[mask]
    total = float(np.sum(y))
    if total > 0:
        y, e = y / total, e / total
    return x, y, e


def add_sphenix_label(ax, *, include_details: bool) -> None:
    ax.text(
        0.035,
        0.945,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11.5,
    )
    if include_details:
        ax.text(
            0.035,
            0.835,
            "Au+Au $\\sqrt{s_{NN}}=200$ GeV\n$15<E_T<35$ GeV, $|\\eta|<0.7$",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=10.5,
            linespacing=1.1,
        )


def draw(cache: dict, outdir: Path) -> list[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(len(VARS), len(STAGES), figsize=(14.6, 10.2), constrained_layout=False)
    fig.patch.set_facecolor("white")
    meta = cache.get("metadata", {})
    complete_runs = meta.get("complete_runs", "?")
    complete_runs_available = meta.get("complete_runs_available")
    total_runs = meta.get("total_runs", "?")
    non_tiny = meta.get("non_tiny_roots", "?")
    root_files = meta.get("root_files", "?")
    fig.suptitle(
        "Interim AuAu b002 shower-shape stage flow from complete runs",
        fontsize=20,
        fontweight="bold",
        y=0.985,
    )
    if complete_runs_available is not None:
        run_text = f"selected complete runs {complete_runs}/{complete_runs_available} currently complete ({total_runs} total b002 runs)"
    else:
        run_text = f"complete runs {complete_runs}/{total_runs}"
    fig.text(
        0.5,
        0.948,
        f"Data only, no active-output merge: {run_text}; non-tiny ROOTs {non_tiny}/{root_files}",
        ha="center",
        va="top",
        fontsize=12.5,
        color="#333333",
    )

    for row, (var, xlabel, xlim) in enumerate(VARS):
        for col, (stage, stage_label) in enumerate(STAGES):
            ax = axes[row][col]
            ax.set_facecolor("white")
            factor = 4 if var in {"e11_to_e33", "weta_cogx"} else 2
            plotted = []
            for cent, cent_label, color in CENTRALITIES:
                arr = arrays(cache, var, cent, stage)
                if arr is None:
                    continue
                x, y, e = select_range(arr, xlim, factor)
                if x.size == 0:
                    continue
                ax.step(x, y, where="mid", color=color, linewidth=1.9, label=cent_label)
                ax.errorbar(
                    x,
                    y,
                    yerr=e,
                    linestyle="none",
                    marker="o",
                    markersize=3.2,
                    color=color,
                    elinewidth=0.8,
                    capsize=0,
                    alpha=0.9,
                )
                plotted.append(y)
            ax.set_xlim(*xlim)
            ax.set_ylim(bottom=0)
            if plotted:
                ymax = max(float(np.nanmax(y)) for y in plotted)
                ax.set_ylim(0, ymax * 1.20 if ymax > 0 else 1.0)
            ax.grid(True, axis="y", color="#e3e7ed", linewidth=0.8)
            ax.tick_params(direction="in", top=True, right=True, labelsize=10.5)
            if row == 0:
                ax.set_title(stage_label, fontsize=14, fontweight="bold", pad=10)
            if col == 0:
                ax.set_ylabel(f"{xlabel}\nnormalized counts", fontsize=12)
            if row == len(VARS) - 1:
                ax.set_xlabel(xlabel, fontsize=12)
            else:
                ax.set_xlabel("")
            add_sphenix_label(ax, include_details=(row == 0 and col == 0))
            if row == 0 and col == 2:
                ax.legend(loc="upper right", frameon=False, fontsize=10.5, handlelength=1.5)

    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.075, top=0.895, wspace=0.24, hspace=0.32)
    out = outdir / "the42_b002_interim_complete_runs_auau_shower_shape_stage_flow.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)

    manifest = {
        "schema": "THE42_B002_INTERIM_SHOWER_SHAPE_PLOT_MANIFEST_V1",
        "cache": str(DEFAULT_CACHE),
        "outputs": [str(out)],
        "source": "AuAu b002 completed-run non-tiny ROOT histogram cache; no active-output merge",
        "variables": [v[0] for v in VARS],
        "stages": [s[0] for s in STAGES],
        "centralities": [c[0] for c in CENTRALITIES],
        "metadata": meta,
    }
    manifest_path = outdir / "the42_b002_interim_complete_runs_auau_shower_shape_stage_flow_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return [out, manifest_path]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = ap.parse_args()
    cache = load_cache(args.cache)
    outputs = draw(cache, args.outdir)
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()

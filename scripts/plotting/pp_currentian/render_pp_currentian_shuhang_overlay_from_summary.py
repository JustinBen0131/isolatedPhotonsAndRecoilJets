#!/usr/bin/env python3
"""Render pp current-IAN BDT overlay using cached this-analysis score histograms."""

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
import json
from pathlib import Path

import numpy as np
import uproot


def unit(hist: np.ndarray) -> np.ndarray:
    total = float(np.sum(hist))
    if total <= 0.0:
        return hist
    return hist / total


def project_shuhang(path: Path, hist_name: str, edges: np.ndarray) -> dict:
    with uproot.open(path) as handle:
        if hist_name not in handle:
            raise SystemExit(f"Missing {hist_name} in {path}")
        values, xedges, yedges = handle[hist_name].to_numpy(flow=False)
    values = np.asarray(values, dtype="float64")
    xedges = np.asarray(xedges, dtype="float64")
    centers = 0.5 * (xedges[:-1] + xedges[1:])
    projection = values.sum(axis=1)
    use = np.isfinite(centers) & np.isfinite(projection) & (centers >= edges[0]) & (centers < edges[-1])
    rebinned = np.histogram(centers[use], bins=edges, weights=projection[use])[0]
    return {
        "hist": unit(rebinned),
        "raw_hist": rebinned.tolist(),
        "source_root": str(path),
        "hist_name": hist_name,
        "source_x_edges": xedges.tolist(),
        "source_y_edges": np.asarray(yedges, dtype="float64").tolist(),
        "source_integral": float(np.sum(projection)),
        "selected_integral": float(np.sum(rebinned)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--shuhang-signal-root", required=True, type=Path)
    parser.add_argument("--shuhang-inclusive-root", required=True, type=Path)
    parser.add_argument("--shuhang-hist", default="h2d_bdt_eta0_pt2_cut1")
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    summary = json.loads(args.summary.read_text())
    edges = np.asarray(summary["bins"], dtype="float64")
    centers = 0.5 * (edges[:-1] + edges[1:])
    signal = unit(np.asarray(summary["signal_raw_hist"], dtype="float64"))
    inclusive = unit(np.asarray(summary["inclusive_raw_hist"], dtype="float64"))
    shuhang_signal = project_shuhang(args.shuhang_signal_root, args.shuhang_hist, edges)
    shuhang_inclusive = project_shuhang(args.shuhang_inclusive_root, args.shuhang_hist, edges)

    args.outdir.mkdir(parents=True, exist_ok=True)

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.2, 5.7))
    ax.step(centers, signal, where="mid", color="#d62728", linewidth=2.1, label="This analysis signal")
    ax.step(centers, inclusive, where="mid", color="#1f5fff", linewidth=2.1, label="This analysis inclusive")
    ax.step(
        centers,
        shuhang_signal["hist"],
        where="mid",
        color="#8b1a1a",
        linewidth=1.9,
        linestyle="--",
        label="PPG12/Shuhang signal",
    )
    ax.step(
        centers,
        shuhang_inclusive["hist"],
        where="mid",
        color="#173a9a",
        linewidth=1.9,
        linestyle="--",
        label="PPG12/Shuhang inclusive",
    )

    ax.text(0.58, 0.92, "sPHENIX Internal", transform=ax.transAxes, fontsize=14, ha="left")
    ax.text(0.58, 0.865, r"$p$+$p$ $\sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=12.5, ha="left")
    ax.text(0.58, 0.812, r"$|\eta|<0.7$, $18<E_T<22$ GeV", transform=ax.transAxes, fontsize=12.5, ha="left")
    ax.text(0.58, 0.762, "current-IAN preselection + NPB > 0.5", transform=ax.transAxes, fontsize=12.0, ha="left")

    ax.set_xlabel("BDT score", fontsize=15)
    ax.set_ylabel("unit-normalized counts", fontsize=15)
    ax.set_xlim(0.0, 1.0)
    ymax = max(float(np.max(signal)), float(np.max(inclusive)), float(np.max(shuhang_signal["hist"])), float(np.max(shuhang_inclusive["hist"])), 0.05)
    ax.set_ylim(0.0, ymax * 1.18)
    ax.tick_params(direction="in", top=True, right=True, labelsize=12)
    ax.legend(frameon=False, fontsize=11, loc="upper left", bbox_to_anchor=(0.02, 0.98), ncol=1)
    fig.tight_layout()

    png = args.outdir / "pp_currentian_basev3e_bdt_score_overlay_vs_shuhang_clean.png"
    fig.savefig(png, dpi=240)
    plt.close(fig)

    out_summary = {
        "plot": str(png),
        "source_summary": str(args.summary),
        "normalization": "each curve divided by its own bin sum",
        "this_analysis_signal_rows_after_cuts": summary.get("signal", {}).get("rows_after_cuts"),
        "this_analysis_inclusive_rows_after_cuts": summary.get("inclusive", {}).get("rows_after_cuts"),
        "shuhang_signal": {k: v for k, v in shuhang_signal.items() if k != "hist"},
        "shuhang_inclusive": {k: v for k, v in shuhang_inclusive.items() if k != "hist"},
        "bins": edges.tolist(),
        "this_analysis_signal_hist": signal.tolist(),
        "this_analysis_inclusive_hist": inclusive.tolist(),
        "shuhang_signal_hist": shuhang_signal["hist"].tolist(),
        "shuhang_inclusive_hist": shuhang_inclusive["hist"].tolist(),
    }
    (args.outdir / "pp_currentian_basev3e_bdt_score_overlay_vs_shuhang_clean_summary.json").write_text(
        json.dumps(out_summary, indent=2) + "\n"
    )
    print(png)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

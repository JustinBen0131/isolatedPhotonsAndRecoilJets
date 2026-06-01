#!/usr/bin/env python3
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

import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_variant_csv(path: Path) -> dict[tuple[str, str], tuple[np.ndarray, np.ndarray]]:
    rows: dict[tuple[str, str], list[tuple[float, float, float]]] = defaultdict(list)
    with path.open() as handle:
        reader = csv.DictReader(line for line in handle if line.startswith(("sample,", "signal,", "inclusive,")))
        for row in reader:
            key = (row["sample"], f"pt{row['pt']} cut{row['cut']} {row['slice']}")
            rows[key].append((float(row["lo"]), float(row["hi"]), float(row["value"])))
    out = {}
    for key, vals in rows.items():
        vals = sorted(vals)
        centers = np.asarray([(lo + hi) * 0.5 for lo, hi, _ in vals])
        hist = np.asarray([v for _, _, v in vals])
        out[key] = (centers, hist)
    return out


def main() -> None:
    root = Path("/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_finalShuhang_20260526_2115/validation/fullsim_shuhang_overlay")
    summary = json.loads((root / "pp_currentian_basev3e_bdt_score_overlay_vs_shuhang_clean_summary.json").read_text())
    variants = read_variant_csv(root / "shuhang_projection_variant_scan.csv")

    bins = np.asarray(summary["bins"], dtype=float)
    centers = 0.5 * (bins[:-1] + bins[1:])
    ours = {
        "signal": np.asarray(summary["this_analysis_signal_hist"], dtype=float),
        "inclusive": np.asarray(summary["this_analysis_inclusive_hist"], dtype=float),
    }

    plot_variants = ["pt2 cut0 allY", "pt2 cut1 allY", "pt2 cut2 allY", "pt2 cut1 Ylt0", "pt2 cut1 Ylt5", "pt2 cut1 Y0to5"]
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), dpi=180, sharey=True)
    colors = ["#8b1a1a", "#1f4cc9", "#2a8c2a", "#7a3dd8", "#e07000", "#555555"]
    for ax, sample, title in [(axes[0], "signal", "Signal"), (axes[1], "inclusive", "Inclusive / background")]:
        ax.step(centers, ours[sample], where="mid", color="black", linewidth=2.6, label="This analysis")
        for variant, color in zip(plot_variants, colors):
            key = (sample, variant)
            if key not in variants:
                continue
            x, h = variants[key]
            ax.step(x, h, where="mid", linewidth=1.55, color=color, label=f"Shuhang {variant}")
        ax.set_title(title, fontsize=15, fontweight="bold")
        ax.set_xlabel("BDT score", fontsize=13)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 0.62)
        ax.grid(True, alpha=0.25)
    axes[0].set_ylabel("unit-normalized counts", fontsize=13)
    axes[1].legend(loc="upper center", bbox_to_anchor=(0.50, 1.03), ncol=2, fontsize=8.5, frameon=False)
    fig.suptitle("Shuhang Histogram Choice / Projection Diagnostic", fontsize=18, fontweight="bold")
    fig.text(
        0.5,
        0.02,
        "All curves use Shuhang's provided ROOT histograms only as references; variants test pt/cut selection and second-axis slicing.",
        ha="center",
        fontsize=11,
    )
    fig.tight_layout(rect=[0.02, 0.06, 0.98, 0.92])
    out = root / "pp_currentian_shuhang_projection_variant_diagnostic.png"
    fig.savefig(out)
    print(out)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Scan all available 0-20% b009 AuAu shower-shape variables."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[5])
SCRIPT_DIR = THIS_FILE.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.append(str(SCRIPT_DIR))

from make_the42_b009_pp_overlay_three_shape_slide import (  # noqa: E402
    DATA_ROOT,
    INCLUSIVE_ROOT,
    SIGNAL_ROOT,
    OUTDIR,
    add_sphenix_internal,
    build_index,
    draw_line,
    draw_open_markers,
    load_aa_curve,
    setup_style,
)


VARIABLES = [
    ("weta", r"$w_{\eta}^{\mathrm{cogX}}$", (0.0, 1.20), 2),
    ("wphi", r"$w_{\phi}^{\mathrm{cogX}}$", (0.0, 1.20), 2),
    ("e11e33", r"$E_{11}/E_{33}$", (0.0, 1.10), 2),
    ("e32e35", r"$E_{32}/E_{35}$", (0.0, 1.10), 2),
    ("et1", r"$e_{T1}$", (0.0, 1.10), 2),
    ("weta33", r"$w_{\eta,33}$", (0.0, 1.20), 2),
    ("weta35", r"$w_{\eta,35}$", (0.0, 1.20), 2),
    ("wphi33", r"$w_{\phi,33}$", (0.0, 1.20), 2),
    ("wphi53", r"$w_{\phi,53}$", (0.0, 1.20), 2),
]
STAGES = [("After preselection", "pre"), ("After tight ID", "tight")]


def render() -> tuple[Path, Path]:
    setup_style()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    indexes = {
        "data": build_index(DATA_ROOT),
        "signal": build_index(SIGNAL_ROOT),
        "inclusive": build_index(INCLUSIVE_ROOT),
    }
    colors = {"signal": "#B91C1C", "inclusive": "#2E63D4", "auau": "#111827"}

    fig, axes = plt.subplots(len(VARIABLES), 2, figsize=(18.0, 31.0), constrained_layout=False)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.925, bottom=0.035, wspace=0.120, hspace=0.360)

    manifest = {
        "output": None,
        "purpose": "0-20% variable scan sheet for choosing shower-shape variables",
        "data_root": str(DATA_ROOT),
        "signal_root": str(SIGNAL_ROOT),
        "inclusive_root": str(INCLUSIVE_ROOT),
        "centrality": "0-20%",
        "curves": [],
    }

    for row, (var, label, x_range, rebin_factor) in enumerate(VARIABLES):
        for col, (stage_label, stage_key) in enumerate(STAGES):
            ax = axes[row, col]
            ax.set_facecolor("#FFFFFF")
            ax.grid(True, axis="y", color="#E5E7EB", linewidth=0.75, alpha=0.95)
            ax.grid(True, axis="x", color="#EEF2F7", linewidth=0.55, alpha=0.70)
            curves = {
                "signal": load_aa_curve(SIGNAL_ROOT, indexes["signal"], "signal MC", f"{stage_key}_sig", var, x_range, rebin_factor),
                "inclusive": load_aa_curve(
                    INCLUSIVE_ROOT,
                    indexes["inclusive"],
                    "inclusive MC",
                    f"{stage_key}_bkg",
                    var,
                    x_range,
                    rebin_factor,
                ),
                "auau": load_aa_curve(DATA_ROOT, indexes["data"], "AuAu data", stage_key, var, x_range, rebin_factor),
            }
            draw_line(ax, curves["signal"], colors["signal"])
            draw_line(ax, curves["inclusive"], colors["inclusive"])
            draw_open_markers(ax, curves["auau"], colors["auau"], "o", 4.8, 6)
            max_y = max(float(c.values.max()) if c.values.size else 0.0 for c in curves.values())
            ax.set_xlim(*x_range)
            ax.set_ylim(0.0, max_y * 1.22 if max_y > 0 else 1.0)
            ax.tick_params(direction="in", top=True, right=True, length=3.5, width=0.72, labelsize=8.0)
            ax.set_xlabel(label, fontsize=10.5, labelpad=1.5)
            if col == 0:
                ax.set_ylabel(f"{label}\nunit-normalized", fontsize=9.8)
            if row == 0:
                ax.set_title(stage_label, fontsize=15.0, fontweight="bold", pad=8)
            if row == 0 and col == 0:
                add_sphenix_internal(ax)
            for key, curve in curves.items():
                manifest["curves"].append(
                    {
                        "variable": var,
                        "stage": stage_key,
                        "sample": key,
                        "integral_before_normalization": curve.integral,
                        "matched_histograms": len(curve.matches),
                        "missing_histogram_count": len(curve.missing),
                        "missing_histograms": curve.missing[:20],
                    }
                )

    handles = [
        Line2D([0], [0], color=colors["signal"], lw=2.6, label="Signal MC"),
        Line2D([0], [0], color=colors["inclusive"], lw=2.6, label="Inclusive MC"),
        Line2D([0], [0], marker="o", color=colors["auau"], mfc="white", mec=colors["auau"], mew=1.5, lw=0, ms=6.2, label="AuAu data"),
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.985, 0.973), ncol=3, frameon=True, fancybox=False, edgecolor="#D1D5DB", facecolor="#F9FAFB", fontsize=12.0)
    fig.text(
        0.075,
        0.982,
        "0-20% AuAu shower-shape variable scan",
        ha="left",
        va="top",
        fontsize=24.0,
        fontweight="bold",
        color="#111827",
    )
    fig.text(
        0.075,
        0.958,
        r"Photon $p_T$: 15-35 GeV; left is after preselection, right is after tight ID. Curves are unit-normalized.",
        ha="left",
        va="top",
        fontsize=12.6,
        color="#374151",
    )

    out = OUTDIR / "the42_b009_0_20_all_shower_shapes_scan_sheet.png"
    manifest_path = OUTDIR / "the42_b009_0_20_all_shower_shapes_scan_sheet_manifest.json"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    manifest["output"] = str(out)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return out, manifest_path


def main() -> int:
    out, manifest = render()
    print(out)
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

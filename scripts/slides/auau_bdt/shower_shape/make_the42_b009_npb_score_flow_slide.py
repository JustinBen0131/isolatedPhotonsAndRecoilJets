#!/usr/bin/env python3
"""0-20% b009 AuAu NPB-score flow before/after NPB preselection."""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[5])
SCRIPT_DIR = THIS_FILE.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.append(str(SCRIPT_DIR))
SCRIPTS_DIR = REPO / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.append(str(SCRIPTS_DIR))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize  # noqa: E402
from make_the42_b009_pp_overlay_three_shape_slide import (  # noqa: E402
    DATA_ROOT,
    INCLUSIVE_ROOT,
    SIGNAL_ROOT,
    OUTDIR,
    add_sphenix_internal,
    build_index,
    draw_line,
    draw_open_markers,
    open_root,
    setup_style,
)


CENT_BINS_0_20 = [(0, 10), (10, 20)]
PT_GROUPS = [
    ("15-18 GeV", [(15, 18)]),
    ("18-24 GeV", [(18, 20), (20, 22), (22, 24)]),
    ("24-35 GeV", [(24, 26), (26, 28), (28, 30), (30, 35)]),
]
STAGES = [
    ("Before NPB preselection", "inclusive"),
    ("After NPB preselection", "npbPass"),
]


@dataclass
class Curve:
    sample: str
    centers: np.ndarray
    edges: np.ndarray
    values: np.ndarray
    errors: np.ndarray
    integral: float
    matches: list[str]
    missing: list[str]


def rebin(hist, factor: int):
    if factor <= 1:
        return hist
    if hist.GetNbinsX() % factor != 0:
        return hist
    out = hist.Rebin(factor, f"{hist.GetName()}_rebin{factor}")
    out.SetDirectory(0)
    return out


def hist_to_curve(hist, sample: str, x_range: tuple[float, float], rebin_factor: int, matches: list[str], missing: list[str]) -> Curve:
    hist = rebin(hist, rebin_factor)
    nb = hist.GetNbinsX()
    edges = np.array([hist.GetXaxis().GetBinLowEdge(i) for i in range(1, nb + 2)], dtype=float)
    centers = np.array([hist.GetXaxis().GetBinCenter(i) for i in range(1, nb + 1)], dtype=float)
    counts = np.array([hist.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    errors = np.array([hist.GetBinError(i) for i in range(1, nb + 1)], dtype=float)
    mask = (centers >= x_range[0]) & (centers <= x_range[1])
    keep = np.flatnonzero(mask)
    if len(keep):
        edges = np.concatenate(([edges[keep[0]]], edges[keep + 1]))
    centers = centers[mask]
    counts = counts[mask]
    errors = errors[mask]
    integral = float(np.sum(counts))
    values = np.zeros_like(counts)
    scaled_errors = np.zeros_like(errors)
    if integral > 0:
        values = counts / integral
        scaled_errors = errors / integral
    return Curve(sample, centers, edges, values, scaled_errors, integral, matches, missing)


def target_name(tag: str, pt: tuple[int, int], cent: tuple[int, int]) -> str:
    return f"h_ss_npbScore_{tag}_pT_{pt[0]}_{pt[1]}_cent_{cent[0]}_{cent[1]}"


def load_curve(
    root_path: Path,
    index: dict[str, list[str]],
    sample: str,
    tag: str,
    pt_bins: list[tuple[int, int]],
    x_range: tuple[float, float] = (0.0, 1.0),
    rebin_factor: int = 2,
) -> Curve:
    f = open_root(root_path)
    acc = None
    matches: list[str] = []
    missing: list[str] = []
    try:
        for pt in pt_bins:
            for cent in CENT_BINS_0_20:
                name = target_name(tag, pt, cent)
                paths = index.get(name, [])
                if not paths:
                    missing.append(name)
                    continue
                obj = f.Get(paths[0])
                if not obj or not obj.InheritsFrom("TH1"):
                    missing.append(paths[0])
                    continue
                h = obj.Clone(f"{name}_{sample}_{tag}")
                h.SetDirectory(0)
                matches.append(paths[0])
                if acc is None:
                    acc = h.Clone(f"npbScore_{sample}_{tag}_sum")
                    acc.SetDirectory(0)
                else:
                    acc.Add(h)
    finally:
        f.Close()
    if acc is None:
        raise RuntimeError(f"no NPB-score histograms for {sample} {tag}")
    return hist_to_curve(acc, sample, x_range, rebin_factor, matches, missing)


def render() -> tuple[Path, Path, Path]:
    setup_style()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    indexes = {
        "data": build_index(DATA_ROOT),
        "signal": build_index(SIGNAL_ROOT),
        "inclusive": build_index(INCLUSIVE_ROOT),
    }
    colors = {"signal": "#B91C1C", "inclusive": "#2E63D4", "auau": "#111827"}
    fig, axes = plt.subplots(3, 2, figsize=slide_figsize(), constrained_layout=False)
    fig.subplots_adjust(left=0.076, right=0.985, top=0.790, bottom=0.090, wspace=0.110, hspace=0.340)

    manifest = {
        "output": None,
        "data_root": str(DATA_ROOT),
        "signal_root": str(SIGNAL_ROOT),
        "inclusive_root": str(INCLUSIVE_ROOT),
        "centrality": "0-20%",
        "pt_groups": [{"label": label, "bins": bins} for label, bins in PT_GROUPS],
        "stages": STAGES,
        "note": "This uses stored h_ss_npbScore_inclusive and h_ss_npbScore_npbPass histograms; final tight-BDT score histograms were not present in the merged b009 ROOT.",
        "curves": [],
    }

    for row, (pt_label, pt_bins) in enumerate(PT_GROUPS):
        for col, (stage_label, tag) in enumerate(STAGES):
            ax = axes[row, col]
            ax.set_facecolor("#FFFFFF")
            ax.grid(True, axis="y", color="#E5E7EB", linewidth=0.75, alpha=0.95)
            ax.grid(True, axis="x", color="#EEF2F7", linewidth=0.55, alpha=0.70)
            curves = {
                "signal": load_curve(SIGNAL_ROOT, indexes["signal"], "signal MC", tag, pt_bins),
                "inclusive": load_curve(INCLUSIVE_ROOT, indexes["inclusive"], "inclusive MC", tag, pt_bins),
                "auau": load_curve(DATA_ROOT, indexes["data"], "AuAu data", tag, pt_bins),
            }
            draw_line(ax, curves["signal"], colors["signal"])
            draw_line(ax, curves["inclusive"], colors["inclusive"])
            draw_open_markers(ax, curves["auau"], colors["auau"], "o", 5.4, 6)
            ax.axvline(0.5, color="#6B7280", lw=1.2, linestyle="--", alpha=0.86)
            max_y = max(float(c.values.max()) if c.values.size else 0.0 for c in curves.values())
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(0.0, max_y * 1.30 if max_y > 0 else 1.0)
            ax.tick_params(direction="in", top=True, right=True, length=4.0, width=0.75)
            ax.set_xlabel("NPB score", fontsize=11.2, labelpad=2.0)
            if col == 0:
                ax.set_ylabel(f"{pt_label}\nunit-normalized", fontsize=10.2)
            if row == 0:
                ax.set_title(stage_label, fontsize=15.2, fontweight="bold", pad=10)
            ax.text(
                0.970,
                0.890,
                pt_label,
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=10.9,
                fontweight="bold",
                bbox={"boxstyle": "round,pad=0.22", "facecolor": "#F9FAFB", "edgecolor": "#E5E7EB", "linewidth": 0.6},
            )
            if row == 0 and col == 0:
                add_sphenix_internal(ax)
            for key, curve in curves.items():
                manifest["curves"].append(
                    {
                        "pt_group": pt_label,
                        "stage": tag,
                        "sample": key,
                        "integral_before_normalization": curve.integral,
                        "matched_histograms": len(curve.matches),
                        "missing_histogram_count": len(curve.missing),
                        "missing_histograms": curve.missing[:20],
                    }
                )

    handles = [
        Line2D([0], [0], color=colors["signal"], lw=2.7, label="Signal MC"),
        Line2D([0], [0], color=colors["inclusive"], lw=2.7, label="Inclusive MC"),
        Line2D([0], [0], marker="o", color=colors["auau"], mfc="white", mec=colors["auau"], mew=1.6, lw=0, ms=6.8, label="AuAu data"),
        Line2D([0], [0], color="#6B7280", lw=1.3, linestyle="--", label="NPB cut"),
    ]
    fig.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(0.985, 0.910),
        ncol=4,
        frameon=True,
        fancybox=False,
        framealpha=0.94,
        edgecolor="#D1D5DB",
        facecolor="#F9FAFB",
        fontsize=11.2,
        handlelength=2.2,
        columnspacing=1.2,
    )
    fig.text(
        0.075,
        0.955,
        "0-20% AuAu BDT-score flow through photon preselection",
        ha="left",
        va="top",
        fontsize=24.2,
        fontweight="bold",
        color="#111827",
    )
    fig.text(
        0.075,
        0.912,
        r"NPB score, photon $p_T$: 15-35 GeV split into three rows. Curves are unit-normalized in each panel.",
        ha="left",
        va="top",
        fontsize=13.0,
        color="#374151",
    )
    fig.text(
        0.075,
        0.878,
        "Left: all scored candidates before the NPB cut. Right: candidates passing the NPB-gated preselection.",
        ha="left",
        va="top",
        fontsize=12.1,
        color="#4B5563",
    )

    out = OUTDIR / "the42_b009_0_20_npb_score_preselection_flow_slide.png"
    manifest_path = OUTDIR / "the42_b009_0_20_npb_score_preselection_flow_manifest.json"
    script_path = OUTDIR / "the42_b009_0_20_npb_score_preselection_flow_speaker_script.md"
    fig.savefig(out, dpi=SLIDE_DPI)
    plt.close(fig)
    manifest["output"] = str(out)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    script_path.write_text(
        "\n".join(
            [
                "# Speaker Notes",
                "",
                "This slide shows the NPB-score distribution in 0-20% AuAu before and after the NPB-gated photon preselection.",
                "",
                "Each row is a photon-pT range. The left column is all scored candidates before applying the NPB cut. The right column is the NPB-pass population after the preselection gate. The dashed line marks the NPB cut at 0.5.",
                "",
                "The plot is intentionally unit-normalized so the visual comparison is the BDT-score shape and not the raw sample size.",
                "",
                "Important caveat: the merged b009 ROOT does not contain the final tight-BDT score distribution under the expected h_tightBDTScore_preselected names, so this is specifically the NPB-preselection score flow rather than a final tight-BDT-score flow.",
                "",
            ]
        )
    )
    return out, manifest_path, script_path


def main() -> int:
    out, manifest, script = render()
    print(out)
    print(manifest)
    print(script)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

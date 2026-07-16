#!/usr/bin/env python3
"""Render the THE-95 versus THE-101 centrality score-shape slide."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.ticker import MaxNLocator
import numpy as np


THIS_FILE = Path(__file__).resolve()
REPO = THIS_FILE.parents[4]
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize


DEFAULT_INPUT = (
    REPO
    / "dataOutput/auau/canonical_bdt_the101_20260714/slides"
    / "the101_minbias_before_after_score_payload.json"
)
DEFAULT_OUTPUT = (
    REPO
    / "dataOutput/auau/canonical_bdt_the101_20260714/slides"
    / "auau_bdt_before_after_minbias_classifier.png"
)

INK = "#111827"
MUTED = "#526071"
GRID = "#dce5f0"
EDGE = "#b9c6d8"
RED = "#d13b35"
BLUE = "#1f78b4"
GREEN = "#1b8a5a"
GOLD = "#c88a00"
ROW_COLORS = ("#8b5a2b", "#147d64")
FONT = "Times New Roman"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def configure_style() -> None:
    for name in (
        "Times New Roman.ttf",
        "Times New Roman Bold.ttf",
        "Times New Roman Italic.ttf",
        "Times New Roman Bold Italic.ttf",
    ):
        path = Path("/System/Library/Fonts/Supplemental") / name
        if path.exists():
            font_manager.fontManager.addfont(str(path))
    plt.rcParams.update(
        {
            "font.family": FONT,
            "font.serif": [FONT, "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.unicode_minus": False,
            "savefig.facecolor": "white",
        }
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def step_xy(edges: list[float], density: list[float]) -> tuple[list[float], list[float]]:
    x: list[float] = []
    y: list[float] = []
    for idx, value in enumerate(density):
        x.extend((edges[idx], edges[idx + 1]))
        y.extend((value, value))
    return x, y


def fmt_count(value: int) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f}M"
    if value >= 1000:
        return f"{value / 1000:.0f}k"
    return str(value)


def draw_metrics(ax: plt.Axes, metrics: dict) -> None:
    ax.add_patch(
        Rectangle(
            (0.018, 0.625),
            0.765,
            0.350,
            transform=ax.transAxes,
            facecolor="white",
            edgecolor=EDGE,
            linewidth=0.8,
            zorder=6,
        )
    )
    rows = (
        ("AUC", f"{metrics['auc']:.3f}", "Gap", f"{metrics['median_gap']:.2f}"),
        ("Logloss", f"{metrics['logloss']:.3f}", "Fake rate", f"{100.0 * metrics['wp80_inclusive_fake']:.1f}%"),
    )
    for row, values in enumerate(rows):
        y = 0.876 - 0.174 * row
        for text, x, ha, bold, size in (
            (values[0], 0.038, "left", False, 8.2),
            (values[1], 0.360, "right", True, 9.4),
            (values[2], 0.420, "left", False, 7.8),
            (values[3], 0.755, "right", True, 9.4),
        ):
            ax.text(
                x,
                y,
                text,
                transform=ax.transAxes,
                ha=ha,
                va="center",
                fontsize=size,
                fontweight="bold" if bold else "normal",
                color=INK,
                zorder=7,
            )


def draw_pad(
    ax: plt.Axes,
    cell: dict,
    edges: list[float],
    ymax: float,
    *,
    show_y: bool,
    show_x: bool,
) -> None:
    sx, sy = step_xy(edges, cell["signal"]["density"])
    ix, iy = step_xy(edges, cell["inclusive"]["density"])
    ax.plot(sx, sy, color=RED, linewidth=2.1, zorder=3)
    ax.fill_between(sx, sy, step="pre", color=RED, alpha=0.08, zorder=2)
    ax.plot(ix, iy, color=BLUE, linewidth=2.1, zorder=3)
    ax.fill_between(ix, iy, step="pre", color=BLUE, alpha=0.08, zorder=2)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, ymax)
    ax.grid(True, color=GRID, linewidth=0.7)
    ax.tick_params(
        axis="both",
        direction="in",
        top=True,
        right=True,
        labelsize=8.7,
        length=4,
        width=0.8,
        colors="#374151",
    )
    ax.yaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
    if show_y:
        ax.set_ylabel("Unit-area density", fontsize=10.3, labelpad=3)
    else:
        ax.set_yticklabels([])
    if show_x:
        ax.set_xlabel("BDT score", fontsize=10.7, labelpad=2)
    else:
        ax.set_xticklabels([])
    for spine in ax.spines.values():
        spine.set_color("#4b5563")
        spine.set_linewidth(0.9)
    draw_metrics(ax, cell["metrics"])
    ax.text(
        0.975,
        0.945,
        f"S {fmt_count(cell['signal']['entries'])}\nI {fmt_count(cell['inclusive']['entries'])}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.3,
        color=MUTED,
        fontweight="bold",
        linespacing=0.88,
        zorder=7,
    )


def draw_auc_table(fig: plt.Figure, payload: dict) -> list[dict]:
    left, bottom, width, height = 0.535, 0.540, 0.420, 0.235
    fig.text(
        left,
        bottom + height + 0.014,
        r"AUC by $p_T$ and centrality (before $\rightarrow$ after)",
        ha="left",
        va="bottom",
        fontsize=18.0,
        fontweight="bold",
        color=INK,
    )
    ax = fig.add_axes([left, bottom, width, height])
    ax.set_axis_off()
    centralities = ((0.0, 20.0), (20.0, 50.0), (50.0, 80.0))
    pts = ((15.0, 20.0), (20.0, 25.0), (25.0, 35.0))
    before = payload["models"]["before"]["auc_by_pt_centrality"]
    after = payload["models"]["after"]["auc_by_pt_centrality"]
    col_w = 1.0 / 3.0
    row_h = 0.25
    for col, (clo, chi) in enumerate(centralities):
        x = col * col_w
        ax.add_patch(Rectangle((x, 0.75), col_w, 0.25, facecolor="#eef3f8", edgecolor=EDGE, linewidth=1.0))
        ax.text(x + col_w / 2, 0.875, f"{clo:g}-{chi:g}%", ha="center", va="center", fontsize=12.2, fontweight="bold", color=INK)
    for row, (plo, phi) in enumerate(pts):
        y = 0.50 - row * row_h
        for col, (clo, chi) in enumerate(centralities):
            x = col * col_w
            key = f"pt_{plo:g}_{phi:g}_cent_{clo:g}_{chi:g}"
            b = float(before[key])
            a = float(after[key])
            delta = a - b
            display_delta = 0.0 if abs(delta) < 5.0e-5 else delta
            fill = "#eff8f3" if display_delta >= 0 else "#fff3ed"
            ax.add_patch(Rectangle((x, y), col_w, row_h, facecolor=fill, edgecolor=EDGE, linewidth=1.0))
            ax.text(x + 0.018, y + row_h * 0.73, f"{plo:g}-{phi:g} GeV", ha="left", va="center", fontsize=9.8, color=MUTED)
            ax.text(x + col_w / 2, y + row_h * 0.40, f"{b:.3f} → {a:.3f}", ha="center", va="center", fontsize=11.6, fontweight="bold", color=INK)
            ax.text(
                x + col_w - 0.016,
                y + row_h * 0.14,
                f"Δ {display_delta:+.4f}",
                ha="right",
                va="center",
                fontsize=8.8,
                color=GREEN if display_delta >= 0 else "#b42318",
            )
    return [{"name": "auc_summary", "kind": "plot", "role": "plot", "bbox": [1370, 302, 2445, 662]}]


def write_csv(payload: dict, path: Path) -> None:
    before = payload["models"]["before"]["auc_by_pt_centrality"]
    after = payload["models"]["after"]["auc_by_pt_centrality"]
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["pt_min", "pt_max", "centrality_min", "centrality_max", "before_auc", "after_auc", "delta_auc"])
        for plo, phi in payload["pt_bins"]:
            for clo, chi in payload["centrality_bins"]:
                key = f"pt_{plo:g}_{phi:g}_cent_{clo:g}_{chi:g}"
                writer.writerow([plo, phi, clo, chi, before[key], after[key], float(after[key]) - float(before[key])])


def main() -> int:
    args = parse_args()
    payload = json.loads(args.input.read_text())
    if payload.get("status") != "READY":
        raise ValueError("score payload is not READY")
    configure_style()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")

    fig.text(
        0.043,
        0.958,
        "Canonical Au+Au photon-ID BDT: minimum-bias classifier comparison",
        ha="left",
        va="top",
        fontsize=28.5,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.043,
        0.878,
        "Same 14-feature contract, six embedded samples, 15-35 GeV and PPG12 weighting; neither model uses a low-energy pretraining cut",
        ha="left",
        va="top",
        fontsize=14.7,
        color=MUTED,
    )

    bullets = (
        ("Before", "unrestricted embedded events; no low-energy cut"),
        ("After", r"require isAuAuMinimumBias() before training-tree filling"),
        ("Comparison", "identical classifier-pass validation sample"),
    )
    y = 0.776
    for idx, (head, body) in enumerate(bullets):
        yy = y - idx * 0.063
        fig.text(0.047, yy, "▸", ha="left", va="center", fontsize=21, color=ROW_COLORS[min(idx, 1)], fontfamily="DejaVu Sans")
        fig.text(0.068, yy, f"{head}:", ha="left", va="center", fontsize=16.2, fontweight="bold", color=INK)
        fig.text(0.139 if idx < 2 else 0.178, yy, body, ha="left", va="center", fontsize=14.6, color=INK)

    layout_nodes = [
        {"name": "title", "kind": "text", "role": "title", "text": "Canonical Au+Au photon-ID BDT: minimum-bias classifier comparison", "font_px": 79, "bbox": [110, 52, 2445, 145], "title_anchor": True, "colon_style_exception": True},
        {"name": "subtitle", "kind": "text", "role": "audience", "text": "Same 14-feature contract, six embedded samples, 15-35 GeV and PPG12 weighting; neither model uses a low-energy pretraining cut", "font_px": 41, "bbox": [110, 178, 2445, 226], "title_axis_align": "left"},
        {"name": "before_contract", "kind": "text", "role": "audience", "text": "Before: unrestricted embedded events; no low-energy cut", "font_px": 40, "bbox": [174, 300, 1260, 347], "text_runs": [{"text": "Before:", "bold": True}, {"text": " unrestricted embedded events; no low-energy cut", "bold": False}]},
        {"name": "after_contract", "kind": "text", "role": "audience", "text": "After: require isAuAuMinimumBias() before training-tree filling", "font_px": 40, "bbox": [174, 391, 1325, 438], "text_runs": [{"text": "After:", "bold": True}, {"text": " require isAuAuMinimumBias() before training-tree filling", "bold": False}]},
        {"name": "comparison_contract", "kind": "text", "role": "audience", "text": "Comparison: identical classifier-pass validation sample", "font_px": 40, "bbox": [174, 482, 1325, 529], "text_runs": [{"text": "Comparison:", "bold": True}, {"text": " identical classifier-pass validation sample", "bold": False}]},
    ]
    layout_nodes.extend(draw_auc_table(fig, payload))

    centralities = (("0_20", "0-20%"), ("20_50", "20-50%"), ("50_80", "50-80%"))
    axes_left = (0.115, 0.397, 0.679)
    ax_width = 0.260
    ax_height = 0.188
    row_bottoms = (0.287, 0.060)
    all_densities = []
    for model in ("before", "after"):
        for key, _ in centralities:
            cell = payload["models"][model]["cells"][key]
            all_densities.extend(cell["signal"]["density"])
            all_densities.extend(cell["inclusive"]["density"])
    ymax = max(all_densities) * 1.10
    edges = payload["bin_edges"]
    for row, model in enumerate(("before", "after")):
        bottom = row_bottoms[row]
        for col, (key, title) in enumerate(centralities):
            ax = fig.add_axes([axes_left[col], bottom, ax_width, ax_height])
            draw_pad(
                ax,
                payload["models"][model]["cells"][key],
                edges,
                ymax,
                show_y=col == 0,
                show_x=row == 1,
            )
            if row == 0:
                ax.set_title(title, fontsize=15.2, fontweight="bold", pad=7, color=INK)
            layout_nodes.append(
                {
                    "name": f"{model}_{key}",
                    "kind": "plot",
                    "role": "plot",
                    "bbox": [int(axes_left[col] * 2560), int((1.0 - bottom - ax_height) * 1440), int((axes_left[col] + ax_width) * 2560), int((1.0 - bottom) * 1440)],
                }
            )

    for row, (label, color) in enumerate(
        (
            ("Before", ROW_COLORS[0]),
            ("After", ROW_COLORS[1]),
        )
    ):
        center = row_bottoms[row] + ax_height / 2
        fig.patches.append(Rectangle((0.008, row_bottoms[row] + 0.020), 0.006, ax_height - 0.040, transform=fig.transFigure, facecolor=color, edgecolor=color))
        fig.text(0.019, center, label, ha="left", va="center", fontsize=16.8, fontweight="bold", color=INK)

    fig.text(0.255, 0.525, "Signal MC: embedded Photon12+20", ha="right", va="center", fontsize=11.8, fontweight="bold", color=RED)
    fig.text(0.275, 0.525, "Inclusive MC: Jet12+20+30+40, no truth filter", ha="left", va="center", fontsize=11.8, fontweight="bold", color=BLUE)

    fig.savefig(output, dpi=SLIDE_DPI, facecolor="white")
    plt.close(fig)

    layout = output.with_suffix(".layout_nodes.json")
    manifest = output.with_suffix(".manifest.json")
    csv_path = output.with_suffix(".csv")
    speaker = output.with_name(output.stem + "_speaker_script.md")
    audit = output.with_suffix(".audit.json")
    write_csv(payload, csv_path)
    layout.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "title_axis_x": 110,
                "minimum_audience_font_px": 40,
                "minimum_title_font_px": 72,
                "minimum_plot_annotation_font_px": 22,
                "nodes": layout_nodes,
            },
            indent=2,
        )
        + "\n"
    )
    speaker.write_text(
        "\n".join(
            (
                "# Au+Au minimum-bias classifier BDT comparison",
                "",
                "The upper-left text states the controlled comparison: both models use the same 14 features, six embedded samples, 15-35 GeV window, and PPG12 display weighting, with no low-energy pretraining veto.",
                "The only training-policy change is the embedded MinimumBiasClassifier requirement before training-tree filling.",
                "For a fair model comparison, both models are evaluated on the identical classifier-pass extraction.",
                "The lower six panels show source-defined truth-photon signal versus the canonical four-sample inclusive-jet background in the three centrality intervals.",
                "The score shapes and pad-level AUC, logloss, median gap, and WP80 fake rate are nearly unchanged.",
                "The upper-right table checks the same conclusion in three photon-pT intervals inside each centrality interval; every AUC change is small.",
                "The conclusion is that the classifier removes the anomalous embedded event population before training without materially changing the photon-ID ranking on accepted events.",
                "",
            )
        )
    )
    manifest.write_text(
        json.dumps(
            {
                "schema": "AUAU_BDT_MINBIAS_BEFORE_AFTER_SLIDE_V1",
                "output_png": str(output),
                "input_payload": {"path": str(args.input.resolve()), "sha256": sha256(args.input.resolve())},
                "models": {
                    "before": "THE-95 no-veto baseline14, trained without MinimumBiasClassifier requirement",
                    "after": "THE-101 canonical baseline14, trained with MinimumBiasClassifier pass requirement",
                },
                "comparison_scope": payload["comparison_scope"],
                "selection": payload["selection"],
                "plot_class_definition": payload["plot_class_definition"],
                "weighting": payload["weight_report"]["weight_mode"],
                "layout_nodes": str(layout),
                "auc_csv": str(csv_path),
                "speaker_script": str(speaker),
                "audit": str(audit),
                "source_slide": "Google Slides object g3efdf2cd1ae_0_446",
            },
            indent=2,
        )
        + "\n"
    )
    audit_cmd = [
        sys.executable,
        str(REPO / "scripts/slides/common/post_render_slide_audit.py"),
        "--png",
        str(output),
        "--layout-nodes",
        str(layout),
        "--output",
        str(audit),
        "--require-pass",
    ]
    subprocess.run(audit_cmd, check=True)
    print(output)
    print(manifest)
    print(csv_path)
    print(speaker)
    print(audit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

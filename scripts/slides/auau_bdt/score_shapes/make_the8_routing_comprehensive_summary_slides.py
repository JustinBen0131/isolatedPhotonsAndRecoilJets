#!/usr/bin/env python3
"""Render comprehensive corrected THE8 routed-BDT summary slide PNGs."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch, Rectangle

THIS_FILE = Path(__file__).resolve()
SCRIPTS_DIR = next((p for p in THIS_FILE.parents if p.name == "scripts"), THIS_FILE.parent)
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.append(str(SCRIPTS_DIR))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
BASE_OUT = REPO / "dataOutput/auauTightBDTValidation/THE59_model_comparison_remakes_20260620"
OUT_DIR = BASE_OUT / "routing_comprehensive_20260622"
PAYLOAD = OUT_DIR / "the8_routed_bdt_comprehensive_baseline_holdout_metrics_weighted.json"

TIMES_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_FONTS = [
    TIMES_DIR / "Times New Roman.ttf",
    TIMES_DIR / "Times New Roman Bold.ttf",
    TIMES_DIR / "Times New Roman Italic.ttf",
    TIMES_DIR / "Times New Roman Bold Italic.ttf",
]
FONT = "Times New Roman"

INK = "#111827"
MUTED = "#64748b"
GRID = "#d8e1ee"
EDGE = "#c8d3e2"
PANEL = "#f8fafc"
READOUT_FILL = "#fff7e6"
READOUT_EDGE = "#e6b85d"
GOOD = "#10713b"
BAD = "#b42318"
NEUTRAL = "#f8fafc"
ROUTE_LABELS = [
    ("cent3", "Cent3"),
    ("cent7", "Cent7"),
    ("et", "ET"),
    ("etcent3", "ET×C3"),
    ("etcent7", "ET×C7"),
]
PRIMARY_SURFACE_ROUTES = [("etcent3", "ET×C3"), ("etcent7", "ET×C7")]
METRICS = [
    ("auc", "ΔAUC ×10³", 1000.0, True, "{:+.1f}"),
    ("logloss", "Δlogloss ×10³", 1000.0, False, "{:+.1f}"),
    ("wp80_inclusive_fake", "Δfake rate pp", 100.0, False, "{:+.1f}"),
    ("median_gap", "Δmedian gap ×10²", 100.0, True, "{:+.1f}"),
]
CMAP = LinearSegmentedColormap.from_list("route_delta", ["#f4b6ae", "#fffdfa", "#bfe6c7"])


def configure_fonts() -> None:
    for font in TIMES_FONTS:
        if font.exists():
            font_manager.fontManager.addfont(str(font))
    plt.rcParams.update(
        {
            "font.family": FONT,
            "font.serif": [FONT, "Times", "DejaVu Serif"],
            "mathtext.fontset": "custom",
            "mathtext.rm": FONT,
            "mathtext.it": f"{FONT}:italic",
            "mathtext.bf": f"{FONT}:bold",
            "axes.unicode_minus": False,
        }
    )


def font_px(points: float) -> float:
    return points * SLIDE_DPI / 72.0


def bbox_px(fig: plt.Figure, artist) -> list[float]:
    fig.canvas.draw()
    bbox = artist.get_window_extent(renderer=fig.canvas.get_renderer())
    scale_x = fig.bbox.width / 2560.0
    scale_y = fig.bbox.height / 1440.0
    height = fig.bbox.height
    return [
        float(bbox.x0 / scale_x),
        float((height - bbox.y1) / scale_y),
        float(bbox.x1 / scale_x),
        float((height - bbox.y0) / scale_y),
    ]


def add_text(
    fig: plt.Figure,
    nodes: list[tuple[str, object, float, str, dict]],
    name: str,
    x: float,
    y: float,
    text: str,
    *,
    size: float,
    role: str = "audience",
    ha: str = "left",
    va: str = "center",
    weight: str = "normal",
    color: str = INK,
    linespacing: float = 1.0,
    zorder: int = 8,
    **flags,
) -> object:
    artist = fig.text(
        x,
        y,
        text,
        ha=ha,
        va=va,
        fontsize=size,
        family=FONT,
        fontweight=weight,
        color=color,
        linespacing=linespacing,
        zorder=zorder,
    )
    nodes.append((name, artist, size, role, flags))
    return artist


def rounded_box(fig: plt.Figure, x: float, y: float, w: float, h: float, *, face: str, edge: str, lw: float = 1.0, radius: float = 0.010, zorder: int = -10) -> None:
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            transform=fig.transFigure,
            boxstyle=f"round,pad=0.006,rounding_size={radius}",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=zorder,
        )
    )


def cell_lookup(payload: dict) -> dict[tuple[str, str, str], dict]:
    return {(c["group_kind"], c["group_label"], c["product"]): c for c in payload["cells"]}


def group_labels(payload: dict, kind: str) -> list[str]:
    labels: list[str] = []
    for cell in payload["cells"]:
        if cell["group_kind"] == kind and cell["product"] == "global" and cell["group_label"] not in labels:
            labels.append(cell["group_label"])
    return labels


def delta_value(cells: dict[tuple[str, str, str], dict], kind: str, label: str, route: str, metric: str) -> float:
    base = cells[(kind, label, "global")]["metrics"][metric]
    cur = cells[(kind, label, route)]["metrics"][metric]
    return float(cur - base)


def benefit_value(delta: float, good_positive: bool) -> float:
    return delta if good_positive else -delta


def metric_matrix(cells: dict[tuple[str, str, str], dict], kind: str, labels: list[str], routes: list[tuple[str, str]], metric: str, scale: float, good_positive: bool) -> tuple[np.ndarray, np.ndarray]:
    raw = np.full((len(labels), len(routes)), np.nan, dtype="float64")
    benefit = np.full_like(raw, np.nan)
    for r, row in enumerate(labels):
        for c, (route, _route_label) in enumerate(routes):
            d = delta_value(cells, kind, row, route, metric)
            raw[r, c] = d * scale
            benefit[r, c] = benefit_value(d, good_positive) * scale
    return raw, benefit


def draw_heatmap(
    ax: plt.Axes,
    values: np.ndarray,
    benefit: np.ndarray,
    row_labels: list[str],
    col_labels: list[str],
    *,
    title: str,
    fmt: str,
    show_ylabels: bool,
    text_size: float,
    title_size: float,
    separator_after: int | None = None,
) -> None:
    finite = np.abs(benefit[np.isfinite(benefit)])
    vmax = float(np.nanmax(finite)) if finite.size else 1.0
    vmax = max(vmax, 1.0)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    ax.imshow(benefit, aspect="auto", cmap=CMAP, norm=norm, interpolation="nearest")
    ax.set_title(title, loc="left", fontsize=title_size, fontweight="bold", pad=8)
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, fontsize=10.5, fontweight="bold")
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=5)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels if show_ylabels else ["" for _ in row_labels], fontsize=9.7 if len(row_labels) > 8 else 11.0, color=INK)
    ax.tick_params(axis="y", length=0, pad=4)
    ax.set_xticks(np.arange(-0.5, len(col_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_labels), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.6)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_edgecolor(EDGE)
        spine.set_linewidth(1.0)
    if separator_after is not None:
        ax.axhline(separator_after + 0.5, color="#334155", lw=1.7)
    for r in range(values.shape[0]):
        for c in range(values.shape[1]):
            val = values[r, c]
            if not np.isfinite(val):
                text = "n/a"
                color = MUTED
            else:
                text = fmt.format(val)
                color = INK
            ax.text(c, r, text, ha="center", va="center", fontsize=text_size, fontweight="bold", color=color)


def draw_metric_grid_slide(
    payload: dict,
    *,
    output: Path,
    title: str,
    kinds_and_rows: list[tuple[str, str]],
    readout_lines: list[str],
) -> Path:
    configure_fonts()
    cells = cell_lookup(payload)
    labels: list[str] = []
    row_kinds: list[str] = []
    separators: list[int] = []
    for idx, (kind, heading) in enumerate(kinds_and_rows):
        rows = group_labels(payload, kind)
        if idx > 0:
            separators.append(len(labels) - 1)
        for row in rows:
            labels.append(f"{heading} {row}" if len(kinds_and_rows) > 1 else row)
            row_kinds.append(kind)
    routes = ROUTE_LABELS
    col_labels = [label for _key, label in routes]

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[tuple[str, object, float, str, dict]] = []
    add_text(fig, nodes, "title", 0.040, 0.938, title, size=40.0, role="title", weight="bold", title_anchor=True)
    add_text(
        fig,
        nodes,
        "contract line",
        0.043,
        0.865,
        "Same fixed holdout. Values are deltas versus the global 14-feature baseline; green means better.",
        size=14.3,
        color=MUTED,
    )

    positions = [
        [0.060, 0.515, 0.420, 0.270],
        [0.535, 0.515, 0.420, 0.270],
        [0.060, 0.185, 0.420, 0.270],
        [0.535, 0.185, 0.420, 0.270],
    ]
    for i, ((metric, metric_title, scale, good_positive, fmt), pos) in enumerate(zip(METRICS, positions)):
        raw_rows = []
        benefit_rows = []
        for kind, row_label in zip(row_kinds, labels):
            clean_label = row_label.split(" ", 1)[1] if len(kinds_and_rows) > 1 else row_label
            raw, benefit = metric_matrix(cells, kind, [clean_label], routes, metric, scale, good_positive)
            raw_rows.append(raw[0])
            benefit_rows.append(benefit[0])
        values = np.asarray(raw_rows)
        benefit = np.asarray(benefit_rows)
        ax = fig.add_axes(pos)
        draw_heatmap(
            ax,
            values,
            benefit,
            labels,
            col_labels,
            title=metric_title,
            fmt=fmt,
            show_ylabels=i in (0, 2),
            text_size=8.8 if len(labels) > 8 else 10.5,
            title_size=14.0,
            separator_after=separators[0] if separators else None,
        )

    rounded_box(fig, 0.050, 0.055, 0.910, 0.085, face=READOUT_FILL, edge=READOUT_EDGE, lw=1.0)
    add_text(fig, nodes, "readout title", 0.066, 0.101, "Readout", size=17.8, weight="bold")
    add_text(fig, nodes, "readout body 1", 0.158, 0.112, readout_lines[0], size=12.8, color=INK, colon_style_exception=True)
    add_text(fig, nodes, "readout body 2", 0.158, 0.083, readout_lines[1], size=12.8, color=INK, colon_style_exception=True)

    return save_slide(fig, nodes, output, payload, extra={"kinds": [k for k, _h in kinds_and_rows], "routes": dict(routes), "readout": readout_lines})


def save_slide(fig: plt.Figure, nodes: list[tuple[str, object, float, str, dict]], output: Path, payload: dict, *, extra: dict) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=SLIDE_DPI)
    layout_path = output.with_suffix(".layout_nodes.json")
    layout_nodes = []
    for name, artist, points, role, flags in nodes:
        layout_nodes.append(
            {
                "name": name,
                "kind": "text",
                "role": role,
                "text": artist.get_text(),
                "font_px": font_px(points),
                "bbox": bbox_px(fig, artist),
                "title_anchor": role == "title",
                **flags,
            }
        )
    layout_path.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "title_axis_x": 0.040 * 2560,
                "minimum_audience_font_px": font_px(10.5),
                "minimum_title_font_px": font_px(35.0),
                "minimum_plot_annotation_font_px": font_px(8.5),
                "nodes": layout_nodes,
            },
            indent=2,
        )
        + "\n"
    )
    manifest_path = output.with_suffix(".manifest.json")
    manifest = {
        "schema": "the8_routing_comprehensive_summary_slide_v1",
        "png": str(output),
        "payload": str(PAYLOAD),
        "layout_nodes": str(layout_path),
        "source": payload.get("source"),
        "selection": payload.get("selection"),
        "labels": payload.get("labels"),
        "metric_weighting": payload.get("metric_weighting"),
        **extra,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    plt.close(fig)
    return output


def write_datapoints(payload: dict, path: Path) -> None:
    cells = cell_lookup(payload)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["group_kind", "group_label", "route", "metric", "global", "route_value", "delta", "scaled_delta"])
        for cell in payload["cells"]:
            if cell["product"] == "global":
                continue
            base = cells[(cell["group_kind"], cell["group_label"], "global")]["metrics"]
            cur = cell["metrics"]
            for metric, _title, scale, _good, _fmt in METRICS:
                writer.writerow(
                    [
                        cell["group_kind"],
                        cell["group_label"],
                        cell["product"],
                        metric,
                        base.get(metric),
                        cur.get(metric),
                        cur.get(metric) - base.get(metric),
                        (cur.get(metric) - base.get(metric)) * scale,
                    ]
                )


def summarize_range(payload: dict, kind: str, routes: list[str], metric: str, *, good_positive: bool, scale: float) -> tuple[float, float, int, int]:
    cells = cell_lookup(payload)
    labels = group_labels(payload, kind)
    benefits = []
    for row in labels:
        for route in routes:
            d = delta_value(cells, kind, row, route, metric)
            benefits.append(benefit_value(d, good_positive) * scale)
    arr = np.asarray(benefits, dtype="float64")
    return float(np.nanmin(arr)), float(np.nanmax(arr)), int(np.sum(arr > 0)), int(np.sum(np.isfinite(arr)))


def draw_surface_slide(payload: dict, output: Path) -> Path:
    configure_fonts()
    cells = cell_lookup(payload)
    et_labels = group_labels(payload, "et")
    fine_labels = [label.replace("%", "") for label in group_labels(payload, "fine_cent")]
    coarse_labels = [label.replace("%", "") for label in group_labels(payload, "coarse_cent")]
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[tuple[str, object, float, str, dict]] = []
    add_text(fig, nodes, "title", 0.040, 0.938, "ET × centrality routing surface", size=40.0, role="title", weight="bold", title_anchor=True)
    add_text(
        fig,
        nodes,
        "contract line",
        0.043,
        0.865,
        "Same fixed holdout. Values are deltas versus the global 14-feature baseline; green means better.",
        size=14.3,
        color=MUTED,
    )

    panels = [
        ("et_coarse_cent", "etcent3", "ET×C3 ΔAUC ×10³", "auc", 1000.0, True, coarse_labels, [0.055, 0.535, 0.405, 0.245]),
        ("et_coarse_cent", "etcent3", "ET×C3 Δfake rate pp", "wp80_inclusive_fake", 100.0, False, coarse_labels, [0.055, 0.230, 0.405, 0.245]),
        ("et_fine_cent", "etcent7", "ET×C7 ΔAUC ×10³", "auc", 1000.0, True, fine_labels, [0.540, 0.535, 0.405, 0.245]),
        ("et_fine_cent", "etcent7", "ET×C7 Δfake rate pp", "wp80_inclusive_fake", 100.0, False, fine_labels, [0.540, 0.230, 0.405, 0.245]),
    ]
    for kind, route, title, metric, scale, good_positive, xlabels, pos in panels:
        ylabels = [label.replace(" GeV", "") for label in et_labels]
        values = np.full((len(ylabels), len(xlabels)), np.nan, dtype="float64")
        benefit = np.full_like(values, np.nan)
        for i, et_label in enumerate(et_labels):
            for j, xlab in enumerate(xlabels):
                cent_label = xlab + "%"
                group_label = f"{et_label} × {cent_label}"
                if (kind, group_label, route) not in cells:
                    continue
                d = delta_value(cells, kind, group_label, route, metric)
                values[i, j] = d * scale
                benefit[i, j] = benefit_value(d, good_positive) * scale
        ax = fig.add_axes(pos)
        draw_heatmap(
            ax,
            values,
            benefit,
            ylabels,
            xlabels,
            title=title,
            fmt="{:+.1f}",
            show_ylabels=pos[0] < 0.1,
            text_size=8.4 if len(xlabels) <= 3 else 7.0,
            title_size=14.0,
        )

    rounded_box(fig, 0.050, 0.055, 0.910, 0.095, face=READOUT_FILL, edge=READOUT_EDGE, lw=1.0)
    all_auc = summarize_range(payload, "et_fine_cent", ["etcent7"], "auc", good_positive=True, scale=1000.0)
    all_fake = summarize_range(payload, "et_fine_cent", ["etcent7"], "wp80_inclusive_fake", good_positive=False, scale=100.0)
    add_text(fig, nodes, "readout title", 0.066, 0.107, "Readout", size=17.8, weight="bold")
    add_text(
        fig,
        nodes,
        "readout body 1",
        0.158,
        0.119,
        f"ET×fine is the broadest route: AUC benefit is positive in {all_auc[2]}/{all_auc[3]} ET×fine cells, spanning {all_auc[0]:+.1f} to {all_auc[1]:+.1f} ×10³.",
        size=12.5,
        colon_style_exception=True,
    )
    add_text(
        fig,
        nodes,
        "readout body 2",
        0.158,
        0.087,
        f"Fake-rate benefit is positive in {all_fake[2]}/{all_fake[3]} ET×fine cells; this surface shows where routing reduces inclusive-MC leakage.",
        size=12.5,
        colon_style_exception=True,
    )
    return save_slide(fig, nodes, output, payload, extra={"surface_routes": dict(PRIMARY_SURFACE_ROUTES)})


def render_all() -> list[Path]:
    payload = json.loads(PAYLOAD.read_text())
    write_datapoints(payload, OUT_DIR / "the8_routing_comprehensive_metric_deltas.csv")
    centrality_auc = summarize_range(payload, "fine_cent", ["cent3", "cent7", "et", "etcent3", "etcent7"], "auc", good_positive=True, scale=1000.0)
    centrality_fake = summarize_range(payload, "fine_cent", ["cent3", "cent7", "et", "etcent3", "etcent7"], "wp80_inclusive_fake", good_positive=False, scale=100.0)
    et_auc = summarize_range(payload, "et", ["cent3", "cent7", "et", "etcent3", "etcent7"], "auc", good_positive=True, scale=1000.0)
    et_fake = summarize_range(payload, "et", ["cent3", "cent7", "et", "etcent3", "etcent7"], "wp80_inclusive_fake", good_positive=False, scale=100.0)
    outputs = [
        draw_metric_grid_slide(
            payload,
            output=OUT_DIR / "the8_routing_comprehensive_centrality_metrics_v1.png",
            title="Routing summary versus centrality",
            kinds_and_rows=[("coarse_cent", "C3"), ("fine_cent", "C7")],
            readout_lines=[
                f"Fine-centrality view: AUC benefit is positive in {centrality_auc[2]}/{centrality_auc[3]} cells; range {centrality_auc[0]:+.1f} to {centrality_auc[1]:+.1f} ×10³.",
                f"Fine-centrality fake-rate benefit is favorable in {centrality_fake[2]}/{centrality_fake[3]} cells versus global.",
            ],
        ),
        draw_metric_grid_slide(
            payload,
            output=OUT_DIR / "the8_routing_comprehensive_et_metrics_v1.png",
            title="Routing summary versus cluster ET",
            kinds_and_rows=[("et", "")],
            readout_lines=[
                f"ET view: AUC benefit is positive in {et_auc[2]}/{et_auc[3]} cells; ET×centrality gives the largest coherent gains.",
                f"Fake-rate benefit is positive in {et_fake[2]}/{et_fake[3]} ET cells; negative numbers in the table mean fewer inclusive-MC fakes.",
            ],
        ),
        draw_surface_slide(payload, OUT_DIR / "the8_routing_comprehensive_et_cent_surface_v1.png"),
    ]
    script = OUT_DIR / "the8_routing_comprehensive_summary_speaker_script.md"
    script.write_text(
        "\n".join(
            [
                "# Corrected THE8 routed-BDT summary",
                "",
                "All panels compare a routed model against the global 14-feature baseline on the same reconstructed baseline holdout.",
                "The class definition is Signal MC from embeddedPhoton rows with is_signal==1 versus inclusive embeddedJet MC with no truth-background filter.",
                "Green always means better than global: higher AUC, lower logloss, lower fake rate, or larger signal-inclusive median gap.",
                "",
                "The centrality slide shows whether centrality-only and ET-aware routing help uniformly in coarse and fine centrality bins.",
                "The ET slide asks the same question in cluster-ET bins.",
                "The surface slide keeps the ET and centrality axes visible at the same time for ET×coarse and ET×fine routing.",
                "",
            ]
        )
    )
    return outputs


def main() -> None:
    for output in render_all():
        print(output)
        print(output.with_suffix(".manifest.json"))
        print(output.with_suffix(".layout_nodes.json"))


if __name__ == "__main__":
    main()

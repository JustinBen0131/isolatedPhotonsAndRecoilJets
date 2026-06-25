#!/usr/bin/env python3
"""Render delta-bar summary slides for corrected THE8 routing diagnostics."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, Rectangle

THIS_FILE = Path(__file__).resolve()
SCRIPTS_DIR = next((p for p in THIS_FILE.parents if p.name == "scripts"), THIS_FILE.parent)
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.append(str(SCRIPTS_DIR))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
BASE_OUT = REPO / "dataOutput/auauTightBDTValidation/THE59_model_comparison_remakes_20260620"
OUT_DIR = BASE_OUT / "routing_metric_summary_20260622"
COARSE_CENT_PAYLOAD = BASE_OUT / "centrality_routing_20260622/the8_centrality_routing_score_payload_baseline_holdout.json"
FINE_CENT_PAYLOAD = BASE_OUT / "centrality_routing_20260622/the8_centrality_routing_fine_rows_payload_baseline_holdout.json"
ET_PAYLOAD = BASE_OUT / "et_routing_20260622/the8_et_routing_score_payload_baseline_holdout.json"
DEFAULT_OUTPUT = OUT_DIR / "the8_routing_delta_bars_baseline_holdout_v1.png"

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
GRID = "#dfe7f1"
EDGE = "#c8d3e2"
PANEL = "#f8fafc"
GOOD = "#eaf7ed"
COARSE = "#0f9f8e"
FINE = "#8a62c4"
ET = "#d97706"
GLOBAL = "#5265d6"
READOUT_FILL = "#fff7e6"
READOUT_EDGE = "#e7b85d"

ROUTE_STYLE = {
    "coarse": ("3 cent BDTs", COARSE),
    "fine": ("7 cent BDTs", FINE),
    "et": ("8 ET BDTs", ET),
}


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


def rounded_box(
    fig: plt.Figure,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    face: str,
    edge: str,
    lw: float = 1.0,
    radius: float = 0.010,
    zorder: int = -10,
) -> None:
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


def cell_map(payload: dict) -> dict[tuple[str, str], dict]:
    return {(cell["row"], cell["col"]): cell for cell in payload["cells"]}


def row_order(payload: dict) -> list[str]:
    if payload.get("row_order"):
        return list(payload["row_order"])
    rows: list[str] = []
    for cell in payload["cells"]:
        row = str(cell["row"])
        if row not in rows:
            rows.append(row)
    return rows


def deltas(payload: dict, routes: list[str]) -> dict[str, dict[str, list[float]]]:
    cells = cell_map(payload)
    rows = row_order(payload)
    out: dict[str, dict[str, list[float]]] = {}
    for route in routes:
        metrics = {"auc": [], "fake_pp": [], "logloss": [], "median_gap": []}
        for row in rows:
            base = cells[(row, "global")]["metrics"]
            cur = cells[(row, route)]["metrics"]
            metrics["auc"].append(float(cur["auc"] - base["auc"]))
            metrics["fake_pp"].append(float(100.0 * (cur["wp80_inclusive_fake"] - base["wp80_inclusive_fake"])))
            metrics["logloss"].append(float(cur["logloss"] - base["logloss"]))
            metrics["median_gap"].append(float(cur["median_gap"] - base["median_gap"]))
        out[route] = metrics
    return out


def axis_limits(values: list[float], *, good_positive: bool) -> tuple[float, float]:
    lo, hi = min(values), max(values)
    span = max(hi - lo, 1e-6)
    pad = span * 0.22
    if good_positive:
        lo = min(lo - pad, -0.00025 if hi < 0.004 else lo - pad)
        hi = hi + pad
    else:
        lo = lo - pad
        hi = max(hi + pad, 0.08)
    return lo, hi


def annotate_bar(ax: plt.Axes, x: float, y: float, label: str, *, xlim: tuple[float, float], size: float) -> None:
    span = xlim[1] - xlim[0]
    dx = span * 0.012
    if x >= 0:
        xpos = min(x + dx, xlim[1] - span * 0.02)
        ha = "left"
    else:
        xpos = max(x - dx, xlim[0] + span * 0.02)
        ha = "right"
    ax.text(xpos, y, label, va="center", ha=ha, fontsize=size, color=INK, fontweight="bold")


def draw_delta_axis(
    fig: plt.Figure,
    payload: dict,
    routes: list[str],
    *,
    ax_pos: list[float],
    metric: str,
    title: str,
    good_positive: bool,
    show_ylabels: bool,
) -> plt.Axes:
    ax = fig.add_axes(ax_pos, zorder=4)
    rows = row_order(payload)
    delta = deltas(payload, routes)
    y_centers = list(range(len(rows)))
    bar_h = 0.64 / max(1, len(routes))
    all_values = [value for route in routes for value in delta[route][metric]]
    xlim = axis_limits(all_values, good_positive=good_positive)
    good_lo, good_hi = (0, xlim[1]) if good_positive else (xlim[0], 0)
    ax.axvspan(good_lo, good_hi, color=GOOD, zorder=0)
    ax.axvline(0, color="#334155", lw=1.1, zorder=2)
    for ridx, route in enumerate(routes):
        label, color = ROUTE_STYLE[route]
        offset = (ridx - (len(routes) - 1) / 2.0) * bar_h
        ys = [y + offset for y in y_centers]
        vals = delta[route][metric]
        ax.barh(ys, vals, height=bar_h * 0.86, color=color, edgecolor=color, alpha=0.96, label=label, zorder=3)
        for y, val in zip(ys, vals):
            if metric == "auc":
                text = f"{val:+.3f}"
            elif metric == "fake_pp":
                text = f"{val:+.2f}"
            elif metric == "logloss":
                text = f"{val:+.3f}"
            else:
                text = f"{val:+.2f}"
            annotate_bar(ax, val, y, text, xlim=xlim, size=7.9 if len(rows) > 3 else 8.5)
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.75, len(rows) - 0.25)
    ax.invert_yaxis()
    ax.set_title(title, fontsize=12.0, fontweight="bold", loc="left", pad=4)
    ax.set_facecolor("white")
    ax.grid(True, axis="x", color=GRID, linewidth=0.9)
    ax.set_yticks(y_centers)
    ax.set_yticklabels(rows if show_ylabels else ["" for _ in rows], fontsize=8.8 if len(rows) > 3 else 10.0, color=MUTED)
    ax.tick_params(axis="x", labelsize=8.8, colors=MUTED, pad=2)
    ax.tick_params(axis="y", length=0, pad=2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#94a3b8")
    if metric == "fake_pp":
        ax.set_xlabel("percentage points; negative is better", fontsize=8.8, color=MUTED, labelpad=1)
    elif metric == "auc":
        ax.set_xlabel("positive is better", fontsize=8.8, color=MUTED, labelpad=1)
    return ax


def write_datapoints(payloads: list[tuple[str, dict, list[str]]], output_csv: Path) -> None:
    rows_out: list[dict[str, object]] = []
    for axis, payload, routes in payloads:
        delta = deltas(payload, routes)
        rows = row_order(payload)
        for bin_label in rows:
            idx = rows.index(bin_label)
            for route in routes:
                rows_out.append(
                    {
                        "axis": axis,
                        "bin": bin_label,
                        "route": route,
                        "route_label": ROUTE_STYLE[route][0],
                        "delta_auc": delta[route]["auc"][idx],
                        "delta_fake_rate_pp": delta[route]["fake_pp"][idx],
                        "delta_logloss": delta[route]["logloss"][idx],
                        "delta_median_gap": delta[route]["median_gap"][idx],
                    }
                )
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows_out[0].keys()))
        writer.writeheader()
        writer.writerows(rows_out)


def delta_ranges(payload: dict, routes: list[str]) -> tuple[tuple[float, float], tuple[float, float]]:
    delta = deltas(payload, routes)
    auc = [v for route in routes for v in delta[route]["auc"]]
    fake = [v for route in routes for v in delta[route]["fake_pp"]]
    return (min(auc), max(auc)), (min(fake), max(fake))


def readout(coarse: dict, fine: dict, et_payload: dict) -> tuple[str, str]:
    fine_auc, fine_fake = delta_ranges(fine, ["coarse", "fine"])
    et_delta = deltas(et_payload, ["et"])["et"]
    wins = sum(1 for auc, fake, logloss in zip(et_delta["auc"], et_delta["fake_pp"], et_delta["logloss"]) if auc > 0 and fake < 0 and logloss < 0)
    return (
        f"Centrality routing is subtle. Fine-bin ΔAUC {fine_auc[0]:+.3f} to {fine_auc[1]:+.3f}; Δfake {fine_fake[0]:+.2f} to {fine_fake[1]:+.2f} pp.",
        f"ET routing is more coherent. ΔAUC positive and Δfake/Δlogloss lower in {wins}/8 ET bins.",
    )


def render(output_path: Path = DEFAULT_OUTPUT) -> Path:
    configure_fonts()
    coarse = json.loads(COARSE_CENT_PAYLOAD.read_text())
    fine = json.loads(FINE_CENT_PAYLOAD.read_text())
    et_payload = json.loads(ET_PAYLOAD.read_text())
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[tuple[str, object, float, str, dict]] = []

    add_text(
        fig,
        nodes,
        "title",
        0.040,
        0.955,
        "Routing deltas versus the global 14-feature BDT",
        size=34.5,
        role="title",
        va="top",
        weight="bold",
        title_anchor=True,
    )
    rounded_box(fig, 0.040, 0.836, 0.920, 0.052, face="white", edge=EDGE, lw=0.95, radius=0.006)
    add_text(fig, nodes, "holdout contract", 0.058, 0.864, "Same fixed holdout in every bar = 10% row holdout, 15-35 GeV, 0-80%", size=12.2, weight="bold", title_band_exception=True)
    add_text(fig, nodes, "sample contract", 0.548, 0.864, "Signal Photon12+20; inclusive Jet12+20+30+40, no truth filter", size=11.4, color=MUTED, title_band_exception=True)

    plot_x0 = 0.183
    auc_x = plot_x0
    fake_x = 0.580
    ax_w = 0.364
    row_h = 0.178
    row_ys = [0.618, 0.382, 0.160]
    row_specs = [
        ("Coarse centrality bins", "0-20, 20-50, 50-80%", coarse, ["coarse", "fine"], FINE),
        ("Fine centrality bins", "0-10 through 60-80%", fine, ["coarse", "fine"], FINE),
        ("ET bins", "15-35 GeV", et_payload, ["et"], ET),
    ]
    for idx, (title, subtitle, payload, routes, accent) in enumerate(row_specs):
        y = row_ys[idx]
        rounded_box(fig, 0.040, y - 0.010, 0.920, row_h + 0.026, face=PANEL, edge=EDGE, lw=0.9, radius=0.010)
        fig.patches.append(Rectangle((0.054, y + row_h - 0.025), 0.008, 0.044, transform=fig.transFigure, facecolor=accent, edgecolor=accent, zorder=2))
        add_text(fig, nodes, f"{title} row title", 0.070, y + row_h - 0.004, title, size=14.5, weight="bold")
        add_text(fig, nodes, f"{title} row subtitle", 0.070, y + row_h - 0.029, subtitle, size=10.7, color=MUTED)
        draw_delta_axis(
            fig,
            payload,
            routes,
            ax_pos=[auc_x, y, ax_w, row_h],
            metric="auc",
            title="ΔAUC",
            good_positive=True,
            show_ylabels=True,
        )
        ax_fake = draw_delta_axis(
            fig,
            payload,
            routes,
            ax_pos=[fake_x, y, ax_w, row_h],
            metric="fake_pp",
            title="Δ fake rate",
            good_positive=False,
            show_ylabels=False,
        )
        if idx == 0:
            leg = ax_fake.legend(
                loc="lower right",
                bbox_to_anchor=(1.0, 1.08),
                ncol=3,
                frameon=True,
                facecolor="white",
                edgecolor=EDGE,
                fontsize=9.0,
                handlelength=1.5,
                columnspacing=1.0,
            )
            leg.get_frame().set_alpha(0.96)

    read1, read2 = readout(coarse, fine, et_payload)
    rounded_box(fig, 0.040, 0.045, 0.920, 0.084, face=READOUT_FILL, edge=READOUT_EDGE, lw=1.05, radius=0.008)
    add_text(fig, nodes, "readout title", 0.058, 0.091, "Readout", size=16.6, weight="bold")
    add_text(fig, nodes, "readout one", 0.142, 0.102, read1, size=11.7, color=INK)
    add_text(fig, nodes, "readout two", 0.142, 0.075, read2, size=11.7, color=INK)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    datapoints_path = output_path.with_name(output_path.stem + "_datapoints.csv")
    write_datapoints(
        [
            ("coarse_centrality", coarse, ["coarse", "fine"]),
            ("fine_centrality", fine, ["coarse", "fine"]),
            ("et", et_payload, ["et"]),
        ],
        datapoints_path,
    )
    fig.savefig(output_path, dpi=SLIDE_DPI)

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
    layout_path = output_path.with_suffix(".layout_nodes.json")
    layout_path.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "title_axis_x": 0.040 * 2560,
                "minimum_audience_font_px": font_px(10.2),
                "minimum_title_font_px": font_px(30.0),
                "minimum_plot_annotation_font_px": font_px(7.8),
                "nodes": layout_nodes,
            },
            indent=2,
        )
        + "\n"
    )
    manifest_path = output_path.with_suffix(".manifest.json")
    manifest_path.write_text(
        json.dumps(
            {
                "schema": "the8_routing_delta_bar_slide_v1",
                "png": str(output_path),
                "coarse_centrality_payload": str(COARSE_CENT_PAYLOAD),
                "fine_centrality_payload": str(FINE_CENT_PAYLOAD),
                "et_payload": str(ET_PAYLOAD),
                "datapoints_csv": str(datapoints_path),
                "layout_nodes": str(layout_path),
                "selection": fine.get("selection"),
                "labels": fine.get("labels"),
                "coarse_model_columns": coarse.get("model_columns"),
                "fine_model_columns": fine.get("model_columns"),
                "et_model_columns": et_payload.get("model_columns"),
                "readout": list(readout(coarse, fine, et_payload)),
                "scope_note": "Bars show routed-model metric deltas relative to centAsFeatBase3x3_pt15to35 on the same corrected baseline holdout. ET x centrality routed variants are not included until scored on this same holdout contract.",
            },
            indent=2,
        )
        + "\n"
    )
    script_path = output_path.with_name(output_path.stem + "_speaker_script.md")
    script_path.write_text(
        "\n".join(
            [
                "# Routing deltas versus the global 14-feature BDT",
                "",
                "This is the same corrected THE8 holdout contract as the score-shape tables, but shown as deltas against the global 14-feature BDT.",
                "",
                "Positive AUC bars are good. Negative fake-rate bars are good. This makes the small routing effects visible without making the audience compare nearly overlapping absolute curves.",
                "",
                read1,
                "",
                read2,
                "",
                "The plotted input is source-class signal Photon12+20 and source-class inclusive Jet12+20+30+40 with no truth-background filter.",
                "",
            ]
        )
    )
    plt.close(fig)
    return output_path


def main() -> None:
    out = render()
    print(out)
    print(out.with_suffix(".manifest.json"))
    print(out.with_suffix(".layout_nodes.json"))
    print(out.with_name(out.stem + "_datapoints.csv"))
    print(out.with_name(out.stem + "_speaker_script.md"))


if __name__ == "__main__":
    main()

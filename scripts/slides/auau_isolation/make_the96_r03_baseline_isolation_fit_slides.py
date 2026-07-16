#!/usr/bin/env python3
"""Build slide-ready R=0.3 or R=0.4 AuAu isolation-fit candidates."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.ticker import MaxNLocator


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))

from slides.common.slide_defaults import (  # noqa: E402
    SLIDE_DPI,
    SLIDE_HEIGHT_PX,
    SLIDE_WIDTH_PX,
    slide_figsize,
)


CONE = "R30"
RADIUS_TEXT = "0.3"
RADIUS_TAG = "r03"
ALGORITHM = "baseline"
PT_SCHEME = "nominal"
PT_RANGE_TEXT = "15-35"
PT_BIN_COUNT = 6
CAMPAIGN_LABEL = "THE96"


def configure_paths(
    cone: str,
    algorithm: str = "baseline",
    pt_scheme: str = "nominal",
    *,
    source_dir_override: Path | None = None,
    output_dir_override: Path | None = None,
    source_prefix_root: str = "the96",
    campaign_label: str = "THE96",
) -> None:
    global CONE, RADIUS_TEXT, RADIUS_TAG, ALGORITHM
    global PT_SCHEME, PT_RANGE_TEXT, PT_BIN_COUNT, CAMPAIGN_LABEL
    global SOURCE_DIR, POINTS_CSV, FLAT_CSV, LINEAR_CSV, SOURCE_MANIFEST
    global OUT_DIR, GRID_PNG, GRID_SCRIPT, GRID_LAYOUT, GRID_AUDIT
    global CENT_PNG, CENT_SCRIPT, CENT_LAYOUT, CENT_AUDIT, OUT_MANIFEST

    if cone not in {"R30", "R40"}:
        raise ValueError(f"unsupported cone: {cone}")
    if algorithm not in {"baseline", "topocluster"}:
        raise ValueError(f"unsupported algorithm: {algorithm}")
    if pt_scheme not in {"nominal", "expanded"}:
        raise ValueError(f"unsupported pT scheme: {pt_scheme}")
    CONE = cone
    ALGORITHM = algorithm
    PT_SCHEME = pt_scheme
    PT_RANGE_TEXT = "15-35" if pt_scheme == "nominal" else "5-40"
    PT_BIN_COUNT = 6 if pt_scheme == "nominal" else 13
    CAMPAIGN_LABEL = campaign_label.strip().upper().replace("-", "")
    if not CAMPAIGN_LABEL or not CAMPAIGN_LABEL.isalnum():
        raise ValueError(f"invalid campaign label: {campaign_label!r}")
    source_prefix_root = source_prefix_root.strip()
    if not source_prefix_root or any(
        character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-"
        for character in source_prefix_root
    ):
        raise ValueError(f"invalid source prefix root: {source_prefix_root!r}")
    RADIUS_TEXT = "0.3" if cone == "R30" else "0.4"
    RADIUS_TAG = "r03" if cone == "R30" else "r04"
    histogram_tag = "r30" if cone == "R30" else "r40"

    SOURCE_DIR = (
        source_dir_override.expanduser().resolve()
        if source_dir_override is not None
        else REPO
        / "dataOutput/auau/isolation_baseline5pct_the96_20260709"
        / f"{RADIUS_TAG}_{ALGORITHM}_{PT_SCHEME}_{'15to35' if PT_SCHEME == 'nominal' else '5to40'}"
    )
    source_prefix = f"{source_prefix_root}_auau_iso_{histogram_tag}_{ALGORITHM}_{PT_SCHEME}"
    POINTS_CSV = SOURCE_DIR / f"{source_prefix}_cutoff_points.csv"
    FLAT_CSV = SOURCE_DIR / f"{source_prefix}_flat_fits.csv"
    LINEAR_CSV = SOURCE_DIR / f"{source_prefix}_linear_fits.csv"
    SOURCE_MANIFEST = SOURCE_DIR / f"{source_prefix}_manifest.json"

    OUT_DIR = (
        output_dir_override.expanduser().resolve()
        if output_dir_override is not None
        else SOURCE_DIR / "slide_candidates"
    )
    slide_prefix = f"auau_{RADIUS_TAG}_{ALGORITHM}_isolation"
    GRID_PNG = OUT_DIR / f"{slide_prefix}_pt_flat_fits_slide.png"
    GRID_SCRIPT = OUT_DIR / f"{slide_prefix}_pt_flat_fits_speaker_script.md"
    GRID_LAYOUT = OUT_DIR / f"{slide_prefix}_pt_flat_fits_layout_nodes.json"
    GRID_AUDIT = OUT_DIR / f"{slide_prefix}_pt_flat_fits_slide_audit.json"
    CENT_PNG = OUT_DIR / f"{slide_prefix}_centrality_fit_test_slide.png"
    CENT_SCRIPT = OUT_DIR / f"{slide_prefix}_centrality_fit_test_speaker_script.md"
    CENT_LAYOUT = OUT_DIR / f"{slide_prefix}_centrality_fit_test_layout_nodes.json"
    CENT_AUDIT = OUT_DIR / f"{slide_prefix}_centrality_fit_test_slide_audit.json"
    OUT_MANIFEST = OUT_DIR / f"{slide_prefix}_fit_slides_manifest.json"


configure_paths("R30")

TIMES_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_FONTS = [
    TIMES_DIR / "Times New Roman.ttf",
    TIMES_DIR / "Times New Roman Bold.ttf",
    TIMES_DIR / "Times New Roman Italic.ttf",
    TIMES_DIR / "Times New Roman Bold Italic.ttf",
]

INK = "#111827"
MUTED = "#475569"
GRID = "#D7DEE8"
EDGE = "#B7C3D2"
BLUE = "#0072B2"
GREEN = "#009E73"
ORANGE = "#D55E00"
RED = "#B91C1C"
NOTE_FILL = "#FFF7ED"
NOTE_EDGE = "#F59E0B"
COLORS = {0.7: BLUE, 0.8: GREEN, 0.9: ORANGE}
MARKERS = {0.7: "o", 0.8: "s", 0.9: "^"}
LABELS = {0.7: "70%", 0.8: "80%", 0.9: "90%"}


def setup_style() -> None:
    for font in TIMES_FONTS:
        if font.exists():
            font_manager.fontManager.addfont(str(font))
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "custom",
            "mathtext.rm": "Times New Roman",
            "mathtext.it": "Times New Roman:italic",
            "mathtext.bf": "Times New Roman:bold",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": INK,
            "axes.linewidth": 1.0,
            "xtick.color": INK,
            "ytick.color": INK,
            "axes.unicode_minus": False,
        }
    )


def font_px(points: float) -> float:
    return points * SLIDE_DPI / 72.0


def text_bbox_px(fig: plt.Figure, artist) -> list[float]:
    fig.canvas.draw()
    bbox = artist.get_window_extent(renderer=fig.canvas.get_renderer())
    scale_x = fig.bbox.width / SLIDE_WIDTH_PX
    scale_y = fig.bbox.height / SLIDE_HEIGHT_PX
    height = fig.bbox.height
    return [
        float(bbox.x0 / scale_x),
        float((height - bbox.y1) / scale_y),
        float(bbox.x1 / scale_x),
        float((height - bbox.y0) / scale_y),
    ]


def add_text(
    fig: plt.Figure,
    nodes: list[dict],
    name: str,
    x: float,
    y: float,
    text: str,
    *,
    size: float,
    role: str = "audience",
    weight: str = "normal",
    color: str = INK,
    ha: str = "left",
    va: str = "center",
    linespacing: float = 1.0,
    rotation: float = 0.0,
    **flags,
) -> None:
    artist = fig.text(
        x,
        y,
        text,
        fontsize=size,
        fontweight=weight,
        color=color,
        ha=ha,
        va=va,
        linespacing=linespacing,
        rotation=rotation,
        zorder=10,
    )
    nodes.append(
        {
            "name": name,
            "kind": "text",
            "role": role,
            "text": text,
            "font_px": font_px(size),
            "bbox": text_bbox_px(fig, artist),
            "title_anchor": role == "title",
            **flags,
        }
    )


def add_artist_node(
    fig: plt.Figure,
    nodes: list[dict],
    name: str,
    artist,
    *,
    size: float,
    text: str,
) -> None:
    nodes.append(
        {
            "name": name,
            "kind": "text",
            "role": "plot_annotation",
            "text": text,
            "font_px": font_px(size),
            "bbox": text_bbox_px(fig, artist),
            "title_anchor": False,
        }
    )


def require_vertical_symmetry(
    nodes: list[dict],
    *,
    upper_name: str,
    middle_name: str,
    lower_name: str,
    tolerance_px: float = 3.0,
) -> None:
    by_name = {node.get("name"): node for node in nodes}
    upper = by_name[upper_name]["bbox"]
    middle = by_name[middle_name]["bbox"]
    lower = by_name[lower_name]["bbox"]
    upper_gap = float(middle[1]) - float(upper[3])
    lower_gap = float(lower[1]) - float(middle[3])
    delta = abs(upper_gap - lower_gap)
    if upper_gap < 0 or lower_gap < 0 or delta > tolerance_px:
        raise RuntimeError(
            f"vertical symmetry failed for {middle_name}: "
            f"upper_gap={upper_gap:.2f}px lower_gap={lower_gap:.2f}px delta={delta:.2f}px"
        )


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
) -> None:
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            transform=fig.transFigure,
            boxstyle="round,pad=0.006,rounding_size=0.006",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=4,
        )
    )


def load_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    points = pd.read_csv(POINTS_CSV)
    flat = pd.read_csv(FLAT_CSV)
    linear = pd.read_csv(LINEAR_CSV)
    expected_points = PT_BIN_COUNT * 16 * 3
    if len(points) != expected_points or len(flat) != 48 or len(linear) != 3:
        raise RuntimeError(
            f"unexpected table sizes: points={len(points)}/{expected_points} "
            f"flat={len(flat)}/48 linear={len(linear)}/3"
        )
    return points, flat, linear


def padded_y_limits(
    *value_groups,
    padding: float = 0.15,
    nonnegative: bool = False,
) -> tuple[float, float]:
    arrays = [np.asarray(values, dtype=float).ravel() for values in value_groups]
    finite = np.concatenate([values[np.isfinite(values)] for values in arrays if values.size])
    if finite.size == 0:
        raise RuntimeError("cannot derive y-axis limits from empty data")
    lower = np.floor((float(finite.min()) - padding) * 2.0) / 2.0
    if nonnegative:
        lower = max(0.0, lower)
    upper = np.ceil((float(finite.max()) + padding) * 2.0) / 2.0
    if upper <= lower:
        upper = lower + 1.0
    return lower, upper


def write_layout(path: Path, nodes: list[dict], *, min_plot_pt: float) -> None:
    path.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "title_axis_x": 0.045 * SLIDE_WIDTH_PX,
                "minimum_audience_font_px": font_px(14.0),
                "minimum_title_font_px": font_px(31.0),
                "minimum_plot_annotation_font_px": font_px(min_plot_pt),
                "nodes": nodes,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def render_grid_slide(points: pd.DataFrame, flat: pd.DataFrame) -> None:
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[dict] = []
    title = (
        f"R = {RADIUS_TEXT} isolation thresholds per-cent fits, {PT_RANGE_TEXT} GeV"
        if ALGORITHM == "baseline"
        else f"R = {RADIUS_TEXT} topo isolation thresholds per-cent fits, {PT_RANGE_TEXT} GeV"
    )

    add_text(
        fig,
        nodes,
        "title",
        0.045,
        0.955,
        title,
        size=32.0,
        role="title",
        weight="bold",
        va="top",
    )
    handles = [
        Line2D(
            [0],
            [0],
            color=COLORS[eff],
            marker=MARKERS[eff],
            linestyle="--",
            linewidth=1.8,
            markersize=6.5,
            label=f"{LABELS[eff]} efficiency",
        )
        for eff in (0.7, 0.8, 0.9)
    ]
    legend = fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.50, 0.893),
        ncol=3,
        frameon=False,
        fontsize=14.5,
        handlelength=2.0,
        columnspacing=2.1,
    )

    centrality_bins = (
        points[["cent_min", "cent_max"]]
        .drop_duplicates()
        .sort_values(["cent_min", "cent_max"])
        .itertuples(index=False, name=None)
    )
    y_min, y_max = padded_y_limits(
        points["threshold_gev"] - points["threshold_error_gev"],
        points["threshold_gev"] + points["threshold_error_gev"],
        flat["value_gev"] - flat["error_gev"],
        flat["value_gev"] + flat["error_gev"],
        nonnegative=True,
    )
    x_min = float(points["pt_min"].min())
    x_max = float(points["pt_max"].max())
    x_ticks = (
        [15, 20, 25, 30, 35]
        if PT_SCHEME == "nominal"
        else [5, 10, 15, 20, 25, 30, 35, 40]
    )
    gs = fig.add_gridspec(
        4,
        4,
        left=0.066,
        right=0.958,
        bottom=0.165,
        top=0.785,
        wspace=0.13,
        hspace=0.28,
    )
    title_artists = []
    representative_tick = None
    for idx, (cent_min, cent_max) in enumerate(centrality_bins):
        row, col = divmod(idx, 4)
        ax = fig.add_subplot(gs[row, col])
        for eff in (0.7, 0.8, 0.9):
            subset = points[
                (points["cent_min"] == cent_min)
                & (points["cent_max"] == cent_max)
                & np.isclose(points["efficiency"], eff)
            ].sort_values("pt_center")
            fit = flat[
                (flat["cent_min"] == cent_min)
                & (flat["cent_max"] == cent_max)
                & np.isclose(flat["efficiency"], eff)
            ].iloc[0]
            ax.errorbar(
                subset["pt_center"],
                subset["threshold_gev"],
                yerr=subset["threshold_error_gev"],
                color=COLORS[eff],
                marker=MARKERS[eff],
                markersize=3.7,
                linestyle="none",
                elinewidth=1.0,
                capsize=1.6,
                zorder=3,
            )
            ax.axhline(
                float(fit["value_gev"]),
                color=COLORS[eff],
                linestyle="--",
                linewidth=1.25,
                zorder=2,
            )
        title_artists.append(
            ax.set_title(
                f"{int(cent_min)}-{int(cent_max)}%",
                fontsize=14.5,
                fontweight="bold",
                pad=3,
            )
        )
        ax.set_xlim(x_min - 0.3, x_max + 0.3)
        ax.set_ylim(y_min, y_max)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6, min_n_ticks=4))
        ax.set_xticks(x_ticks)
        ax.grid(color=GRID, linewidth=0.65, alpha=0.8)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=13.0, length=2.5, pad=2)
        if col != 0:
            ax.set_yticklabels([])
        if row != 3:
            ax.set_xticklabels([])
        for spine in ax.spines.values():
            spine.set_linewidth(0.9)
            spine.set_color(INK)
        if representative_tick is None and row == 3:
            fig.canvas.draw()
            labels = [label for label in ax.get_xticklabels() if label.get_text()]
            representative_tick = labels[0] if labels else None

    add_text(
        fig,
        nodes,
        "shared x axis label",
        0.515,
        0.112,
        r"photon $p_T$ [GeV]",
        size=17.0,
        role="plot_annotation",
        ha="center",
    )
    add_text(
        fig,
        nodes,
        "shared y axis label",
        0.027,
        0.475,
        r"$E_T^{\mathrm{iso,cut}}$ [GeV]",
        size=17.0,
        role="plot_annotation",
        ha="center",
        rotation=90,
    )
    add_text(
        fig,
        nodes,
        "interpretation",
        0.068,
        0.047,
        r"Many centrality bins show residual $p_T$ dependence; the flat fits are a compact approximation.",
        size=16.0,
        weight="bold",
        color=INK,
        title_band_exception=True,
    )

    fig.canvas.draw()
    add_artist_node(fig, nodes, "grid panel title", title_artists[0], size=14.5, text="0-5%")
    if representative_tick is not None:
        add_artist_node(fig, nodes, "grid tick labels", representative_tick, size=13.0, text="representative axis tick")
    legend_text = legend.get_texts()[0]
    add_artist_node(fig, nodes, "grid legend", legend_text, size=14.5, text="70%, 80%, 90% efficiency")
    require_vertical_symmetry(
        nodes,
        upper_name="title",
        middle_name="grid legend",
        lower_name="grid panel title",
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(GRID_PNG, dpi=SLIDE_DPI)
    write_layout(GRID_LAYOUT, nodes, min_plot_pt=13.0)
    plt.close(fig)


def render_centrality_slide(flat: pd.DataFrame, linear: pd.DataFrame) -> None:
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    nodes: list[dict] = []
    title = (
        f"R = {RADIUS_TEXT}, {PT_RANGE_TEXT} GeV, per cent isolation fits"
        if ALGORITHM == "baseline"
        else f"R = {RADIUS_TEXT}, {PT_RANGE_TEXT} GeV, topo per cent isolation fits"
    )
    algorithm_label = "tower-cone baseline" if ALGORITHM == "baseline" else "topocluster isolation"

    add_text(
        fig,
        nodes,
        "title",
        0.045,
        0.955,
        title,
        size=32.0,
        role="title",
        weight="bold",
        va="top",
    )
    ax = fig.add_axes([0.068, 0.145, 0.620, 0.675])
    centrality = np.linspace(0.0, 80.0, 300)
    line_endpoints = []
    for fit in linear.itertuples():
        line_endpoints.extend(
            [
                float(fit.intercept_gev),
                float(fit.intercept_gev) + 80.0 * float(fit.slope_gev_per_percent),
            ]
        )
    y_min, y_max = padded_y_limits(
        flat["value_gev"] - flat["error_gev"],
        flat["value_gev"] + flat["error_gev"],
        line_endpoints,
    )
    for eff in (0.7, 0.8, 0.9):
        subset = flat[np.isclose(flat["efficiency"], eff)].sort_values("cent_center")
        fit = linear[np.isclose(linear["efficiency"], eff)].iloc[0]
        ax.errorbar(
            subset["cent_center"],
            subset["value_gev"],
            yerr=subset["error_gev"],
            color=COLORS[eff],
            marker=MARKERS[eff],
            markersize=7.0,
            linestyle="none",
            elinewidth=1.35,
            capsize=2.2,
            label=f"{LABELS[eff]} efficiency",
            zorder=3,
        )
        ax.plot(
            centrality,
            float(fit["intercept_gev"]) + float(fit["slope_gev_per_percent"]) * centrality,
            color=COLORS[eff],
            linewidth=2.0,
            zorder=2,
        )
    ax.set_xlim(0, 80)
    ax.set_ylim(y_min, y_max)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=7, min_n_ticks=5))
    ax.set_xlabel("Centrality [%]", fontsize=16.0, labelpad=7)
    ax.set_ylabel(r"Flat-fit $E_T^{\mathrm{iso,cut}}$ [GeV]", fontsize=16.0, labelpad=9)
    ax.tick_params(labelsize=13.0)
    ax.grid(color=GRID, linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)
    ax.legend(loc="upper right", frameon=False, fontsize=12.8, handlelength=1.4)
    ax.set_title(
        f"Embedded photons | R={RADIUS_TEXT} {algorithm_label} | {PT_SCHEME} {PT_RANGE_TEXT} GeV",
        loc="left",
        fontsize=14.5,
        pad=10,
        color=MUTED,
    )
    internal_word = ax.text(
        0.972,
        0.790,
        "Internal",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontfamily="Times New Roman",
        fontstyle="normal",
        fontweight="normal",
        fontsize=13.5,
        color=INK,
    )
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    internal_bbox = internal_word.get_window_extent(renderer=renderer)
    sphenix_right_x = ax.transAxes.inverted().transform((internal_bbox.x0 - 6.0, internal_bbox.y0))[0]
    sphenix_word = ax.text(
        sphenix_right_x,
        0.790,
        "sPHENIX",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontfamily="Times New Roman",
        fontstyle="italic",
        fontweight="bold",
        fontsize=13.5,
        color=INK,
    )
    sample_label = ax.text(
        0.972,
        0.742,
        "Embedded Photon+Jet 12+20",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=12.8,
        color=INK,
    )
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)
        spine.set_color(INK)

    fig.patches.append(
        Rectangle(
            (0.724, 0.145),
            0.0025,
            0.675,
            transform=fig.transFigure,
            facecolor=EDGE,
            edgecolor="none",
            zorder=3,
        )
    )
    add_text(
        fig,
        nodes,
        "fit test heading",
        0.742,
        0.755,
        "Per centrality fit functions",
        size=19.7,
        weight="bold",
    )
    add_text(
        fig,
        nodes,
        "fit definition",
        0.755,
        0.695,
        r"$E_T^{\mathrm{iso,cut}} = a + bC$",
        size=17.0,
        color=MUTED,
    )

    for separator_y in (0.485, 0.305):
        fig.patches.append(
            Rectangle(
                (0.755, separator_y),
                0.190,
                0.0012,
                transform=fig.transFigure,
                facecolor=EDGE,
                edgecolor="none",
                zorder=3,
            )
        )

    y_positions = [0.610, 0.430, 0.250]
    for eff, y in zip((0.7, 0.8, 0.9), y_positions):
        fit = linear[np.isclose(linear["efficiency"], eff)].iloc[0]
        equation = (
            rf"{float(fit['intercept_gev']):.3f} "
            rf"{float(fit['slope_gev_per_percent']):+.4f}$C$"
        )
        fit_quality = rf"$\chi^2$/ndf = {float(fit['chi2']):.1f}/{int(fit['ndf'])}"
        add_text(
            fig,
            nodes,
            f"{LABELS[eff]} efficiency label",
            0.755,
            y,
            f"{LABELS[eff]} efficiency",
            size=16.8,
            weight="bold",
            color=COLORS[eff],
        )
        add_text(
            fig,
            nodes,
            f"{LABELS[eff]} equation",
            0.755,
            y - 0.045,
            equation,
            size=15.8,
            color=INK,
        )
        add_text(
            fig,
            nodes,
            f"{LABELS[eff]} fit quality",
            0.755,
            y - 0.085,
            fit_quality,
            size=14.2,
            color=MUTED,
        )

    fig.canvas.draw()
    add_artist_node(
        fig,
        nodes,
        "sPHENIX plot label",
        sphenix_word,
        size=13.5,
        text="sPHENIX",
    )
    add_artist_node(
        fig,
        nodes,
        "internal plot label",
        internal_word,
        size=13.5,
        text="Internal",
    )
    add_artist_node(
        fig,
        nodes,
        "embedded sample plot label",
        sample_label,
        size=12.8,
        text="Embedded Photon+Jet 12+20",
    )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(CENT_PNG, dpi=SLIDE_DPI)
    write_layout(CENT_LAYOUT, nodes, min_plot_pt=12.0)
    plt.close(fig)


def write_sidecars(linear: pd.DataFrame) -> None:
    algorithm_spoken = "tower-cone baseline" if ALGORITHM == "baseline" else "topocluster isolation"
    GRID_SCRIPT.write_text(
        "\n".join(
            [
                f"# Speaker Script: R={RADIUS_TEXT} isolation thresholds across photon pT",
                "",
                f"This slide shows the R equals {RADIUS_TEXT} {algorithm_spoken} threshold derived from embedded truth-matched photons.",
                f"Each panel is one five-percent centrality interval. The {PT_BIN_COUNT} points span the {PT_SCHEME} photon-pT range from {PT_RANGE_TEXT.replace('-', ' to ')} GeV, and the dashed lines are constant fits for 70, 80, and 90 percent signal efficiency.",
                "The dominant trend is centrality: the required isolation threshold falls steadily from central to peripheral events.",
                "Within an individual centrality bin the pT variation is modest. Many, but not all, of the constant fits show statistically significant residual pT dependence.",
                "So these flat fits are a useful compression of each bin, but they should be treated as an approximation rather than proof of exact pT independence.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    qualities = {
        LABELS[float(row.efficiency)]: f"{float(row.chi2):.1f}/{int(row.ndf)}"
        for row in linear.itertuples()
    }
    CENT_SCRIPT.write_text(
        "\n".join(
            [
                "# Speaker Script: Centrality parameterization test",
                "",
                "This slide takes the flat-fit threshold from each five-percent centrality bin and tests the simplest possible global parameterization: one straight line from zero to eighty percent centrality.",
                "The thresholds clearly decrease toward peripheral events, but the large fit residuals show that one straight line is not an adequate global description.",
                f"The chi-squared per degree-of-freedom values are {qualities['70%']} for 70 percent efficiency, {qualities['80%']} for 80 percent, and {qualities['90%']} for 90 percent.",
                "Those values make the conclusion straightforward: the global linear form is not adequate for the final baseline.",
                "The next comparison should test a smooth curved function, a piecewise form, and direct lookup values in each five-percent centrality bin. We should freeze the prescription only after that comparison.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    source_manifest = json.loads(SOURCE_MANIFEST.read_text(encoding="utf-8"))
    OUT_MANIFEST.write_text(
        json.dumps(
            {
                "schema": f"{CAMPAIGN_LABEL}_AUAU_{CONE}_{ALGORITHM.upper()}_ISOLATION_SLIDES_V1",
                "campaign_tag": source_manifest.get("campaign_tag"),
                "sample": "embedded Photon12+20 truth-matched signal",
                "algorithm": (
                    f"R={RADIUS_TEXT} tower-cone baseline without topocluster isolation"
                    if ALGORITHM == "baseline"
                    else f"R={RADIUS_TEXT} topocluster isolation"
                ),
                "pt_range": f"{PT_SCHEME} {PT_RANGE_TEXT} GeV in {PT_BIN_COUNT} bins",
                "centrality": "16 bins of width 5% spanning 0-80%",
                "efficiencies": [0.70, 0.80, 0.90],
                "inputs": {
                    "cutoff_points": str(POINTS_CSV),
                    "flat_fits": str(FLAT_CSV),
                    "linear_fits": str(LINEAR_CSV),
                    "source_manifest": str(SOURCE_MANIFEST),
                    "source_root": source_manifest.get("input_root"),
                },
                "slides": [
                    {
                        "png": str(GRID_PNG),
                        "speaker_script": str(GRID_SCRIPT),
                        "layout_nodes": str(GRID_LAYOUT),
                        "audit": str(GRID_AUDIT),
                        "claim": "centrality drives the large threshold shift; within-bin pT variation is modest but statistically visible",
                    },
                    {
                        "png": str(CENT_PNG),
                        "speaker_script": str(CENT_SCRIPT),
                        "layout_nodes": str(CENT_LAYOUT),
                        "audit": str(CENT_AUDIT),
                        "claim": "one global linear centrality model is statistically inadequate",
                    },
                ],
                "linear_fit_chi2_ndf": {
                    LABELS[float(row.efficiency)]: [float(row.chi2), int(row.ndf)]
                    for row in linear.itertuples()
                },
                "mutation_boundary": "local PNG candidates only; no Google Slides mutation",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def run_audits() -> None:
    audit = SCRIPTS / "slides/common/post_render_slide_audit.py"
    for png, layout, output in (
        (GRID_PNG, GRID_LAYOUT, GRID_AUDIT),
        (CENT_PNG, CENT_LAYOUT, CENT_AUDIT),
    ):
        subprocess.run(
            [
                sys.executable,
                str(audit),
                "--png",
                str(png),
                "--layout-nodes",
                str(layout),
                "--output",
                str(output),
                "--require-pass",
            ],
            check=True,
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cone", choices=("R30", "R40"), default="R30")
    parser.add_argument("--algorithm", choices=("baseline", "topocluster"), default="baseline")
    parser.add_argument("--pt-scheme", choices=("nominal", "expanded"), default="nominal")
    parser.add_argument("--source-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--source-prefix-root", default="the96")
    parser.add_argument("--campaign-label", default="THE96")
    args = parser.parse_args()
    configure_paths(
        args.cone,
        args.algorithm,
        args.pt_scheme,
        source_dir_override=args.source_dir,
        output_dir_override=args.output_dir,
        source_prefix_root=args.source_prefix_root,
        campaign_label=args.campaign_label,
    )
    setup_style()
    points, flat, linear = load_tables()
    render_grid_slide(points, flat)
    render_centrality_slide(flat, linear)
    write_sidecars(linear)
    run_audits()
    print(f"Wrote {GRID_PNG}")
    print(f"Wrote {CENT_PNG}")
    print(f"Wrote {GRID_SCRIPT}")
    print(f"Wrote {CENT_SCRIPT}")
    print(f"Wrote {OUT_MANIFEST}")


if __name__ == "__main__":
    main()

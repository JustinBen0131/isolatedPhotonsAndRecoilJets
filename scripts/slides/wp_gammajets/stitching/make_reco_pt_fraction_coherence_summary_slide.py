#!/usr/bin/env python3
"""Build a one-slide summary of reco-cluster pT fraction coherence checks."""

from __future__ import annotations

import importlib.util
import json
from io import BytesIO
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
from PIL import Image, ImageDraw, ImageFont


REPO = Path(__file__).resolve().parents[4]
OUT = REPO / "outputs/manual-20260612-working-point-slide2-3-reco-pt-fraction-summary"

EMBEDDED_SCRIPT = REPO / "scripts/slides/wp_gammajets/stitching/make_reco_cluster_et_leakage_followup_slide.py"
PP_SCRIPT = REPO / "scripts/slides/pp_currentian/leakage/make_pp_reco_cluster_et_leakage_followup_slide.py"
PP_ABCD_CAP_BASE = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_currentIAN_recoEt50NPBMax50_canonicalTree_20260606_152337"
    / "validation/reco_cluster_et_leakage_abcdMatchedFixedIso2_sampleCap"
)
PP_ABCD_CAP_CSV = PP_ABCD_CAP_BASE / "pp_reco_cluster_et_leakage_components.csv"
PP_ABCD_CAP_SUMMARY = PP_ABCD_CAP_BASE / "pp_reco_cluster_et_leakage_summary.json"

W, H = 2560, 1440
FONT_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES = FONT_DIR / "Times New Roman.ttf"
TIMES_BOLD = FONT_DIR / "Times New Roman Bold.ttf"

INK = (18, 24, 38)
MUTED = (78, 86, 101)
LINE = (207, 215, 228)
PANEL_LINE = (196, 190, 226)
PANEL_FILL = (255, 255, 255)
SOFT = (248, 249, 252)
TEAL = (26, 111, 116)
TEAL_SOFT = (235, 248, 248)
TEAL_LINE = (174, 214, 216)
PURPLE = (92, 78, 152)
AU_AU_RED = (218, 61, 50)
PP_BLUE = (36, 74, 210)


def font(size: int, *, bold: bool = False) -> ImageFont.FreeTypeFont:
    path = TIMES_BOLD if bold else TIMES
    try:
        return ImageFont.truetype(str(path), size)
    except OSError:
        return ImageFont.load_default()


F = {
    "title": font(86, bold=True),
    "subtitle": font(37),
    "panel_title": font(44, bold=True),
    "panel_subtitle": font(37),
    "explainer": font(62),
    "panel_label": font(50, bold=True),
    "body": font(56),
    "body_bold": font(56, bold=True),
    "small": font(37),
    "small_bold": font(37, bold=True),
}


def load_module(path: Path, name: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def text_bbox(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, fnt: ImageFont.FreeTypeFont) -> tuple[int, int, int, int]:
    bbox = draw.textbbox(xy, text, font=fnt)
    return tuple(int(v) for v in bbox)


def rounded(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], fill, outline=LINE, width=2, radius=18) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def paste_contained(canvas: Image.Image, img: Image.Image, box: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    x0, y0, x1, y1 = box
    bw, bh = x1 - x0, y1 - y0
    scale = min(bw / img.width, bh / img.height)
    nw, nh = int(img.width * scale), int(img.height * scale)
    resized = img.resize((nw, nh), Image.Resampling.LANCZOS)
    bg = Image.new("RGB", (bw, bh), "white")
    ox, oy = (bw - nw) // 2, (bh - nh) // 2
    bg.paste(resized, (ox, oy))
    canvas.paste(bg, (x0, y0))
    return (x0 + ox, y0 + oy, x0 + ox + nw, y0 + oy + nh)


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    message: str,
    fnt: ImageFont.FreeTypeFont,
    *,
    max_width: int,
    fill=INK,
    line_gap: int = 8,
) -> tuple[int, int, int, int]:
    words = message.split()
    lines: list[str] = []
    current = ""
    for word in words:
        trial = f"{current} {word}".strip()
        if draw.textlength(trial, font=fnt) <= max_width or not current:
            current = trial
        else:
            lines.append(current)
            current = word
    if current:
        lines.append(current)

    top = y
    right = x
    for line in lines:
        draw.text((x, y), line, font=fnt, fill=fill)
        bbox = draw.textbbox((x, y), line, font=fnt)
        right = max(right, int(bbox[2]))
        y += fnt.size + line_gap
    return (x, top, right, y - line_gap)


def draw_rich_line(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    label: str,
    body: str,
    *,
    fill=INK,
    label_fill=INK,
    body_font: ImageFont.FreeTypeFont | None = None,
    label_font: ImageFont.FreeTypeFont | None = None,
) -> tuple[tuple[int, int, int, int], dict[str, Any]]:
    label_font = label_font or F["body_bold"]
    body_font = body_font or F["body"]
    draw.text((x, y), label, font=label_font, fill=label_fill)
    label_w = int(draw.textlength(label, font=label_font))
    draw.text((x + label_w, y), body, font=body_font, fill=fill)
    bbox = draw.textbbox((x, y), label + body, font=body_font)
    bbox = (x, int(bbox[1]), x + label_w + int(draw.textlength(body, font=body_font)), int(bbox[3]))
    node = {
        "text": label + body,
        "text_runs": [
            {"text": label, "bold": True},
            {"text": body, "bold": False},
        ],
    }
    return bbox, node


def draw_panel_header(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    title: str,
    *,
    accent: tuple[int, int, int],
) -> tuple[int, int, int, int]:
    x0, y0, _, _ = box
    draw.rounded_rectangle((x0 + 20, y0 + 20, x0 + 32, y0 + 88), radius=6, fill=accent)
    draw.text((x0 + 52, y0 + 24), title, font=F["panel_title"], fill=INK)
    return text_bbox(draw, (x0 + 52, y0 + 24), title, F["panel_title"])


def draw_top_canvas_label(
    draw: ImageDraw.ImageDraw,
    plot_ink: tuple[int, int, int, int],
    lines: tuple[str, ...],
    *,
    accent: tuple[int, int, int],
) -> tuple[int, int, int, int]:
    x0, y0, _, _ = plot_ink
    x = x0 + 170
    y = y0 + 326
    line_gap = -2
    boxes: list[tuple[int, int, int, int]] = []
    for line in lines:
        bbox = draw.textbbox((x, y), line, font=F["panel_label"])
        draw.text((x, y), line, font=F["panel_label"], fill=accent, stroke_width=4, stroke_fill=(255, 255, 255))
        boxes.append((int(bbox[0] - 4), int(bbox[1] - 4), int(bbox[2] + 4), int(bbox[3] + 4)))
        y += F["panel_label"].size + line_gap
    return (
        min(box[0] for box in boxes),
        min(box[1] for box in boxes),
        max(box[2] for box in boxes),
        max(box[3] for box in boxes),
    )


def draw_plot_border(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], *, outline: tuple[int, int, int]) -> tuple[int, int, int, int]:
    draw.rounded_rectangle(box, radius=10, outline=outline, width=4)
    return box


def measure_arrow_bullet(
    draw: ImageDraw.ImageDraw,
    *,
    x: int,
    y: int,
    text: str,
    fnt: ImageFont.FreeTypeFont,
) -> tuple[int, int, int, int]:
    marker_half = max(11, int(round(fnt.size * 0.25)))
    marker_width = marker_half * 2
    text_x = x + marker_width + 22
    bbox = draw.textbbox((text_x, y), text, font=fnt)
    cy = (bbox[1] + bbox[3]) // 2
    return (x, min(cy - marker_half, bbox[1]), int(bbox[2]), max(cy + marker_half, bbox[3]))


def draw_arrow_bullet(
    draw: ImageDraw.ImageDraw,
    *,
    x: int,
    y: int,
    text: str,
    fnt: ImageFont.FreeTypeFont,
    marker_fill: tuple[int, int, int] = TEAL,
    text_fill: tuple[int, int, int] = INK,
) -> tuple[int, int, int, int]:
    marker_half = max(11, int(round(fnt.size * 0.25)))
    marker_width = marker_half * 2
    text_x = x + marker_width + 22
    bbox = draw.textbbox((text_x, y), text, font=fnt)
    cy = (bbox[1] + bbox[3]) // 2
    draw.polygon([(x, cy - marker_half), (x, cy + marker_half), (x + marker_width, cy)], fill=marker_fill)
    draw.text((text_x, y), text, font=fnt, fill=text_fill)
    return (x, min(cy - marker_half, bbox[1]), int(bbox[2]), max(cy + marker_half, bbox[3]))


def render_composition_plot(
    pivot: list[dict[str, Any]],
    *,
    samples: tuple[str, ...],
    colors: dict[str, tuple[int, int, int]],
    markers: dict[str, str],
    source_line: str,
    windows_text: str | None,
    display_range: tuple[float, float] = (12.0, 50.0),
) -> Image.Image:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogLocator, NullFormatter

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.35,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig = plt.figure(figsize=(10.55, 7.20), dpi=120, facecolor="white")
    gs = fig.add_gridspec(
        2,
        1,
        height_ratios=[3.3, 1.17],
        left=0.112,
        right=0.982,
        bottom=0.112,
        top=0.972,
        hspace=0.060,
    )
    ax = fig.add_subplot(gs[0])
    fax = fig.add_subplot(gs[1], sharex=ax)

    lo, hi = display_range
    display = [row for row in pivot if lo <= float(row["x_center"]) <= hi]
    xs = np.array([r["x_center"] for r in display], dtype=float)
    xerr = np.array([[r["x_err_low"] for r in display], [r["x_err_high"] for r in display]], dtype=float)

    for sample in (*samples, "Sum"):
        y = np.array([r[f"weighted_entries_{sample}"] for r in display], dtype=float)
        ey = np.array([r[f"weighted_error_{sample}"] for r in display], dtype=float)
        mask = y > 0.0
        color = np.array(colors[sample]) / 255.0
        ax.errorbar(
            xs[mask],
            y[mask],
            yerr=ey[mask],
            xerr=xerr[:, mask],
            fmt=markers[sample],
            markersize=5.9 if sample != "Sum" else 6.4,
            markerfacecolor=color if sample != "Sum" else "white",
            markeredgecolor=color,
            markeredgewidth=1.05,
            linestyle="none",
            color=color,
            ecolor=color,
            elinewidth=1.0,
            capsize=1.8,
            label=sample if sample != "Sum" else "Weighted sum",
            zorder=5 if sample == "Sum" else 4,
        )

    for sample in samples:
        color = np.array(colors[sample]) / 255.0
        vals = np.array([r[f"fraction_{sample}"] for r in display], dtype=float)
        mask = np.array([r["weighted_entries_Sum"] for r in display], dtype=float) > 0.0
        fax.plot(
            xs[mask],
            vals[mask],
            marker=markers[sample],
            markersize=6.0,
            linestyle="none",
            color=color,
            label=sample,
            zorder=3,
        )

    positive = [
        float(row[f"weighted_entries_{sample}"])
        for row in display
        for sample in (*samples, "Sum")
        if float(row[f"weighted_entries_{sample}"]) > 0.0
    ]
    ax.set_yscale("log")
    ax.set_xlim(*display_range)
    ax.set_ylim(max(1.0e-5, min(positive) * 0.45), max(positive) * 2.4)
    ax.set_ylabel("Weighted entries / 1 GeV bin", fontsize=22.0)
    ax.tick_params(labelbottom=False, labelsize=18.0, length=6)
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.grid(which="major", color="0.84", linestyle=":", linewidth=0.75)
    ax.grid(which="minor", color="0.91", linestyle=":", linewidth=0.45)

    fax.set_ylim(-0.035, 1.05)
    fax.set_yticks([0.0, 0.25, 0.50, 0.75, 1.0])
    fax.set_ylabel("Fraction of\nweighted sum", fontsize=21.0)
    fax.set_xlabel(r"Reco photon-cluster $E_T$ ($p_T^\gamma$) [GeV]", fontsize=22.0, loc="right")
    fax.tick_params(labelsize=18.0, length=6)
    fax.grid(which="major", color="0.84", linestyle=":", linewidth=0.75)

    ax.text(
        0.970,
        0.985,
        r"$\it{\bf{sPHENIX}}$ Internal" + "\n" + source_line,
        transform=ax.transAxes,
        fontsize=20.0,
        ha="right",
        va="top",
        linespacing=1.05,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 2.0},
        zorder=20,
    )
    if windows_text:
        ax.text(
            0.040,
            0.085,
            windows_text,
            transform=ax.transAxes,
            fontsize=18.0,
            ha="left",
            va="bottom",
            color=np.array(INK) / 255.0,
            linespacing=1.16,
            bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": "0.82", "alpha": 0.90, "linewidth": 0.8},
            zorder=20,
        )
    legend = ax.legend(
        loc="upper right",
        bbox_to_anchor=(0.965, 0.760),
        frameon=False,
        fontsize=17.0,
        handlelength=1.10,
        labelspacing=0.16,
        columnspacing=0.70,
        borderaxespad=0.2,
        markerscale=1.0,
        ncol=3,
    )
    for item in legend.get_texts():
        item.set_fontsize(17.0)

    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=120, facecolor="white")
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGB")


def make_plot_images() -> tuple[Image.Image, Image.Image]:
    embedded = load_module(EMBEDDED_SCRIPT, "embedded_reco_cluster_fraction")
    pp = load_module(PP_SCRIPT, "pp_reco_cluster_fraction")
    embedded_plot = render_composition_plot(
        embedded.pivot_rows(embedded.load_rows()),
        samples=embedded.SAMPLES,
        colors=embedded.COLORS,
        markers=embedded.MARKERS,
        source_line="PYTHIA8 embedded inclusive jet, 0-80%",
        windows_text=None,
    )
    pp_plot = render_composition_plot(
        pp.pivot_rows(pp.load_rows(PP_ABCD_CAP_CSV)),
        samples=pp.SAMPLES,
        colors=pp.COLORS,
        markers=pp.MARKERS,
        source_line="PYTHIA8 pp inclusive jet, 200 GeV",
        windows_text=None,
    )
    return embedded_plot.convert("RGB"), pp_plot.convert("RGB")


def plot_annotation_nodes(prefix: str, ink: tuple[int, int, int, int]) -> list[dict[str, Any]]:
    x0, y0, x1, y1 = ink
    w = x1 - x0
    h = y1 - y0
    return [
        {
            "name": f"{prefix} source label annotation",
            "kind": "text",
            "role": "plot_annotation",
            "font_px": 29,
            "bbox": [x0 + int(0.64 * w), y0 + int(0.02 * h), x1 - int(0.02 * w), y0 + int(0.11 * h)],
            "text": "sPHENIX Internal / PYTHIA8 source label",
        },
        {
            "name": f"{prefix} legend annotation",
            "kind": "text",
            "role": "plot_annotation",
            "font_px": 29,
            "bbox": [x0 + int(0.50 * w), y0 + int(0.25 * h), x1 - int(0.03 * w), y0 + int(0.38 * h)],
            "text": "sample legend",
        },
        {
            "name": f"{prefix} fraction axis annotation",
            "kind": "text",
            "role": "plot_annotation",
            "font_px": 29,
            "bbox": [x0 + int(0.01 * w), y0 + int(0.69 * h), x0 + int(0.13 * w), y0 + int(0.90 * h)],
            "text": "fraction axis label",
        },
        {
            "name": f"{prefix} reco cluster axis annotation",
            "kind": "text",
            "role": "plot_annotation",
            "font_px": 31,
            "bbox": [x0 + int(0.50 * w), y0 + int(0.91 * h), x1, y1],
            "text": "reco photon-cluster axis label",
        },
    ]


def write_speaker_script(path: Path) -> None:
    path.write_text(
        """# Speaker Script: Updated stitching makes reco-cluster pT fractions coherent

This is the follow-up to last week's stitching question.

The important update is that the ownership windows and the pp row scope have been corrected. With those updates, each PYTHIA input sample's share of the weighted reco-cluster pT spectrum now looks consistent between the embedded inclusive-jet check and the pp ABCD-matched check.

The key thing to read is the lower fraction panel in each plot. It shows how much each input sample contributes as reco-cluster pT increases, and the handoff across samples now looks in line rather than obviously mismatched.

The pp legend has one extra low-pT sample because that source set starts with jet8, while the embedded check starts at Jet12; the comparison point is the smooth handoff pattern, not identical sample membership.
""",
        encoding="utf-8",
    )


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    png_path = OUT / "reco_pt_fraction_coherence_summary.png"
    layout_path = OUT / "layout_nodes.json"
    manifest_path = OUT / "manifest.json"
    script_path = OUT / "speaker_script.md"

    embedded_plot, pp_plot = make_plot_images()

    canvas = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(canvas)
    nodes: list[dict[str, Any]] = []

    title = "Updated stitching gives coherent reco-cluster pT fractions"
    draw.text((74, 55), title, font=F["title"], fill=INK)
    title_bbox = text_bbox(draw, (74, 55), title, F["title"])
    nodes.append({"name": "slide title", "kind": "text", "role": "title", "font_px": 86, "bbox": list(title_bbox), "text": title})

    title_axis_x = 74
    left_panel = (74, 255, 1275, 1058)
    right_panel = (1285, 255, 2486, 1058)
    nodes.extend(
        [
            {
                "name": "embedded plot slot",
                "kind": "panel",
                "bbox": list(left_panel),
                "symmetry_group": "main_plot_panels",
                "title_axis_align": "left",
            },
            {
                "name": "pp plot slot",
                "kind": "panel",
                "bbox": list(right_panel),
                "symmetry_group": "main_plot_panels",
            },
        ]
    )

    explainer = "Lower panels show each PYTHIA input sample's share of the weighted reco-cluster pT spectrum"
    draw.text((74, 165), explainer, font=F["explainer"], fill=INK)
    explainer_bbox = text_bbox(draw, (74, 165), explainer, F["explainer"])
    nodes.append(
        {
            "name": "shared source fraction explanation",
            "kind": "text",
            "role": "audience",
            "font_px": 62,
            "bbox": list(explainer_bbox),
            "text": explainer,
            "intentional_alignment": "top",
            "title_axis_align": "left",
        }
    )

    left_plot_box = (94, 260, 1255, 1052)
    right_plot_box = (1305, 260, 2466, 1052)
    left_ink = paste_contained(canvas, embedded_plot, left_plot_box)
    right_ink = paste_contained(canvas, pp_plot, right_plot_box)
    left_frame = left_panel
    right_frame = right_panel
    left_border = draw_plot_border(draw, left_frame, outline=AU_AU_RED)
    right_border = draw_plot_border(draw, right_frame, outline=PP_BLUE)
    auau_label_bbox = draw_top_canvas_label(draw, left_ink, ("AuAu", "(emb PYTHIA)"), accent=AU_AU_RED)
    pp_label_bbox = draw_top_canvas_label(draw, right_ink, ("pp", "(PYTHIA)"), accent=PP_BLUE)
    nodes.extend(
        [
            {"name": "embedded plot usable box", "kind": "figure_box", "bbox": list(left_plot_box)},
            {
                "name": "embedded plot ink",
                "kind": "figure_ink",
                "parent": "embedded plot usable box",
                "outer_frame": "embedded plot red border",
                "max_outer_frame_top_gap_px": 28,
                "bbox": list(left_ink),
            },
            {
                "name": "embedded plot red border",
                "kind": "panel",
                "bbox": list(left_border),
                "documented_exception": "intentional AuAu-vs-pp color identity from Justin sketch",
                "symmetry_group": "plot_color_frames",
                "edge_color": "#da3d32",
                "color_difference_reason": "semantic AuAu-vs-pp identity",
                "title_axis_align": "left",
            },
            {"name": "pp plot usable box", "kind": "figure_box", "bbox": list(right_plot_box)},
            {
                "name": "pp plot ink",
                "kind": "figure_ink",
                "parent": "pp plot usable box",
                "outer_frame": "pp plot blue border",
                "max_outer_frame_top_gap_px": 28,
                "bbox": list(right_ink),
            },
            {
                "name": "pp plot blue border",
                "kind": "panel",
                "bbox": list(right_border),
                "documented_exception": "intentional AuAu-vs-pp color identity from Justin sketch",
                "symmetry_group": "plot_color_frames",
                "edge_color": "#244ad2",
                "color_difference_reason": "semantic AuAu-vs-pp identity",
            },
            {
                "name": "AuAu sample label",
                "kind": "text",
                "role": "audience",
                "font_px": 50,
                "bbox": list(auau_label_bbox),
                "text": "AuAu\n(emb PYTHIA)",
                "documented_exception": "intentional top-canvas identity label, reduced and lowered to avoid data overlap",
            },
            {
                "name": "pp sample label",
                "kind": "text",
                "role": "audience",
                "font_px": 50,
                "bbox": list(pp_label_bbox),
                "text": "pp\n(PYTHIA)",
                "documented_exception": "intentional top-canvas identity label, reduced and lowered to avoid data overlap",
            },
        ]
    )
    nodes.extend(plot_annotation_nodes("embedded plot", left_ink))
    nodes.extend(plot_annotation_nodes("pp plot", right_ink))

    result_box = (74, 1080, 2486, 1366)
    rounded(draw, result_box, (255, 255, 255), outline=TEAL_LINE, width=3, radius=18)
    nodes.append({"name": "takeaway card", "kind": "card", "bbox": list(result_box), "title_axis_align": "left"})
    draw.rounded_rectangle((result_box[0], result_box[1] + 26, result_box[0] + 14, result_box[3] - 26), radius=7, fill=TEAL)

    body_lines = [
        "Follow-up to last week with updated ownership windows fixing the truth-spectrum stitching.",
        "Input-sample shares now look in line across embedded and pp.",
    ]
    card_text_boxes: list[tuple[int, int, int, int]] = []
    line_gap = 62
    bullet_x = result_box[0] + 68
    relative_y = 0
    measured_boxes: list[tuple[int, int, int, int]] = []
    line_heights: list[int] = []
    for line in body_lines:
        bbox = measure_arrow_bullet(draw, x=bullet_x, y=relative_y, text=line, fnt=F["body"])
        measured_boxes.append(bbox)
        line_heights.append(bbox[3] - bbox[1])
        relative_y += line_heights[-1] + line_gap
    measured_block = (
        min(box[0] for box in measured_boxes),
        min(box[1] for box in measured_boxes),
        max(box[2] for box in measured_boxes),
        max(box[3] for box in measured_boxes),
    )
    card_center_y = (result_box[1] + result_box[3]) / 2
    measured_center_y = (measured_block[1] + measured_block[3]) / 2
    text_y = int(round(card_center_y - measured_center_y))
    for idx, line in enumerate(body_lines):
        bbox = draw_arrow_bullet(draw, x=bullet_x, y=text_y, text=line, fnt=F["body"])
        card_text_boxes.append(bbox)
        nodes.append(
            {
                "name": f"takeaway line {idx + 1}",
                "kind": "text",
                "role": "audience",
                "font_px": 56,
                "bbox": list(bbox),
                "text": line,
                "intentional_alignment": "top",
            }
        )
        text_y += line_heights[idx] + line_gap
    text_block = (
        min(box[0] for box in card_text_boxes),
        min(box[1] for box in card_text_boxes),
        max(box[2] for box in card_text_boxes),
        max(box[3] for box in card_text_boxes),
    )
    nodes.append(
        {
            "name": "takeaway text block",
            "kind": "text",
            "role": "audience",
            "parent": "takeaway card",
            "bbox": list(text_block),
            "font_px": 56,
            "text": "\n".join(body_lines),
            "vertical_fill_min_ratio": 0.54,
            "vertical_fill_max_ratio": 0.68,
        }
    )

    canvas.save(png_path)
    write_speaker_script(script_path)

    layout = {
        "schema_version": 1,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "slide_size": [W, H],
        "title_axis_x": title_axis_x,
        "title_axis_tolerance_px": 3,
        "vertical_margin_balance": {
            "top_node": "slide title",
            "bottom_node": "takeaway card",
            "target_gap_px": title_axis_x,
            "tolerance_px": 4,
        },
        "nodes": nodes,
    }
    layout_path.write_text(json.dumps(layout, indent=2) + "\n", encoding="utf-8")

    manifest = {
        "generated_at": layout["generated_at"],
        "artifact": str(png_path.relative_to(REPO)),
        "deck_context": {
            "presentation_id": "167x-He2rOOBO2i4nNS6Pdcqu7Wv03GeFMuWH9tRRx-8",
            "source_slides": [2, 3],
            "mutation": "none; local PNG candidate only",
        },
        "sources": {
            "embedded_script": str(EMBEDDED_SCRIPT.relative_to(REPO)),
            "pp_script": str(PP_SCRIPT.relative_to(REPO)),
            "embedded_csv": "dataOutput/stitchDiagnostics/focus21_clusterEt_leakage_jet1234_20260526/inclusive_jet1234_reco_cluster_et_weighted_components_fine1gev12to50.csv",
            "pp_csv": str(PP_ABCD_CAP_CSV.relative_to(REPO)),
            "pp_summary": str(PP_ABCD_CAP_SUMMARY.relative_to(REPO)),
        },
        "interpretation_notes": [
            "The pp panel intentionally includes jet8 while the embedded panel starts at Jet12; the source sample sets differ.",
            "In-plot ownership-window text from the source slides was omitted here because the combined slide explains the updated ownership/row-scope in the bottom card and needs uncluttered plot panels.",
            "The red AuAu label/border and blue pp label/border are intentional top-canvas identity cues requested in Justin's sketch. Labels were reduced and lowered after full-resolution crop inspection so they do not cover markers or error bars.",
            "The colored frames align to the title axis; the plot images are inset inside those frames to preserve y-axis label buffer.",
        ],
        "checks_expected": ["post_render_slide_audit", "codex_os_guard after-slide-render", "slide_worker"],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(png_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Render a cleaner slide-4 candidate for the THE-32 low-calo veto flow."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[4]
OUTDIR = ROOT / "outputs/manual-20260612-working-point-slide4-cleanup"
PNG = OUTDIR / "event_energy_veto_threshold_flow_slide_clean_v9.png"
LAYOUT = OUTDIR / "layout_nodes.json"
MANIFEST = OUTDIR / "event_energy_veto_threshold_flow_slide_clean_v9_manifest.json"
SCRIPT = OUTDIR / "event_energy_veto_threshold_flow_slide_clean_v9_script.md"

W, H = 2560, 1440
M = 74

FONT = Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf")
BOLD = Path("/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf")

INK = "#172033"
MUTED = "#526173"
BLUE = "#1F77B4"
BLUE_SOFT = "#EAF4FB"
GOLD = "#A66300"
GOLD_EDGE = "#D9A441"
GREEN = "#1F7A4D"
LINE = "#D8E2EC"
CARD_EDGE = "#BBC7D3"
CARD_FILL = "#FFFFFF"
PANEL_FILL = "#F8FAFC"


def font(px: int, *, bold: bool = False) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(BOLD if bold else FONT), px)


def text_bbox(draw: ImageDraw.ImageDraw, xy: tuple[float, float], text: str, fnt: ImageFont.FreeTypeFont) -> list[float]:
    return [float(v) for v in draw.textbbox(xy, text, font=fnt, anchor="lt")]


def multiline_bbox(draw: ImageDraw.ImageDraw, xy: tuple[float, float], lines: list[str], fnt: ImageFont.FreeTypeFont, line_gap: int) -> list[float]:
    x, y = xy
    boxes = []
    cursor = y
    for line in lines:
        boxes.append(text_bbox(draw, (x, cursor), line, fnt))
        cursor += fnt.size + line_gap
    return [min(b[0] for b in boxes), min(b[1] for b in boxes), max(b[2] for b in boxes), max(b[3] for b in boxes)]


def draw_multiline(
    draw: ImageDraw.ImageDraw,
    xy: tuple[float, float],
    lines: list[str],
    fnt: ImageFont.FreeTypeFont,
    *,
    fill: str,
    line_gap: int = 8,
    anchor: str = "lt",
) -> list[float]:
    x, y = xy
    boxes = []
    cursor = y
    for line in lines:
        draw.text((x, cursor), line, font=fnt, fill=fill, anchor=anchor)
        boxes.append(text_bbox(draw, (x, cursor), line, fnt))
        cursor += fnt.size + line_gap
    return [min(b[0] for b in boxes), min(b[1] for b in boxes), max(b[2] for b in boxes), max(b[3] for b in boxes)]


def rounded(draw: ImageDraw.ImageDraw, bbox: tuple[int, int, int, int], *, fill: str, outline: str, width: int = 3, radius: int = 20) -> None:
    draw.rounded_rectangle(bbox, radius=radius, fill=fill, outline=outline, width=width)


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], *, fill: str = "#93A4B8") -> None:
    x0, y0 = start
    x1, y1 = end
    draw.line((x0, y0, x1, y1), fill=fill, width=5)
    draw.polygon([(x1, y1), (x1 - 22, y1 - 12), (x1 - 22, y1 + 12)], fill=fill)


def draw_runs(
    draw: ImageDraw.ImageDraw,
    xy: tuple[float, float],
    runs: list[tuple[str, ImageFont.FreeTypeFont, int]],
    *,
    fill: str,
) -> list[float]:
    x, y = xy
    cursor = x
    boxes = []
    for text, fnt, dy in runs:
        draw.text((cursor, y + dy), text, font=fnt, fill=fill, anchor="lt")
        bbox = text_bbox(draw, (cursor, y + dy), text, fnt)
        boxes.append(bbox)
        cursor = bbox[2]
    return [min(b[0] for b in boxes), min(b[1] for b in boxes), max(b[2] for b in boxes), max(b[3] for b in boxes)]


def draw_centered_text(
    draw: ImageDraw.ImageDraw,
    center: tuple[float, float],
    text: str,
    fnt: ImageFont.FreeTypeFont,
    *,
    fill: str,
) -> list[float]:
    raw = text_bbox(draw, (0, 0), text, fnt)
    w = raw[2] - raw[0]
    h = raw[3] - raw[1]
    x = center[0] - w / 2 - raw[0]
    y = center[1] - h / 2 - raw[1]
    draw.text((x, y), text, font=fnt, fill=fill, anchor="lt")
    return text_bbox(draw, (x, y), text, fnt)


def measure_runs(runs: list[tuple[str, ImageFont.FreeTypeFont, int]]) -> tuple[float, float]:
    scratch = Image.new("RGB", (2800, 300), "white")
    scratch_draw = ImageDraw.Draw(scratch)
    bbox = draw_runs(scratch_draw, (0, 0), runs, fill=INK)
    return bbox[2] - bbox[0], bbox[3] - bbox[1]


def node(name: str, kind: str, bbox: list[float] | tuple[int, int, int, int], **extra: Any) -> dict[str, Any]:
    out: dict[str, Any] = {"name": name, "kind": kind, "bbox": [round(float(v), 3) for v in bbox]}
    out.update(extra)
    return out


def render() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    im = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(im)
    nodes: list[dict[str, Any]] = []

    title_font = font(78, bold=True)
    formula_font = font(54, bold=True)
    formula_sub = font(34, bold=True)
    formula_small = font(38)
    formula_small_sub = font(25)
    card_title = font(42, bold=True)
    card_body = font(41)
    bottom_font = font(44)

    title_text = "Low-calo veto uses one conservative threshold per centrality bin"
    title_bbox = text_bbox(draw, (M, M), title_text, title_font)
    draw.text((M, M), title_text, font=title_font, fill=INK, anchor="lt")
    nodes.append(
        node(
            "slide title",
            "text",
            title_bbox,
            role="title",
            font_px=title_font.size,
            text=title_text,
            title_anchor=True,
        )
    )

    formula_box = (M, 218, W - M, 480)
    rounded(draw, formula_box, fill="#FFFDF7", outline=GOLD_EDGE, width=4, radius=22)
    nodes.append(node("formula decision band", "box", formula_box, title_axis_align="left", fill_color="#FFFDF7", edge_color=GOLD_EDGE))

    label = "Event removed when y falls below the bin threshold"
    label_bbox_raw = text_bbox(draw, (0, 0), label, formula_small)
    label_xy = ((W - (label_bbox_raw[2] - label_bbox_raw[0])) / 2, 257)
    label_bbox = text_bbox(draw, label_xy, label, formula_small)
    draw.text(label_xy, label, font=formula_small, fill=GOLD, anchor="lt")
    nodes.append(
        node(
            "formula decision label",
            "text",
            label_bbox,
            role="audience",
            font_px=formula_small.size,
            text=label,
            parent="formula decision band",
            intentional_top_aligned=True,
        )
    )

    formula_runs = [
        ("T", formula_font, 0),
        ("bin", formula_sub, 28),
        (" = max[ median(y) - 5 × 1.4826 MAD(y), q", formula_font, 0),
        ("0.1%", formula_sub, 28),
        ("(y) ]", formula_font, 0),
    ]
    formula_w, _ = measure_runs(formula_runs)
    formula_bbox = draw_runs(
        draw,
        ((W - formula_w) / 2, 340),
        formula_runs,
        fill=INK,
    )
    nodes.append(
        node(
            "formula text",
            "text",
            formula_bbox,
            role="audience",
            font_px=formula_font.size,
            text="T_bin = max[median(y) - 5 × 1.4826 MAD(y), q_0.1%(y)]",
            parent="formula decision band",
            intentional_top_aligned=True,
        )
    )
    y_runs = [
        ("y = log", formula_small, 0),
        ("10", formula_small_sub, 19),
        ("(E", formula_small, 0),
        ("CEMC", formula_small_sub, 19),
        (" + E", formula_small, 0),
        ("IHCal", formula_small_sub, 19),
        (" + E", formula_small, 0),
        ("OHCal", formula_small_sub, 19),
        (" + 1)", formula_small, 0),
    ]
    y_w, _ = measure_runs(y_runs)
    y_bbox = draw_runs(
        draw,
        ((W - y_w) / 2, 426),
        y_runs,
        fill=MUTED,
    )
    nodes.append(
        node(
            "energy variable definition",
            "text",
            y_bbox,
            role="audience",
            font_px=formula_small.size,
            text="y = log10(E_CEMC + E_IHCal + E_OHCal + 1)",
            parent="formula decision band",
            intentional_top_aligned=True,
        )
    )

    flow_y0, flow_y1 = 555, 954
    card_gap = 37
    card_w = int((W - 2 * M - 3 * card_gap) / 4)
    card_h = flow_y1 - flow_y0
    cards = [
        ("1", "Input", ["Compute y in one", "5% centrality bin"], CARD_FILL, CARD_EDGE),
        ("2", "Center", ["median(y) marks", "the normal band"], CARD_FILL, CARD_EDGE),
        ("3", "Width", ["MAD = median", "absolute deviation", "sets sigma-like width"], CARD_FILL, CARD_EDGE),
        ("4", "Threshold", ["Go 5 widths below", "center, then apply", "the 0.1% floor"], CARD_FILL, CARD_EDGE),
    ]
    card_boxes: list[tuple[int, int, int, int]] = []
    for idx, (num, label_text, body, fill, edge) in enumerate(cards):
        x0 = M + idx * (card_w + card_gap)
        box = (x0, flow_y0, x0 + card_w, flow_y1)
        card_boxes.append(box)
        rounded(draw, box, fill=fill, outline=edge, width=3, radius=18)
        nodes.append(node(f"flow card {num}", "box", box, symmetry_group="threshold flow cards", fill_color=fill, edge_color=edge))

        badge = (x0 + 32, flow_y0 + 30, x0 + 94, flow_y0 + 92)
        draw.ellipse(badge, fill=BLUE, outline=BLUE)
        nodes.append(node(f"flow card {num} badge circle", "box", badge, parent=f"flow card {num}", fill_color=BLUE, edge_color=BLUE))
        badge_f = font(38, bold=True)
        badge_center = ((badge[0] + badge[2]) / 2, (badge[1] + badge[3]) / 2)
        num_bbox = draw_centered_text(draw, badge_center, num, badge_f, fill="white")
        nodes.append(node(f"flow card {num} badge number", "text", num_bbox, role="process_marker", font_px=badge_f.size, text=num, parent=f"flow card {num} badge circle"))

        title_bbox = text_bbox(draw, (x0 + 122, flow_y0 + 38), label_text, card_title)
        draw.text((x0 + 122, flow_y0 + 38), label_text, font=card_title, fill=INK, anchor="lt")
        nodes.append(node(f"flow card {num} title", "text", title_bbox, role="audience", font_px=card_title.size, text=label_text, parent=f"flow card {num}", intentional_top_aligned=True))

        body_bbox = multiline_bbox(draw, (0, 0), body, card_body, 9)
        body_h = body_bbox[3] - body_bbox[1]
        body_w = body_bbox[2] - body_bbox[0]
        body_x = x0 + (card_w - body_w) / 2
        body_y = flow_y0 + 176 + (140 - body_h) / 2
        actual = draw_multiline(draw, (body_x, body_y), body, card_body, fill=INK, line_gap=9)
        nodes.append(
            node(
                f"flow card {num} body",
                "text",
                actual,
                role="audience",
                font_px=card_body.size,
                text=" ".join(body),
                parent=f"flow card {num}",
                vertical_alignment="top",
                vertical_fill_min_ratio=0.19,
                vertical_fill_max_ratio=0.42,
                **(
                    {
                        "requires_acronym_expansion": "MAD",
                        "acronym_expansion_terms": ["median", "absolute", "deviation"],
                        "functional_role_terms": ["width", "sigma-like"],
                    }
                    if num == "3"
                    else {}
                ),
            )
        )

    for left, right in zip(card_boxes, card_boxes[1:]):
        arrow(draw, (left[2] + 8, (left[1] + left[3]) // 2), (right[0] - 12, (right[1] + right[3]) // 2))

    bottom_box = (M, 1055, W - M, H - M)
    rounded(draw, bottom_box, fill=PANEL_FILL, outline="#BFD0DE", width=3, radius=22)
    nodes.append(node("bottom takeaway card", "box", bottom_box, title_axis_align="left", fill_color=PANEL_FILL, edge_color="#BFD0DE"))

    bullets = [
        "Conservative by construction, with the line far below the normal event-energy band.",
        "Blind to labels and BDTs, because only total event-calo energy enters the veto.",
    ]
    bullet_bbox = multiline_bbox(draw, (0, 0), bullets, bottom_font, 28)
    bullet_w = bullet_bbox[2] - bullet_bbox[0]
    bullet_h = bullet_bbox[3] - bullet_bbox[1]
    bullet_x = M + 124
    bullet_y = bottom_box[1] + ((bottom_box[3] - bottom_box[1]) - bullet_h) / 2
    bullet_actual = draw_multiline(draw, (bullet_x, bullet_y), bullets, bottom_font, fill=INK, line_gap=28)
    for y in (bullet_y + 19, bullet_y + bottom_font.size + 28 + 19):
        x = M + 76
        draw.polygon([(x, y - 12), (x, y + 12), (x + 23, y)], fill=BLUE)
    nodes.append(
        node(
            "bottom takeaway text",
            "text",
            bullet_actual,
            role="audience",
            font_px=bottom_font.size,
            text="\n".join(bullets),
            parent="bottom takeaway card",
            vertical_fill_min_ratio=0.34,
            vertical_fill_max_ratio=0.58,
        )
    )

    im.save(PNG)

    LAYOUT.write_text(
        json.dumps(
            {
                "slide_size": [W, H],
                "title_axis_x": M,
                "title_axis_tolerance_px": 3,
                "minimum_audience_font_px": 34,
                "minimum_title_font_px": 66,
                "vertical_margin_balance": {
                    "top_node": "slide title",
                    "bottom_node": "bottom takeaway card",
                    "target_gap_px": M,
                    "tolerance_px": 4,
                },
                "nodes": nodes,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    SCRIPT.write_text(
        "\n".join(
            [
                "# Low-Calo Veto Threshold Slide Script",
                "",
                "This slide explains how the event-energy veto threshold is set before BDT training.",
                "The variable is total event calorimeter energy, compressed as y = log10(E_CEMC + E_IHCal + E_OHCal + 1).",
                "For each 5% centrality bin, the normal band is summarized by the median and a robust MAD width.",
                "The threshold is set five sigma-like widths below the median, with the 0.1% quantile acting as a floor.",
                "The intended takeaway is that this veto is conservative and independent of truth labels, BDT score, and sample weights.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    MANIFEST.write_text(
        json.dumps(
            {
                "png": str(PNG),
                "layout_nodes": str(LAYOUT),
                "speaker_script": str(SCRIPT),
                "purpose": "Clean local PNG candidate for working-point deck slide 4; no Google Slides mutation.",
                "source_slide": "167x-He2rOOBO2i4nNS6Pdcqu7Wv03GeFMuWH9tRRx-8 / g3e8f2d411e1_0_17",
                "formula": "T_bin = max[median(y) - 5*1.4826*MAD(y), q_0.1%(y)]",
                "cut_variable": "y = log10(E_CEMC + E_IHCal + E_OHCal + 1)",
                "design_notes": [
                    "No subtitles; one context line below the title.",
                    "Crisp white cards with minimal tinting.",
                    "Four-sided outer margin balance enforced by post-render audit.",
                    "Details such as the Gaussian MAD conversion stay in speaker script, not on the canvas.",
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    print(PNG)


if __name__ == "__main__":
    render()

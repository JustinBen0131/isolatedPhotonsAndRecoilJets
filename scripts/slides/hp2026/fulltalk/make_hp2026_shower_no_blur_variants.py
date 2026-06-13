#!/usr/bin/env python3
"""Local no-blur strategy prototypes for the HP2026 shower-shape build.

This is a local design comparison only. It imports the current HP2026
full-talk drawing primitives and writes PNG previews; it does not mutate the
Google Slides deck or the primary full-talk generator.
"""

from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw

import make_hp2026_fulltalk_candidates as full


ROOT = full.ROOT
OUT_DIR = ROOT / "outputs/manual-20260608-slides8-10-no-blur-variants"
W, H = full.W, full.H
EXPAND_PANEL_TOP = 306
EXPAND_PANEL_BOTTOM = 1282
SYNTH_FIGURE_LABEL_Y = 132
SYNTH_GRID_Y = 176
SYNTH_GRID_CELL = 58
SYNTH_CAPTION_Y = 490
SYNTH_LOGIC_Y0 = 612
SYNTH_LOGIC_Y1 = 762
SYNTH_TAKEAWAY_Y0 = 802
SYNTH_TAKEAWAY_Y1 = 936

PANELS = [
    {
        "short": "Core compactness",
        "title": "Is the core compact?",
        "next": "next check",
        "takeaway": "Core energy stays concentrated near the seed.",
        "accent": full.PHOTON_DARK,
        "draw": full.draw_compact_core_panel,
    },
    {
        "short": "Shoulder width",
        "title": "Are the shoulders narrow?",
        "next": "next check",
        "takeaway": "Seed-excluded shoulders stay narrow.",
        "accent": full.SPHENIX_BLUE,
        "draw": full.draw_narrow_shoulders_panel,
    },
    {
        "short": "Stretch / split",
        "title": "Is the shower stretched or split?",
        "next": "final check",
        "takeaway": "Energy remains confined rather than elongated or split.",
        "accent": full.TEAL,
        "draw": full.draw_elongation_panel,
    },
]

SUBTITLES = [
    "Start with one local question: is the energy concentrated in the core?",
    "Keep the core check as context, then test whether the surrounding energy stays narrow.",
    "Together, these checks separate compact photon-like deposits from broader or split backgrounds.",
]

BOTTOM_STEP_TEXT = [
    "Core compactness asks whether the candidate starts as one concentrated photon-like deposit.",
    "Adding shoulder width checks whether nearby energy stays narrow after the seed tower is removed.",
    "Together, these variables ask whether the EMCal energy looks like one compact photon, rather than a broad, elongated, or split decay-like shower.",
]


def draw_concept_badge(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    index: int,
    *,
    size: int = 66,
    offset: tuple[int, int] = (36, 34),
    bar: bool = True,
) -> None:
    """Draw a consistent numbered concept badge for the shower-shape sequence."""
    x0, y0, _, y1 = box
    accent = PANELS[index]["accent"]
    if bar:
        draw.rounded_rectangle((x0, y0, x0 + 12, y1), radius=6, fill=(*accent, 220))
    cx0 = x0 + offset[0]
    cy0 = y0 + offset[1]
    chip = (cx0, cy0, cx0 + size, cy0 + size)
    draw.ellipse(chip, fill=(255, 255, 255, 255), outline=(*accent, 235), width=4)
    num = str(index + 1)
    nf = full.font(full.TIMES_BOLD, int(size * 0.55))
    bbox = draw.textbbox((0, 0), num, font=nf)
    cx = (chip[0] + chip[2]) / 2
    cy = (chip[1] + chip[3]) / 2
    draw.text((cx - (bbox[0] + bbox[2]) / 2, cy - (bbox[1] + bbox[3]) / 2), num, font=nf, fill=accent)


def draw_core_fraction_eq(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    size: int,
    *,
    prefix: str = "Core fraction = ",
    fill: tuple[int, int, int] = full.INK,
    prefix_bold: bool = True,
) -> int:
    small = max(15, round(size * 0.58))
    sub = max(14, round(size * 0.64))
    parts: list[tuple[str, int, int, Path]] = []
    if prefix:
        parts.append((prefix, size, 0, full.TIMES_BOLD if prefix_bold else full.TIMES))
    parts.extend(
        [
            ("E", size, 0, full.TIMES_BOLD),
            ("2x2", small, round(size * 0.34), full.TIMES_BOLD),
            (" / ", size, 0, full.TIMES),
            ("E", size, 0, full.TIMES_BOLD),
            ("cluster", small, round(size * 0.34), full.TIMES_BOLD),
        ]
    )
    return full.draw_formula_run(draw, xy, parts, fill=fill)


def draw_e1x1_e3x3_eq(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    size: int,
    *,
    prefix: str = "",
    fill: tuple[int, int, int] = full.INK,
) -> int:
    sub = max(14, round(size * 0.62))
    parts: list[tuple[str, int, int, Path]] = []
    if prefix:
        parts.append((prefix, size, 0, full.TIMES_BOLD))
    parts.extend(
        [
            ("E", size, 0, full.TIMES_BOLD),
            ("1x1", sub, round(size * 0.34), full.TIMES_BOLD),
            (" / E", size, 0, full.TIMES_BOLD),
            ("3x3", sub, round(size * 0.34), full.TIMES_BOLD),
        ]
    )
    return full.draw_formula_run(draw, xy, parts, fill=fill)


def draw_e3x2_e3x5_eq(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    size: int,
    *,
    fill: tuple[int, int, int] = full.INK,
) -> int:
    sub = max(14, round(size * 0.62))
    return full.draw_formula_run(
        draw,
        xy,
        [
            ("E", size, 0, full.TIMES_BOLD),
            ("3x2", sub, round(size * 0.34), full.TIMES_BOLD),
            (" / E", size, 0, full.TIMES_BOLD),
            ("3x5", sub, round(size * 0.34), full.TIMES_BOLD),
            (" = narrow-strip share", size, 0, full.TIMES),
        ],
        fill=fill,
    )


def panel_boxes() -> list[tuple[int, int, int, int]]:
    panel_top = 324
    panel_bottom = 1138
    gap = 38
    panel_w = (W - 2 * 132 - 2 * gap) // 3
    return [
        (132, panel_top, 132 + panel_w, panel_bottom),
        (132 + panel_w + gap, panel_top, 132 + 2 * panel_w + gap, panel_bottom),
        (132 + 2 * (panel_w + gap), panel_top, W - 132, panel_bottom),
    ]


def new_base(frame: int, subtitle: str | None = None) -> Image.Image:
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*full.SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*full.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*full.PHOTON, 255))
    draw.text((132, 76), "Reading photon-like shower shapes", font=full.font(full.TIMES_BOLD, 86), fill=full.INK)
    draw.text((136, 190), subtitle or SUBTITLES[frame], font=full.font(full.TIMES_ITALIC, 56), fill=full.MUTED)
    draw.line((132, 286, W - 132, 286), fill=(221, 226, 232, 255), width=3)
    full.add_top_right_sphenix_logo_like_slide2(img)
    return img


def bottom_takeaway(img: Image.Image, text: str) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    box = (132, 1192, W - 132, 1308)
    draw.rounded_rectangle(
        box,
        radius=8,
        fill=(239, 246, 250, 255),
        outline=(213, 226, 235, 255),
        width=2,
    )
    draw_wrapped_vcenter(draw, text, (176, box[1] + 14, W - 176, box[3] - 14), full.font(full.TIMES_ITALIC, 34), fill=full.BLUE, line_gap=6)


def wrap_lines(draw: ImageDraw.ImageDraw, text: str, max_width: int, fnt) -> list[str]:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        trial = word if not current else f"{current} {word}"
        if full.text_box(draw, trial, fnt)[0] <= max_width:
            current = trial
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    return lines


def draw_wrapped_vcenter(
    draw: ImageDraw.ImageDraw,
    text: str,
    box: tuple[int, int, int, int],
    fnt,
    *,
    fill: tuple[int, int, int] = full.BLUE,
    line_gap: int = 8,
) -> None:
    x0, y0, x1, y1 = box
    lines = wrap_lines(draw, text, x1 - x0, fnt)
    if not lines:
        return
    heights = [full.text_box(draw, line, fnt)[1] for line in lines]
    total_h = sum(heights) + line_gap * (len(lines) - 1)
    y = y0 + ((y1 - y0) - total_h) / 2
    for line, h in zip(lines, heights):
        draw.text((x0, y), line, font=fnt, fill=fill)
        y += h + line_gap


TAKEAWAY_BOX_STYLES = [
    {
        "fill": (255, 249, 234, 255),
        "outline": (230, 195, 132, 255),
        "label": full.PHOTON_DARK,
        "body": full.BLUE,
    },
    {
        "fill": (238, 248, 254, 255),
        "outline": (173, 215, 241, 255),
        "label": full.SPHENIX_BLUE,
        "body": full.BLUE,
    },
    {
        "fill": (238, 249, 250, 255),
        "outline": (166, 212, 220, 255),
        "label": full.TEAL,
        "body": full.BLUE,
    },
]


def wrap_words_by_width(draw: ImageDraw.ImageDraw, words: list[str], max_width: int, fnt) -> list[str]:
    lines: list[str] = []
    current: list[str] = []
    for word in words:
        trial = " ".join([*current, word])
        if not current or full.text_box(draw, trial, fnt)[0] <= max_width:
            current.append(word)
        else:
            lines.append(" ".join(current))
            current = [word]
    if current:
        lines.append(" ".join(current))
    return lines


def draw_photon_like_takeaway(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    index: int,
    body: str,
    *,
    font_size: int = 36,
    pad_x: int = 32,
    line_gap: int = 8,
) -> None:
    """Draw a centered takeaway with only the lead label bolded."""
    style = TAKEAWAY_BOX_STYLES[index]
    draw.rounded_rectangle(box, radius=10, fill=style["fill"], outline=style["outline"], width=2)
    label = "Photon-like:"
    label_font = full.font(full.TIMES_BOLD, font_size)
    body_font = full.font(full.TIMES, font_size)
    label_w, label_h = full.text_box(draw, label, label_font)
    space_w = full.text_box(draw, " ", body_font)[0]
    x0, y0, x1, y1 = box
    left = x0 + pad_x
    right = x1 - pad_x
    max_width = right - left
    body_words = body.split()

    first_line_words: list[str] = []
    remaining = max_width - label_w - space_w
    while body_words:
        trial = " ".join([*first_line_words, body_words[0]])
        if first_line_words and full.text_box(draw, trial, body_font)[0] > remaining:
            break
        if not first_line_words and full.text_box(draw, body_words[0], body_font)[0] > remaining:
            break
        first_line_words.append(body_words.pop(0))

    rest_lines = wrap_words_by_width(draw, body_words, max_width, body_font)
    while first_line_words and len(rest_lines) == 1:
        first_body = " ".join(first_line_words)
        first_w = label_w + space_w + full.text_box(draw, first_body, body_font)[0]
        second_w = full.text_box(draw, rest_lines[0], body_font)[0]
        if second_w >= 0.62 * first_w:
            break
        moved = first_line_words[-1]
        candidate_second = f"{moved} {rest_lines[0]}"
        candidate_first = " ".join(first_line_words[:-1])
        candidate_first_w = label_w + (space_w + full.text_box(draw, candidate_first, body_font)[0] if candidate_first else 0)
        candidate_second_w = full.text_box(draw, candidate_second, body_font)[0]
        if candidate_first and candidate_first_w <= max_width and candidate_second_w <= max_width:
            first_line_words.pop()
            rest_lines[0] = candidate_second
        else:
            break
    lines: list[tuple[bool, str]] = [(True, " ".join(first_line_words))]
    lines.extend((False, line) for line in rest_lines)
    line_heights = [max(label_h if is_first else 0, full.text_box(draw, text, body_font)[1]) for is_first, text in lines]
    total_h = sum(line_heights) + line_gap * (len(lines) - 1)
    y = y0 + ((y1 - y0) - total_h) / 2
    for (is_first, text), line_h in zip(lines, line_heights):
        if is_first:
            body_w = full.text_box(draw, text, body_font)[0] if text else 0
            line_w = label_w + (space_w + body_w if text else 0)
            x = x0 + ((x1 - x0) - line_w) / 2
            draw.text((x, y), label, font=label_font, fill=style["label"])
            if text:
                draw.text((x + label_w + space_w, y), text, font=body_font, fill=style["body"])
        else:
            line_w = full.text_box(draw, text, body_font)[0]
            x = x0 + ((x1 - x0) - line_w) / 2
            draw.text((x, y), text, font=body_font, fill=style["body"])
        y += line_h + line_gap


def draw_centered_text(
    draw: ImageDraw.ImageDraw,
    text: str,
    center: tuple[float, float],
    fnt,
    fill: tuple[int, int, int],
) -> None:
    tw, th = full.text_box(draw, text, fnt)
    draw.text((center[0] - tw / 2, center[1] - th / 2), text, font=fnt, fill=fill)


def draw_arrowhead(
    draw: ImageDraw.ImageDraw,
    tip: tuple[float, float],
    direction: tuple[float, float],
    fill: tuple[int, int, int],
    *,
    size: int = 12,
) -> None:
    dx, dy = direction
    norm = math.hypot(dx, dy) or 1.0
    ux, uy = dx / norm, dy / norm
    px, py = -uy, ux
    base = (tip[0] - ux * size, tip[1] - uy * size)
    points = [
        tip,
        (base[0] + px * size * 0.48, base[1] + py * size * 0.48),
        (base[0] - px * size * 0.48, base[1] - py * size * 0.48),
    ]
    draw.polygon(points, fill=fill)


def draw_double_arrow(
    draw: ImageDraw.ImageDraw,
    start: tuple[float, float],
    end: tuple[float, float],
    fill: tuple[int, int, int],
    *,
    width: int = 4,
    head: int = 13,
) -> None:
    draw.line((start[0], start[1], end[0], end[1]), fill=(*fill, 230), width=width)
    draw_arrowhead(draw, end, (end[0] - start[0], end[1] - start[1]), fill, size=head)
    draw_arrowhead(draw, start, (start[0] - end[0], start[1] - end[1]), fill, size=head)


def draw_width_symbol(
    draw: ImageDraw.ImageDraw,
    xy: tuple[float, float],
    subscript: str,
    size: int,
    fill: tuple[int, int, int],
) -> int:
    """Draw w with eta/phi as a true subscript and return the next x position."""
    base_font = full.font(full.TIMES_BOLD, size)
    sub_font = full.font(full.TIMES_BOLD, max(16, round(size * 0.64)))
    x, y = xy
    base_w, _ = full.text_box(draw, "w", base_font)
    sub_offset = round(size * 0.36)
    draw.text((x, y), "w", font=base_font, fill=fill)
    draw.text((x + base_w + 1, y + sub_offset), subscript, font=sub_font, fill=fill)
    return round(x + base_w + full.text_box(draw, subscript, sub_font)[0] + 3)


def width_symbol_box(
    draw: ImageDraw.ImageDraw,
    subscript: str,
    size: int,
) -> tuple[int, int]:
    base_font = full.font(full.TIMES_BOLD, size)
    sub_font = full.font(full.TIMES_BOLD, max(16, round(size * 0.64)))
    base_w, base_h = full.text_box(draw, "w", base_font)
    sub_w, sub_h = full.text_box(draw, subscript, sub_font)
    return base_w + sub_w + 3, max(base_h, round(size * 0.36) + sub_h)


def draw_centered_width_symbol(
    draw: ImageDraw.ImageDraw,
    subscript: str,
    center: tuple[float, float],
    size: int,
    fill: tuple[int, int, int],
) -> None:
    w, h = width_symbol_box(draw, subscript, size)
    draw_width_symbol(draw, (center[0] - w / 2, center[1] - h / 2), subscript, size, fill)


def draw_width_explanation(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    *,
    size: int,
    suffix: str,
) -> None:
    x, y = xy
    x = draw_width_symbol(draw, (x, y), "η", size, full.SPHENIX_BLUE)
    fnt = full.font(full.TIMES, size)
    gap_text = " and "
    draw.text((x, y), gap_text, font=fnt, fill=full.MUTED)
    x += full.text_box(draw, gap_text, fnt)[0]
    x = draw_width_symbol(draw, (x, y), "φ", size, (230, 70, 45))
    draw.text((x + 2, y), suffix, font=fnt, fill=full.MUTED)


def draw_core_geometry_labels(draw: ImageDraw.ImageDraw, gx: int, gy: int, cell: int, *, large: bool = False) -> None:
    label_size = 26 if large else 20
    seed_size = 30 if large else 20
    line_w = 3 if large else 2
    two_label = (gx + cell - (42 if large else 10), gy + cell - (48 if large else 20))
    three_label = (gx + 3 * cell - (16 if large else 10), gy + 3 * cell + (16 if large else 8))
    draw.line(
        (gx + cell * 1.72, gy + cell * 1.22, gx + cell * 1.18, gy + cell * 1.02),
        fill=(*full.PHOTON_DARK, 205),
        width=line_w,
    )
    draw.line(
        (gx + cell * 3.62, gy + cell * 3.56, gx + cell * 4.0, gy + cell * 4.0),
        fill=(*full.SPHENIX_BLUE, 205),
        width=line_w,
    )
    full.callout_label(draw, two_label, "2x2 core", full.PHOTON_DARK, size=label_size)
    full.callout_label(draw, three_label, "3x3 local core", full.SPHENIX_BLUE, size=label_size)
    seed_text = "seed"
    seed_font = full.font(full.TIMES_BOLD, seed_size)
    stw, sth = full.text_box(draw, seed_text, seed_font)
    sx = gx + 2 * cell + (8 if large else 10)
    sy = gy + 2 * cell + (0.82 * cell if large else 0.62 * cell)
    if large:
        draw.rounded_rectangle((sx - 8, sy - 4, sx + stw + 8, sy + sth + 4), radius=6, fill=(255, 255, 255, 230))
    draw.text((sx, sy), seed_text, font=seed_font, fill=full.INK)


def draw_width_labels(
    draw: ImageDraw.ImageDraw,
    grid: tuple[int, int, int, int],
    *,
    compact: bool,
    label_size: int,
    eta_label_inside: bool = False,
) -> None:
    gx0, gy0, gx1, gy1 = grid
    width_color = full.SPHENIX_BLUE
    bracket_color = (230, 70, 45)
    x = gx0 + 18 if eta_label_inside else gx0 - 22
    draw_double_arrow(draw, (x, gy0 + 58), (x, gy1 - 58), width_color, width=4, head=12)
    eta_center_x = x + 32 if eta_label_inside else x - 26
    draw_centered_width_symbol(
        draw,
        "η",
        (eta_center_x, (gy0 + gy1) / 2),
        label_size,
        width_color,
    )
    draw_centered_width_symbol(
        draw,
        "φ",
        ((gx0 + gx1) / 2, gy0 - (18 if label_size >= 34 else 13)),
        label_size,
        bracket_color,
    )


def draw_region_tags_for_strip_ratio(
    draw: ImageDraw.ImageDraw,
    grid: tuple[int, int, int, int],
    cell: int,
    *,
    label_size: int = 20,
) -> None:
    gx0, gy0, gx1, _ = grid
    full.callout_label(
        draw,
        (gx0 + 0.84 * cell, gy0 + 1.10 * cell),
        "3x5 region",
        full.SPHENIX_BLUE,
        size=label_size,
    )
    full.callout_label(
        draw,
        (gx0 + 1.34 * cell, gy0 + 2.28 * cell),
        "3x2 strip",
        full.PHOTON_DARK,
        size=label_size,
    )
    draw.line(
        (gx0 + 1.58 * cell, gy0 + 1.62 * cell, gx0 + 1.05 * cell, gy0 + 0.14 * cell),
        fill=(*full.SPHENIX_BLUE, 190),
        width=2,
    )
    draw.line(
        (gx0 + 2.14 * cell, gy0 + 2.6 * cell, gx1 - 1.2 * cell, gy0 + 2.5 * cell),
        fill=(*full.PHOTON_DARK, 205),
        width=2,
    )


def draw_stretch_region_labels(
    draw: ImageDraw.ImageDraw,
    left_grid: tuple[int, int, int, int],
    right_grid: tuple[int, int, int, int],
    cell: int,
    *,
    size: int = 24,
) -> None:
    """Label the 3x2 and 3x5 regions beside the maps with short leaders."""
    lgx0, lgy0, _, _ = left_grid
    _, rgy0, rgx1, _ = right_grid
    orange = full.PHOTON_DARK
    blue = full.SPHENIX_BLUE
    orange_fill = (255, 246, 222, 255)
    blue_fill = (231, 246, 255, 255)
    fnt = full.font(full.TIMES_BOLD, size)

    def chip(xy: tuple[float, float], text: str, color: tuple[int, int, int], fill: tuple[int, int, int]) -> tuple[int, int, int, int]:
        tw, th = full.text_box(draw, text, fnt)
        rect = (int(xy[0]), int(xy[1]), int(xy[0] + tw + 28), int(xy[1] + th + 14))
        draw.rounded_rectangle(rect, radius=8, fill=fill, outline=(*color, 255), width=3)
        draw_centered_text(draw, text, ((rect[0] + rect[2]) / 2, (rect[1] + rect[3]) / 2), fnt, color)
        return rect

    left_rect = chip((lgx0 + 6, lgy0 + 2 * cell + 6), "3x2 strip", orange, orange_fill)
    blue_w = full.text_box(draw, "3x5 local region", fnt)[0] + 28
    right_rect = chip((rgx1 - blue_w - 8, rgy0 + 8), "3x5 local region", blue, blue_fill)

    draw.line((left_rect[2], (left_rect[1] + left_rect[3]) / 2, lgx0 + cell, lgy0 + 2.5 * cell), fill=(*orange, 220), width=3)
    draw.line((right_rect[2] - 2, right_rect[3], rgx1 - 0.95 * cell, rgy0 + 0.25 * cell), fill=(*blue, 220), width=3)


def draw_strip_ratio_key(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    *,
    y: int,
    size: int = 23,
) -> None:
    x0, _, x1, _ = box
    items = [
        ("3x2 strip", full.PHOTON_DARK, (255, 246, 222, 255)),
        ("3x5 local region", full.SPHENIX_BLUE, (231, 246, 255, 255)),
    ]
    fnt = full.font(full.TIMES_BOLD, size)
    widths = [full.text_box(draw, text, fnt)[0] + 50 for text, _, _ in items]
    gap = 26
    total = sum(widths) + gap
    x = x0 + ((x1 - x0) - total) / 2
    for (text, color, fill), width in zip(items, widths):
        pill = (int(x), y, int(x + width), y + 46)
        draw.rounded_rectangle(pill, radius=9, fill=fill, outline=(*color, 255), width=3)
        draw_centered_text(draw, text, ((pill[0] + pill[2]) / 2, (pill[1] + pill[3]) / 2), fnt, color)
        x += width + gap


def draw_core_region_key(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    *,
    y: int,
    size: int = 22,
) -> None:
    x0, _, x1, _ = box
    items = [
        ("2x2 core", full.PHOTON_DARK, (255, 246, 222, 255)),
        ("3x3 local core", full.SPHENIX_BLUE, (231, 246, 255, 255)),
    ]
    fnt = full.font(full.TIMES_BOLD, size)
    widths = [full.text_box(draw, text, fnt)[0] + 48 for text, _, _ in items]
    gap = 24
    total = sum(widths) + gap
    x = x0 + ((x1 - x0) - total) / 2
    for (text, color, fill), width in zip(items, widths):
        pill = (int(x), y, int(x + width), y + 46)
        draw.rounded_rectangle(pill, radius=9, fill=fill, outline=(*color, 255), width=3)
        draw_centered_text(draw, text, ((pill[0] + pill[2]) / 2, (pill[1] + pill[3]) / 2), fnt, color)
        x += width + gap


def draw_labeled_takeaway_vcenter(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    label: str,
    text: str,
    label_font,
    text_font,
    *,
    fill: tuple[int, int, int] = full.BLUE,
    line_gap: int = 8,
    pad_x: int = 38,
) -> None:
    x0, y0, x1, y1 = box
    lines = wrap_lines(draw, text, x1 - x0 - 2 * pad_x, text_font)
    label_h = full.text_box(draw, label, label_font)[1]
    body_h = sum(full.text_box(draw, line, text_font)[1] for line in lines) + line_gap * max(0, len(lines) - 1)
    total_h = label_h + 12 + body_h
    y = y0 + ((y1 - y0) - total_h) / 2
    draw.text((x0 + pad_x, y), label, font=label_font, fill=fill)
    y += label_h + 12
    for line in lines:
        h = full.text_box(draw, line, text_font)[1]
        draw.text((x0 + pad_x, y), line, font=text_font, fill=fill)
        y += h + line_gap


def draw_placeholder(
    img: Image.Image,
    box: tuple[int, int, int, int],
    index: int,
    label: str = "next check",
    dim_completed: bool = False,
) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    accent = PANELS[index]["accent"]
    fill = (248, 251, 253, 255) if not dim_completed else (252, 253, 254, 255)
    outline = (213, 226, 235, 255)
    draw.rounded_rectangle(box, radius=12, fill=fill, outline=outline, width=2)
    draw.rounded_rectangle((x0, y0, x0 + 12, y1), radius=6, fill=(*accent, 220))

    chip = (x0 + 42, y0 + 38, x0 + 106, y0 + 102)
    draw.ellipse(chip, fill=(255, 255, 255, 255), outline=(*accent, 235), width=3)
    num = str(index + 1)
    nf = full.font(full.TIMES_BOLD, 38)
    tw, th = full.text_box(draw, num, nf)
    draw.text(
        (chip[0] + (chip[2] - chip[0] - tw) / 2, chip[1] + (chip[3] - chip[1] - th) / 2 - 1),
        num,
        font=nf,
        fill=accent,
    )

    title_font = full.font(full.TIMES_BOLD, 35)
    muted = (81, 91, 105)
    draw.text((x0 + 132, y0 + 45), f"Step {index + 1}", font=title_font, fill=muted)
    body = "Revealed on next click" if label == "next check" else label
    full.draw_wrapped(
        draw,
        body,
        (x0 + 132, y0 + 96),
        x1 - x0 - 190,
        full.font(full.TIMES_ITALIC, 27),
        fill=full.LIGHT_MUTED,
        line_gap=4,
    )

    # A quiet glyph gives the future card structure without leaking content.
    cx, cy = (x0 + x1) // 2, y0 + 392
    draw.rounded_rectangle(
        (cx - 154, cy - 114, cx + 154, cy + 114),
        radius=18,
        fill=(255, 255, 255, 210),
        outline=(222, 230, 238, 255),
        width=2,
    )
    for dx in (-72, 0, 72):
        alpha = 70 if not dim_completed else 52
        draw.ellipse((cx + dx - 28, cy - 28, cx + dx + 28, cy + 28), outline=(*accent, alpha), width=6)
    draw.line((cx - 100, cy + 88, cx + 100, cy + 88), fill=(215, 225, 235, 230), width=3)


def draw_completed_summary(img: Image.Image, box: tuple[int, int, int, int], index: int) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    accent = PANELS[index]["accent"]
    draw.rounded_rectangle(box, radius=12, fill=(251, 253, 254, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((x0, y0, x0 + 12, y1), radius=6, fill=(*accent, 210))

    chip = (x0 + 40, y0 + 36, x0 + 102, y0 + 98)
    draw.ellipse(chip, fill=(255, 255, 255, 255), outline=(*accent, 225), width=3)
    draw.line((chip[0] + 17, chip[1] + 33, chip[0] + 29, chip[1] + 45), fill=(*accent, 255), width=5)
    draw.line((chip[0] + 29, chip[1] + 45, chip[0] + 47, chip[1] + 20), fill=(*accent, 255), width=5)
    draw.text((x0 + 128, y0 + 43), PANELS[index]["short"], font=full.font(full.TIMES_BOLD, 35), fill=full.INK)
    full.draw_wrapped(
        draw,
        PANELS[index]["takeaway"],
        (x0 + 128, y0 + 98),
        x1 - x0 - 172,
        full.font(full.TIMES_ITALIC, 28),
        fill=full.BLUE,
        line_gap=5,
    )
    # Keep the rest intentionally calm; this is a completed memory trace.
    mid_y = y0 + 415
    draw.rounded_rectangle((x0 + 70, mid_y - 54, x1 - 70, mid_y + 54), radius=12, fill=(240, 247, 251, 255), outline=(218, 226, 235, 255), width=2)
    draw.text((x0 + 102, mid_y - 18), "covered", font=full.font(full.TIMES_ITALIC, 29), fill=full.LIGHT_MUTED)


def draw_progress_rail(draw: ImageDraw.ImageDraw, frame: int, y: int = 302) -> None:
    labels = ["Core", "Shoulders", "Stretch / split"]
    xs = [760, 1280, 1800]
    for i in range(2):
        draw.line((xs[i] + 48, y, xs[i + 1] - 48, y), fill=(178, 194, 210, 255), width=5)
        draw.polygon([(xs[i + 1] - 48, y), (xs[i + 1] - 68, y - 12), (xs[i + 1] - 68, y + 12)], fill=(178, 194, 210, 255))
    for i, (x, label) in enumerate(zip(xs, labels)):
        active = i == frame
        done = i < frame
        accent = PANELS[i]["accent"]
        fill = accent if active else ((255, 255, 255) if done else (243, 248, 251))
        outline = accent if (active or done) else full.PANEL_EDGE
        draw.ellipse((x - 42, y - 42, x + 42, y + 42), fill=(*fill, 255), outline=(*outline, 255), width=4)
        if done:
            draw.line((x - 16, y + 2, x - 3, y + 16), fill=(*accent, 255), width=5)
            draw.line((x - 3, y + 16, x + 20, y - 18), fill=(*accent, 255), width=5)
        else:
            text = str(i + 1)
            tf = full.font(full.TIMES_BOLD, 36)
            tw, th = full.text_box(draw, text, tf)
            draw.text((x - tw / 2, y - th / 2 - 2), text, font=tf, fill=(255, 255, 255) if active else accent)
        rail_label = label if i <= frame else f"Step {i + 1}"
        lf = full.font(full.TIMES_BOLD if active else full.TIMES, 25)
        tw, _ = full.text_box(draw, rail_label, lf)
        draw.text((x - tw / 2, y + 54), rail_label, font=lf, fill=full.INK if active else full.MUTED)


def variant_a(frame: int) -> Image.Image:
    img = new_base(frame)
    boxes = panel_boxes()
    revealed = set(range(frame + 1))
    for idx, box in enumerate(boxes):
        if idx in revealed:
            PANELS[idx]["draw"](img, box)
        else:
            draw_placeholder(img, box, idx)
    bottom_takeaway(img, BOTTOM_STEP_TEXT[frame])
    full.draw_hp2026_identity_footer(img)
    return img


def variant_b(frame: int) -> Image.Image:
    img = new_base(frame)
    boxes = panel_boxes()
    for idx, box in enumerate(boxes):
        if idx < frame:
            draw_completed_summary(img, box, idx)
        elif idx == frame:
            PANELS[idx]["draw"](img, box)
        else:
            draw_placeholder(img, box, idx)
    bottom_takeaway(
        img,
        [
            "Current check: core compactness. Upcoming checks stay informationally absent until the click.",
            "Current check: seed-excluded shoulder width. The core check remains as a concise memory trace.",
            "Current check: elongation or splitting. The three checks now support the full photon-like-shower takeaway.",
        ][frame],
    )
    full.draw_hp2026_identity_footer(img)
    return img


def variant_c(frame: int) -> Image.Image:
    img = new_base(frame, subtitle=SUBTITLES[frame])
    draw = ImageDraw.Draw(img, "RGBA")
    draw_progress_rail(draw, frame, y=294)
    # One active card uses the middle/lower canvas so the audience never reads ahead.
    active_box = (348, 396, W - 348, 1138)
    PANELS[frame]["draw"](img, active_box)
    bottom_takeaway(
        img,
        [
            "Step 1 isolates the local-core question before any later variables appear.",
            "Step 2 adds the seed-excluded shoulder-width question without showing the final check yet.",
            "Step 3 completes the sequence: compact core, narrow shoulders, and no elongation or splitting.",
        ][frame],
    )
    full.draw_hp2026_identity_footer(img)
    return img


def draw_simple_opaque_cover(img: Image.Image, box: tuple[int, int, int, int], index: int) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    accent = PANELS[index]["accent"]
    draw.rounded_rectangle(
        box,
        radius=12,
        fill=(248, 251, 253, 255),
        outline=(213, 226, 235, 255),
        width=2,
    )
    draw.rounded_rectangle((x0, y0, x0 + 12, y1), radius=6, fill=(*accent, 215))
    chip = (x0 + 40, y0 + 34, x0 + 96, y0 + 90)
    draw.ellipse(chip, fill=(255, 255, 255, 255), outline=(*accent, 220), width=3)
    step = str(index + 1)
    sf = full.font(full.TIMES_BOLD, 32)
    tw, th = full.text_box(draw, step, sf)
    draw.text(
        (chip[0] + (chip[2] - chip[0] - tw) / 2, chip[1] + (chip[3] - chip[1] - th) / 2 - 1),
        step,
        font=sf,
        fill=accent,
    )


def variant_d_simple_opaque(frame: int) -> Image.Image:
    # This intentionally preserves the current full slide language and layout;
    # only the inactive-panel treatment changes from blur/veil to opaque cover.
    img = full.build_slide05_shower_shape_variables(set())
    boxes = panel_boxes()
    hidden = {1, 2} if frame == 0 else ({2} if frame == 1 else set())
    for idx in sorted(hidden):
        draw_simple_opaque_cover(img, boxes[idx], idx)
    return img.convert("RGB")


def draw_core_summary_card(img: Image.Image, box: tuple[int, int, int, int], *, show_takeaway: bool = True) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw_concept_badge(draw, box, 0, size=60, offset=(34, 34), bar=True)
    title = "Is the core compact?"
    tf = full.font(full.TIMES_BOLD, 44)
    tw, _ = full.text_box(draw, title, tf)
    draw.text((x0 + (x1 - x0 - tw) / 2 + 18, y0 + 38), title, font=tf, fill=full.INK)

    values = [
        [0.03, 0.05, 0.07, 0.05, 0.03],
        [0.05, 0.18, 0.32, 0.18, 0.05],
        [0.07, 0.36, 0.95, 0.42, 0.07],
        [0.05, 0.20, 0.38, 0.19, 0.05],
        [0.03, 0.05, 0.07, 0.05, 0.03],
    ]
    diffuse = [
        [0.12, 0.20, 0.28, 0.20, 0.12],
        [0.22, 0.38, 0.52, 0.40, 0.22],
        [0.30, 0.55, 0.78, 0.58, 0.32],
        [0.22, 0.40, 0.54, 0.42, 0.22],
        [0.12, 0.20, 0.30, 0.20, 0.12],
    ]
    cell = SYNTH_GRID_CELL
    gap = 46
    total = 2 * 5 * cell + gap
    left_x = x0 + (x1 - x0 - total) // 2
    gy = y0 + SYNTH_GRID_Y
    label = "Local EMCal energy map"
    lf = full.font(full.TIMES_BOLD, 37)
    ltw, _ = full.text_box(draw, label, lf)
    draw.text((x0 + (x1 - x0 - ltw) / 2 + 18, y0 + SYNTH_FIGURE_LABEL_Y), label, font=lf, fill=full.BLUE)
    grids = [
        full.draw_small_tower_grid(
            draw,
            (left_x, gy),
            cell,
            values,
            highlight=(1, 1, 3, 3),
            highlight_color=full.PHOTON_DARK,
            secondary=(1, 1, 4, 4),
            secondary_color=full.SPHENIX_BLUE,
        ),
        full.draw_small_tower_grid(
            draw,
            (left_x + 5 * cell + gap, gy),
            cell,
            diffuse,
            highlight=(1, 1, 3, 3),
            highlight_color=full.PHOTON_DARK,
            secondary=(1, 1, 4, 4),
            secondary_color=full.SPHENIX_BLUE,
        ),
    ]
    for grid in grids:
        gx0, gy0, _, _ = grid
        draw.rectangle((gx0 + cell, gy0 + cell, gx0 + 4 * cell, gy0 + 4 * cell), outline=(*full.SPHENIX_BLUE, 255), width=5)
        draw.rectangle((gx0 + cell, gy0 + cell, gx0 + 3 * cell, gy0 + 3 * cell), outline=(*full.PHOTON_DARK, 255), width=5)
        draw.ellipse((gx0 + 2.5 * cell - 8, gy0 + 2.5 * cell - 8, gx0 + 2.5 * cell + 8, gy0 + 2.5 * cell + 8), fill=(0, 0, 0, 255))
    for caption, grid in (("photon-like", grids[0]), ("background-like", grids[1])):
        gx0, _, gx1, _ = grid
        cf = full.font(full.TIMES_BOLD, 34)
        ctw, _ = full.text_box(draw, caption, cf)
        draw.text((gx0 + (gx1 - gx0 - ctw) / 2, y0 + SYNTH_CAPTION_Y), caption, font=cf, fill=full.BLUE if caption == "photon-like" else full.MUTED)
    draw_core_region_key(draw, (x0 + 54, y0, x1 - 54, y0 + 1), y=y0 + 548, size=25)

    logic = (x0 + 56, y0 + SYNTH_LOGIC_Y0, x1 - 56, y0 + SYNTH_LOGIC_Y1)
    draw.rounded_rectangle(logic, radius=10, fill=(249, 251, 253, 255), outline=(219, 227, 236, 255), width=2)
    draw_core_fraction_eq(draw, (logic[0] + 28, logic[1] + 30), 32, prefix="Core fraction = ")
    draw.line((logic[0] + 28, logic[1] + 78, logic[2] - 28, logic[1] + 78), fill=(216, 225, 235, 255), width=2)
    draw_e1x1_e3x3_eq(draw, (logic[0] + 28, logic[1] + 104), 32, prefix="")
    draw.text((logic[0] + 196, logic[1] + 104), " = center-tower dominance", font=full.font(full.TIMES, 32), fill=full.INK)
    if show_takeaway:
        summary = (x0 + 42, y0 + SYNTH_TAKEAWAY_Y0, x1 - 42, y0 + SYNTH_TAKEAWAY_Y1)
        draw_photon_like_takeaway(
            draw,
            summary,
            0,
            "energy stays concentrated in the core.",
            font_size=39,
            line_gap=7,
        )


def draw_large_core_focus_panel(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw_concept_badge(draw, box, 0, size=72, offset=(76, 48), bar=True)
    title = "Is the core compact?"
    tf = full.font(full.TIMES_BOLD, 64)
    tw, _ = full.text_box(draw, title, tf)
    draw.text((x0 + (x1 - x0 - tw) / 2 + 24, y0 + 24), title, font=tf, fill=full.INK)

    values = [
        [0.03, 0.05, 0.07, 0.05, 0.03],
        [0.05, 0.18, 0.32, 0.18, 0.05],
        [0.07, 0.36, 0.95, 0.42, 0.07],
        [0.05, 0.20, 0.38, 0.19, 0.05],
        [0.03, 0.05, 0.07, 0.05, 0.03],
    ]
    diffuse = [
        [0.12, 0.20, 0.28, 0.20, 0.12],
        [0.22, 0.38, 0.52, 0.40, 0.22],
        [0.30, 0.55, 0.78, 0.58, 0.32],
        [0.22, 0.40, 0.54, 0.42, 0.22],
        [0.12, 0.20, 0.30, 0.20, 0.12],
    ]
    cell = 68
    gap = 92
    total_grid_w = 2 * 5 * cell + gap
    gx = x0 + 80 + ((1012 - 80) - total_grid_w) // 2
    gy = y0 + 170
    label = "Local EMCal energy map"
    lf = full.font(full.TIMES_BOLD, 46)
    ltw, _ = full.text_box(draw, label, lf)
    label_x = x0 + 80 + (932 - ltw) / 2
    label_y = y0 + 108
    draw.text((label_x, label_y), label, font=lf, fill=full.BLUE)
    grids = [
        full.draw_small_tower_grid(
            draw,
            (gx, gy),
            cell,
            values,
            highlight=(1, 1, 3, 3),
            highlight_color=full.PHOTON_DARK,
            secondary=(1, 1, 4, 4),
            secondary_color=full.SPHENIX_BLUE,
        ),
        full.draw_small_tower_grid(
            draw,
            (gx + 5 * cell + gap, gy),
            cell,
            diffuse,
            highlight=(1, 1, 3, 3),
            highlight_color=full.PHOTON_DARK,
            secondary=(1, 1, 4, 4),
            secondary_color=full.SPHENIX_BLUE,
        ),
    ]
    for grid in grids:
        gx0, gy0, _, _ = grid
        draw.rectangle((gx0 + cell, gy0 + cell, gx0 + 4 * cell, gy0 + 4 * cell), outline=(*full.SPHENIX_BLUE, 255), width=6)
        draw.rectangle((gx0 + cell, gy0 + cell, gx0 + 3 * cell, gy0 + 3 * cell), outline=(*full.PHOTON_DARK, 255), width=6)
        draw.ellipse((gx0 + 2.5 * cell - 10, gy0 + 2.5 * cell - 10, gx0 + 2.5 * cell + 10, gy0 + 2.5 * cell + 10), fill=(0, 0, 0, 255))
    for caption, grid in (("photon-like", grids[0]), ("background-like", grids[1])):
        gx0, _, gx1, _ = grid
        cf = full.font(full.TIMES_BOLD, 43)
        ctw, _ = full.text_box(draw, caption, cf)
        draw.text((gx0 + (gx1 - gx0 - ctw) / 2, y0 + 528), caption, font=cf, fill=full.BLUE if caption == "photon-like" else full.MUTED)
    draw_core_region_key(draw, (x0 + 112, y0, x0 + 980, y0 + 1), y=y0 + 578, size=28)

    logic = (x0 + 80, y0 + 650, x0 + 1012, y0 + 802)
    draw.rounded_rectangle(logic, radius=10, fill=(249, 251, 253, 255), outline=(219, 227, 236, 255), width=2)
    draw_core_fraction_eq(draw, (logic[0] + 30, logic[1] + 30), 38, prefix="Core fraction = ")
    draw.line((logic[0] + 30, logic[1] + 78, logic[2] - 30, logic[1] + 78), fill=(216, 225, 235, 255), width=2)
    draw_e1x1_e3x3_eq(draw, (logic[0] + 30, logic[1] + 104), 38)
    draw.text((logic[0] + 226, logic[1] + 104), " = center-tower dominance in the local core", font=full.font(full.TIMES, 38), fill=full.INK)

    summary = (x0 + 80, y0 + 834, x0 + 1012, y0 + 966)
    draw_photon_like_takeaway(
        draw,
        summary,
        0,
        "energy stays concentrated in the core.",
        font_size=40,
        line_gap=8,
    )

    note = (x0 + 1124, y0 + 168, x1 - 76, y0 + 712)
    draw.rounded_rectangle(note, radius=12, fill=(248, 251, 254, 255), outline=(218, 226, 235, 255), width=2)
    draw.text((note[0] + 38, note[1] + 28), "What this first check establishes", font=full.font(full.TIMES_BOLD, 50), fill=full.INK)
    draw.line((note[0] + 38, note[1] + 92, note[2] - 38, note[1] + 92), fill=(214, 224, 234, 255), width=2)
    rows = [
        ("local core", "energy concentrated near the cluster seed"),
        ("core fraction", "2x2 core energy divided by total cluster energy"),
        ("E1x1/E3x3", "center-tower dominance inside the local 3x3"),
    ]
    y = note[1] + 128
    for label, body in rows:
        row_box = (note[0] + 36, y - 20, note[2] - 36, y + 82)
        draw.rounded_rectangle(row_box, radius=8, fill=(255, 255, 255, 245), outline=(226, 233, 241, 255), width=1)
        if label == "E1x1/E3x3":
            draw_e1x1_e3x3_eq(draw, (row_box[0] + 24, y + 2), 36, fill=full.PHOTON_DARK)
        else:
            draw.text((row_box[0] + 24, y), label, font=full.font(full.TIMES_BOLD, 37), fill=full.PHOTON_DARK)
        full.draw_wrapped(draw, body, (row_box[0] + 296, y + 1), row_box[2] - row_box[0] - 324, full.font(full.TIMES, 35), fill=full.INK, line_gap=6)
        y += 118
    takeaway = (note[0], note[3] + 22, note[2], note[3] + 216)
    draw.rounded_rectangle(takeaway, radius=10, fill=(242, 248, 253, 255), outline=(197, 220, 239, 255), width=2)
    takeaway_text = "A prompt-photon-like cluster should begin as one concentrated EMCal deposit."
    tf = full.font(full.TIMES_BOLD, 46)
    lines = wrap_lines(draw, takeaway_text, takeaway[2] - takeaway[0] - 92, tf)
    line_gap = 9
    heights = [full.text_box(draw, line, tf)[1] for line in lines]
    total_h = sum(heights) + line_gap * max(0, len(lines) - 1)
    ty = takeaway[1] + ((takeaway[3] - takeaway[1]) - total_h) / 2
    for line, lh in zip(lines, heights):
        tw, _ = full.text_box(draw, line, tf)
        draw.text((takeaway[0] + ((takeaway[2] - takeaway[0]) - tw) / 2, ty), line, font=tf, fill=full.BLUE)
        ty += lh + line_gap


def draw_large_shoulders_panel(img: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw_concept_badge(draw, box, 1, size=66, offset=(40, 36), bar=True)
    title = "Are the shoulders narrow?"
    tf = full.font(full.TIMES_BOLD, 60)
    tw, _ = full.text_box(draw, title, tf)
    draw.text((x0 + (x1 - x0 - tw) / 2, y0 + 30), title, font=tf, fill=full.INK)

    compact = [
        [0.01, 0.02, 0.04, 0.02, 0.01],
        [0.02, 0.08, 0.20, 0.08, 0.02],
        [0.04, 0.18, 0.92, 0.18, 0.04],
        [0.02, 0.08, 0.20, 0.08, 0.02],
        [0.01, 0.02, 0.04, 0.02, 0.01],
    ]
    broad = [
        [0.10, 0.18, 0.25, 0.18, 0.10],
        [0.18, 0.34, 0.46, 0.34, 0.18],
        [0.25, 0.48, 0.80, 0.48, 0.25],
        [0.18, 0.34, 0.46, 0.34, 0.18],
        [0.10, 0.18, 0.25, 0.18, 0.10],
    ]
    cell = 78
    total_grid_w = 2 * 5 * cell + 172
    left_x = x0 + (x1 - x0 - total_grid_w) // 2
    gy = y0 + 190
    label = "Energy map with center tower omitted"
    lf = full.font(full.TIMES_BOLD, 46)
    ltw, _ = full.text_box(draw, label, lf)
    draw.text((x0 + (x1 - x0 - ltw) / 2, y0 + 108), label, font=lf, fill=full.BLUE)
    grids = [
        full.draw_small_tower_grid(draw, (left_x, gy), cell, compact),
        full.draw_small_tower_grid(draw, (left_x + 5 * cell + 172, gy), cell, broad),
    ]
    for idx, grid in enumerate(grids):
        gx0, gy0, gx1, _ = grid
        if idx == 0:
            start, end = gx0 + 138, gx1 - 138
        else:
            start, end = gx0 + 34, gx1 - 34
        draw_double_arrow(draw, (start, gy0 + 12), (end, gy0 + 12), (230, 70, 45), width=5, head=14)
        center_cell = (gx0 + 2 * cell, gy0 + 2 * cell, gx0 + 3 * cell, gy0 + 3 * cell)
        draw.rectangle(center_cell, fill=(255, 255, 255, 215), outline=(160, 168, 176, 215), width=2)
        draw.line((center_cell[0] + 10, center_cell[1] + 10, center_cell[2] - 10, center_cell[3] - 10), fill=(160, 168, 176, 190), width=3)
        draw.line((center_cell[0] + 10, center_cell[3] - 10, center_cell[2] - 10, center_cell[1] + 10), fill=(160, 168, 176, 190), width=3)
        draw.ellipse((gx0 + 2.5 * cell - 10, gy0 + 2.5 * cell - 10, gx0 + 2.5 * cell + 10, gy0 + 2.5 * cell + 10), fill=(0, 0, 0, 255))
        draw_width_labels(draw, grid, compact=idx == 0, label_size=42)
    for label, grid in (("photon-like", grids[0]), ("background-like", grids[1])):
        gx0, _, gx1, _ = grid
        lf = full.font(full.TIMES_BOLD, 43)
        tw, _ = full.text_box(draw, label, lf)
        draw.text((gx0 + (gx1 - gx0 - tw) / 2, y0 + 612), label, font=lf, fill=full.BLUE if label == "photon-like" else full.MUTED)
    logic = (x0 + 86, y0 + 674, x1 - 86, y0 + 800)
    draw.rounded_rectangle(logic, radius=10, fill=(249, 251, 253, 255), outline=(219, 227, 236, 255), width=2)
    full.draw_formula_run(
        draw,
        (logic[0] + 34, logic[1] + 28),
        [("Shoulder width", 42, 0, full.TIMES_BOLD), (" = seed-excluded spread around the core", 42, 0, full.TIMES)],
        fill=full.INK,
    )
    draw_width_explanation(draw, (logic[0] + 34, logic[1] + 82), size=38, suffix=" describe the surrounding energy")
    takeaway = (x0 + 86, y0 + 824, x1 - 86, y0 + 948)
    draw_photon_like_takeaway(
        draw,
        takeaway,
        1,
        "narrow shoulders around the core.",
        font_size=43,
        line_gap=8,
    )


def draw_shoulders_summary_card(img: Image.Image, box: tuple[int, int, int, int], *, show_takeaway: bool = True) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw_concept_badge(draw, box, 1, size=58, offset=(30, 30), bar=True)
    title = "Are the shoulders narrow?"
    tf = full.font(full.TIMES_BOLD, 38)
    tw, _ = full.text_box(draw, title, tf)
    draw.text((x0 + (x1 - x0 - tw) / 2 + 18, y0 + 40), title, font=tf, fill=full.INK)

    compact = [
        [0.01, 0.02, 0.04, 0.02, 0.01],
        [0.02, 0.08, 0.20, 0.08, 0.02],
        [0.04, 0.18, 0.92, 0.18, 0.04],
        [0.02, 0.08, 0.20, 0.08, 0.02],
        [0.01, 0.02, 0.04, 0.02, 0.01],
    ]
    broad = [
        [0.10, 0.18, 0.25, 0.18, 0.10],
        [0.18, 0.34, 0.46, 0.34, 0.18],
        [0.25, 0.48, 0.80, 0.48, 0.25],
        [0.18, 0.34, 0.46, 0.34, 0.18],
        [0.10, 0.18, 0.25, 0.18, 0.10],
    ]
    cell = SYNTH_GRID_CELL
    gap = 46
    total = 2 * 5 * cell + gap
    left_x = x0 + (x1 - x0 - total) // 2
    gy = y0 + SYNTH_GRID_Y
    grids = [
        full.draw_small_tower_grid(draw, (left_x, gy), cell, compact),
        full.draw_small_tower_grid(draw, (left_x + 5 * cell + gap, gy), cell, broad),
    ]
    for idx, grid in enumerate(grids):
        gx0, gy0, gx1, _ = grid
        start, end = (gx0 + 106, gx1 - 106) if idx == 0 else (gx0 + 26, gx1 - 26)
        draw_double_arrow(draw, (start, gy0 + 10), (end, gy0 + 10), (230, 70, 45), width=4, head=11)
        center_cell = (gx0 + 2 * cell, gy0 + 2 * cell, gx0 + 3 * cell, gy0 + 3 * cell)
        draw.rectangle(center_cell, fill=(255, 255, 255, 215), outline=(160, 168, 176, 215), width=2)
        draw.line((center_cell[0] + 8, center_cell[1] + 8, center_cell[2] - 8, center_cell[3] - 8), fill=(160, 168, 176, 190), width=3)
        draw.line((center_cell[0] + 8, center_cell[3] - 8, center_cell[2] - 8, center_cell[1] + 8), fill=(160, 168, 176, 190), width=3)
        draw.ellipse((gx0 + 2.5 * cell - 8, gy0 + 2.5 * cell - 8, gx0 + 2.5 * cell + 8, gy0 + 2.5 * cell + 8), fill=(0, 0, 0, 255))
        draw_width_labels(draw, grid, compact=idx == 0, label_size=31, eta_label_inside=True)
    for label, grid in (("photon-like", grids[0]), ("background-like", grids[1])):
        gx0, _, gx1, _ = grid
        lf = full.font(full.TIMES_BOLD, 34)
        tw, _ = full.text_box(draw, label, lf)
        draw.text((gx0 + (gx1 - gx0 - tw) / 2, y0 + SYNTH_CAPTION_Y), label, font=lf, fill=full.BLUE if label == "photon-like" else full.MUTED)
    logic = (x0 + 54, y0 + SYNTH_LOGIC_Y0, x1 - 54, y0 + SYNTH_LOGIC_Y1)
    draw.rounded_rectangle(logic, radius=10, fill=(249, 251, 253, 255), outline=(219, 227, 236, 255), width=2)
    full.draw_formula_run(
        draw,
        (logic[0] + 28, logic[1] + 28),
        [("Shoulder width", 33, 0, full.TIMES_BOLD), (" = seed-excluded spread", 33, 0, full.TIMES)],
        fill=full.INK,
    )
    draw_width_explanation(draw, (logic[0] + 28, logic[1] + 74), size=32, suffix=" describe surrounding energy")
    if show_takeaway:
        summary = (x0 + 42, y0 + SYNTH_TAKEAWAY_Y0, x1 - 42, y0 + SYNTH_TAKEAWAY_Y1)
        draw_photon_like_takeaway(
            draw,
            summary,
            1,
            "narrow shoulders around the core.",
            font_size=38,
            line_gap=7,
        )


def draw_stretch_summary_card(img: Image.Image, box: tuple[int, int, int, int], *, show_takeaway: bool = True) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw_concept_badge(draw, box, 2, size=58, offset=(30, 30), bar=True)
    tf = full.font(full.TIMES_BOLD, 35)
    for line, yoff in (("Is the shower", 34), ("stretched or split?", 76)):
        tw, _ = full.text_box(draw, line, tf)
        draw.text((x0 + (x1 - x0 - tw) / 2 + 22, y0 + yoff), line, font=tf, fill=full.INK)

    compact = [
        [0.02, 0.04, 0.06, 0.04, 0.02],
        [0.03, 0.10, 0.20, 0.10, 0.03],
        [0.05, 0.32, 0.95, 0.34, 0.05],
        [0.03, 0.10, 0.20, 0.10, 0.03],
        [0.02, 0.04, 0.06, 0.04, 0.02],
    ]
    split = [
        [0.02, 0.06, 0.10, 0.08, 0.04],
        [0.06, 0.18, 0.28, 0.24, 0.10],
        [0.18, 0.62, 0.98, 0.78, 0.44],
        [0.16, 0.54, 0.72, 0.58, 0.34],
        [0.04, 0.12, 0.18, 0.14, 0.06],
    ]
    cell = SYNTH_GRID_CELL
    gap = 46
    total = 2 * 5 * cell + gap
    left_x = x0 + (x1 - x0 - total) // 2
    gy = y0 + SYNTH_GRID_Y
    grids = [
        full.draw_small_tower_grid(draw, (left_x, gy), cell, compact),
        full.draw_small_tower_grid(draw, (left_x + 5 * cell + gap, gy), cell, split),
    ]
    for grid in grids:
        gx0, gy0, _, _ = grid
        draw.rectangle((gx0 + cell, gy0, gx0 + 4 * cell, gy0 + 5 * cell), outline=(*full.SPHENIX_BLUE, 235), width=5)
        draw.rectangle((gx0 + cell, gy0 + 2 * cell, gx0 + 4 * cell, gy0 + 4 * cell), outline=(*full.PHOTON_DARK, 245), width=5)
        draw.ellipse((gx0 + 2.5 * cell - 8, gy0 + 2.5 * cell - 8, gx0 + 2.5 * cell + 8, gy0 + 2.5 * cell + 8), fill=(0, 0, 0, 255))
    plot_label = "Narrow strip compared with wider local region"
    plf = full.font(full.TIMES_BOLD, 32)
    pltw, _ = full.text_box(draw, plot_label, plf)
    draw.text((x0 + (x1 - x0 - pltw) / 2 + 14, y0 + SYNTH_FIGURE_LABEL_Y), plot_label, font=plf, fill=full.BLUE)
    for label, grid in (("photon-like", grids[0]), ("background-like", grids[1])):
        gx0, _, gx1, _ = grid
        lf = full.font(full.TIMES_BOLD, 34)
        tw, _ = full.text_box(draw, label, lf)
        draw.text((gx0 + (gx1 - gx0 - tw) / 2, y0 + SYNTH_CAPTION_Y), label, font=lf, fill=full.BLUE if label == "photon-like" else full.MUTED)
    draw_strip_ratio_key(draw, (x0 + 54, y0, x1 - 54, y0 + 1), y=y0 + 548, size=25)
    logic = (x0 + 54, y0 + SYNTH_LOGIC_Y0, x1 - 54, y0 + SYNTH_LOGIC_Y1)
    draw.rounded_rectangle(logic, radius=10, fill=(249, 251, 253, 255), outline=(219, 227, 236, 255), width=2)
    eq_size = 38
    eq_text_w = full.text_box(draw, "E3x2 / E3x5 = narrow-strip share", full.font(full.TIMES, eq_size))[0]
    draw_e3x2_e3x5_eq(draw, (logic[0] + max(28, ((logic[2] - logic[0]) - eq_text_w) // 2), logic[1] + 34), eq_size)
    if show_takeaway:
        summary = (x0 + 42, y0 + SYNTH_TAKEAWAY_Y0, x1 - 42, y0 + SYNTH_TAKEAWAY_Y1)
        draw_photon_like_takeaway(
            draw,
            summary,
            2,
            "little elongation or splitting in the local region.",
            font_size=36,
            line_gap=7,
        )


def draw_clean_full_synthesis() -> Image.Image:
    img = new_base(
        2,
        subtitle="Together, these checks separate compact photon-like deposits from broader or split backgrounds.",
    )
    panel_top = EXPAND_PANEL_TOP
    panel_bottom = EXPAND_PANEL_BOTTOM
    gap = 38
    panel_w = (W - 2 * 132 - 2 * gap) // 3
    boxes = [
        (132, panel_top, 132 + panel_w, panel_bottom),
        (132 + panel_w + gap, panel_top, 132 + 2 * panel_w + gap, panel_bottom),
        (132 + 2 * (panel_w + gap), panel_top, W - 132, panel_bottom),
    ]
    draw_core_summary_card(img, boxes[0], show_takeaway=True)
    draw_shoulders_summary_card(img, boxes[1], show_takeaway=True)
    draw_stretch_summary_card(img, boxes[2], show_takeaway=True)
    full.draw_hp2026_identity_footer(img)
    return img.convert("RGB")


def variant_e_expand_collapse(frame: int) -> Image.Image:
    if frame == 2:
        return draw_clean_full_synthesis()

    subtitle = [
        "Start with one local question: is the energy concentrated in the core?",
        "Keep the core check as context, then test whether the surrounding energy stays narrow.",
    ][frame]
    img = new_base(frame, subtitle=subtitle)
    if frame == 0:
        draw_large_core_focus_panel(img, (132, EXPAND_PANEL_TOP, W - 132, EXPAND_PANEL_BOTTOM))
    else:
        draw_core_summary_card(img, (132, EXPAND_PANEL_TOP, 885, EXPAND_PANEL_BOTTOM))
        draw_large_shoulders_panel(img, (935, EXPAND_PANEL_TOP, W - 132, EXPAND_PANEL_BOTTOM))
    full.draw_hp2026_identity_footer(img)
    return img.convert("RGB")


VARIANTS = [
    ("A", "Opaque placeholders", "Minimal replacement: same three-card geometry; future cards are present but informationally absent.", variant_a),
    ("B", "Current / completed / upcoming", "Most guided delivery: current card stays detailed; covered cards collapse to one-line memory traces.", variant_b),
    ("C", "Single active card + progress rail", "Cleanest one-idea-at-a-time reveal; strongest redesign and least future-reading risk.", variant_c),
]


def save_variant_frames() -> dict[str, list[Path]]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    paths: dict[str, list[Path]] = {}
    for key, _, _, builder in VARIANTS:
        paths[key] = []
        for frame in range(3):
            img = builder(frame).convert("RGB")
            path = OUT_DIR / f"variant_{key}_frame{frame + 1}_slide{8 + frame:02d}.png"
            img.save(path, "PNG")
            paths[key].append(path)
    return paths


def make_contact_sheet(paths: dict[str, list[Path]]) -> Path:
    thumb_w, thumb_h = 768, 432
    left_label_w = 410
    top_label_h = 88
    row_h = thumb_h + 88
    sheet = Image.new("RGB", (left_label_w + 3 * thumb_w, top_label_h + len(VARIANTS) * row_h), (246, 249, 252))
    draw = ImageDraw.Draw(sheet, "RGBA")
    title_font = full.font(full.TIMES_BOLD, 36)
    sub_font = full.font(full.TIMES, 24)
    small_font = full.font(full.TIMES_ITALIC, 21)
    col_font = full.font(full.TIMES_BOLD, 28)

    draw.text((26, 22), "Slides 8-10 no-blur strategy comparison", font=title_font, fill=full.INK)
    for col, label in enumerate(("Slide 8 / first click", "Slide 9 / second click", "Slide 10 / full reveal")):
        x = left_label_w + col * thumb_w
        tw, _ = full.text_box(draw, label, col_font)
        draw.text((x + (thumb_w - tw) / 2, 30), label, font=col_font, fill=full.BLUE)

    for row, (key, title, desc, _) in enumerate(VARIANTS):
        y = top_label_h + row * row_h
        draw.rounded_rectangle((18, y + 20, left_label_w - 22, y + row_h - 20), radius=14, fill=(255, 255, 255, 255), outline=(218, 226, 235, 255), width=2)
        draw.text((42, y + 52), f"Variant {key}", font=full.font(full.TIMES_BOLD, 34), fill=full.INK)
        draw.text((42, y + 104), title, font=sub_font, fill=full.BLUE)
        full.draw_wrapped(draw, desc, (42, y + 148), left_label_w - 84, small_font, fill=full.MUTED, line_gap=5)
        verdict = {
            "A": "Recommended default",
            "B": "Best guided narration",
            "C": "Cleanest but bigger redesign",
        }[key]
        draw.rounded_rectangle((42, y + row_h - 82, left_label_w - 58, y + row_h - 40), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=1)
        draw.text((62, y + row_h - 75), verdict, font=full.font(full.TIMES_BOLD, 22), fill=full.TEAL if key == "A" else full.BLUE)

        for col, path in enumerate(paths[key]):
            thumb = Image.open(path).convert("RGB").resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
            x = left_label_w + col * thumb_w
            sheet.paste(thumb, (x, y + 48))
            draw.rectangle((x, y + 48, x + thumb_w - 1, y + 48 + thumb_h - 1), outline=(204, 216, 228, 255), width=2)

    contact = OUT_DIR / "slides8_10_no_blur_strategy_comparison.png"
    sheet.save(contact, "PNG")
    return contact


def save_simple_opaque_option() -> tuple[list[Path], Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    paths = []
    for frame in range(3):
        img = variant_d_simple_opaque(frame)
        path = OUT_DIR / f"variant_D_simple_opaque_frame{frame + 1}_slide{8 + frame:02d}.png"
        img.save(path, "PNG")
        paths.append(path)

    thumb_w, thumb_h = 768, 432
    label_h = 82
    left_w = 460
    sheet = Image.new("RGB", (left_w + 3 * thumb_w, label_h + thumb_h + 40), (246, 249, 252))
    draw = ImageDraw.Draw(sheet, "RGBA")
    draw.rounded_rectangle((18, 22, left_w - 24, label_h + thumb_h + 18), radius=14, fill=(255, 255, 255, 255), outline=(218, 226, 235, 255), width=2)
    draw.text((42, 54), "Variant D", font=full.font(full.TIMES_BOLD, 38), fill=full.INK)
    draw.text((42, 112), "Simple opaque cover", font=full.font(full.TIMES, 28), fill=full.BLUE)
    full.draw_wrapped(
        draw,
        "Closest to the current blurred build: keep the same slide language and geometry; cover unrevealed cards with an opaque clean panel.",
        (42, 164),
        left_w - 88,
        full.font(full.TIMES_ITALIC, 24),
        fill=full.MUTED,
        line_gap=5,
    )
    draw.rounded_rectangle((42, label_h + thumb_h - 54, left_w - 58, label_h + thumb_h - 12), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=1)
    draw.text((62, label_h + thumb_h - 47), "Lowest-risk replacement", font=full.font(full.TIMES_BOLD, 22), fill=full.TEAL)
    for col, (path, label) in enumerate(zip(paths, ("Slide 8 / first click", "Slide 9 / second click", "Slide 10 / full reveal"))):
        x = left_w + col * thumb_w
        lf = full.font(full.TIMES_BOLD, 28)
        tw, _ = full.text_box(draw, label, lf)
        draw.text((x + (thumb_w - tw) / 2, 30), label, font=lf, fill=full.BLUE)
        thumb = Image.open(path).convert("RGB").resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        sheet.paste(thumb, (x, label_h))
        draw.rectangle((x, label_h, x + thumb_w - 1, label_h + thumb_h - 1), outline=(204, 216, 228, 255), width=2)
    contact = OUT_DIR / "variant_D_simple_opaque_cover_contact_sheet.png"
    sheet.save(contact, "PNG")
    return paths, contact


def save_expand_collapse_option() -> tuple[list[Path], Path]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    paths = []
    for frame in range(3):
        img = variant_e_expand_collapse(frame)
        path = OUT_DIR / f"variant_E_expand_collapse_frame{frame + 1}_slide{8 + frame:02d}.png"
        img.save(path, "PNG")
        paths.append(path)

    thumb_w, thumb_h = 768, 432
    label_h = 82
    left_w = 470
    sheet = Image.new("RGB", (left_w + 3 * thumb_w, label_h + thumb_h + 40), (246, 249, 252))
    draw = ImageDraw.Draw(sheet, "RGBA")
    draw.rounded_rectangle((18, 22, left_w - 24, label_h + thumb_h + 18), radius=14, fill=(255, 255, 255, 255), outline=(218, 226, 235, 255), width=2)
    draw.text((42, 54), "Variant E", font=full.font(full.TIMES_BOLD, 38), fill=full.INK)
    draw.text((42, 112), "Expand → collapse → synthesize", font=full.font(full.TIMES, 27), fill=full.BLUE)
    full.draw_wrapped(
        draw,
        "Uses space for one large readable active idea, then compresses covered content into context before the final three-card synthesis.",
        (42, 164),
        left_w - 88,
        full.font(full.TIMES_ITALIC, 24),
        fill=full.MUTED,
        line_gap=5,
    )
    draw.rounded_rectangle((42, label_h + thumb_h - 54, left_w - 58, label_h + thumb_h - 12), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=1)
    draw.text((62, label_h + thumb_h - 47), "Best readability-first option", font=full.font(full.TIMES_BOLD, 22), fill=full.TEAL)
    for col, (path, label) in enumerate(zip(paths, ("Slide 8 / large core", "Slide 9 / core summary + shoulders", "Slide 10 / full synthesis"))):
        x = left_w + col * thumb_w
        lf = full.font(full.TIMES_BOLD, 28)
        tw, _ = full.text_box(draw, label, lf)
        draw.text((x + (thumb_w - tw) / 2, 30), label, font=lf, fill=full.BLUE)
        thumb = Image.open(path).convert("RGB").resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        sheet.paste(thumb, (x, label_h))
        draw.rectangle((x, label_h, x + thumb_w - 1, label_h + thumb_h - 1), outline=(204, 216, 228, 255), width=2)
    contact = OUT_DIR / "variant_E_expand_collapse_contact_sheet.png"
    sheet.save(contact, "PNG")
    return paths, contact


def main() -> None:
    paths = save_variant_frames()
    contact = make_contact_sheet(paths)
    simple_paths, simple_contact = save_simple_opaque_option()
    expand_paths, expand_contact = save_expand_collapse_option()
    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "google_slides_mutation": False,
        "source_generator": str((ROOT / "scripts/slides/hp2026/fulltalk/make_hp2026_fulltalk_candidates.py").relative_to(ROOT)),
        "source_functions": [
            "build_slide05_shower_shape_variables()",
            "slide05_focus_core_compactness()",
            "slide05_focus_shoulders()",
            "slide04()",
            "soften_focus_region() is the current blur/fade helper being replaced in these prototypes",
        ],
        "design_goal": "Replace blurred/faded unrevealed shower-shape cards with staged disclosure treatments where future content is structurally present but unreadable.",
        "variants": [
            {
                "key": key,
                "title": title,
                "description": desc,
                "outputs": [str(p.relative_to(ROOT)) for p in paths[key]],
            }
            for key, title, desc, _ in VARIANTS
        ],
        "contact_sheet": str(contact.relative_to(ROOT)),
        "simple_opaque_cover_contact_sheet": str(simple_contact.relative_to(ROOT)),
        "simple_opaque_cover_outputs": [str(p.relative_to(ROOT)) for p in simple_paths],
        "expand_collapse_contact_sheet": str(expand_contact.relative_to(ROOT)),
        "expand_collapse_outputs": [str(p.relative_to(ROOT)) for p in expand_paths],
        "notes": [
            "No Google Slides deck was changed.",
            "The primary full-talk generator was imported, not edited.",
            "These frames are local strategy prototypes for Justin to choose a no-blur pattern before propagating across HPslides_v1.",
        ],
    }
    manifest_path = OUT_DIR / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(contact)
    print(simple_contact)
    print(expand_contact)
    print(manifest_path)
    for key in sorted(paths):
        for path in paths[key]:
            print(path)


if __name__ == "__main__":
    main()

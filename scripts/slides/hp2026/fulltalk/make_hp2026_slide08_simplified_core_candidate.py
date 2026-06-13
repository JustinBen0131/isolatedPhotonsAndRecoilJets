#!/usr/bin/env python3
"""Local PNG-only candidate for simplified HP2026 Slide 8.

This renders a review candidate only. It does not mutate Google Slides.
"""

from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw

import make_hp2026_fulltalk_candidates as full


ROOT = full.ROOT
OUT_DIR = ROOT / "outputs/manual-20260610-slide8-compact-core-simplified"
PNG_PATH = OUT_DIR / "slide08_core_compact_simplified_candidate.png"
SCRIPT_PATH = OUT_DIR / "slide08_core_compact_simplified_candidate_script.md"
MANIFEST_PATH = OUT_DIR / "manifest.json"

W, H = full.W, full.H
PHOTON_LABEL_RED = (196, 45, 39)
GEOM_CORE_PURPLE = (107, 76, 154)
GEOM_TOWER_GREEN = (27, 127, 90)
CENTER_TOWER_RED = PHOTON_LABEL_RED


def draw_centered(draw: ImageDraw.ImageDraw, text: str, center: tuple[float, float], font, fill) -> None:
    bbox = draw.textbbox((0, 0), text, font=font)
    draw.text(
        (center[0] - (bbox[0] + bbox[2]) / 2, center[1] - (bbox[1] + bbox[3]) / 2),
        text,
        font=font,
        fill=fill,
    )


def draw_arrow(
    draw: ImageDraw.ImageDraw,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    fill: tuple[int, int, int],
    width: int = 4,
    head: int = 15,
) -> None:
    draw.line((start[0], start[1], end[0], end[1]), fill=(*fill, 230), width=width)
    angle = math.atan2(end[1] - start[1], end[0] - start[0])
    pts = [
        end,
        (end[0] - head * math.cos(angle - 0.42), end[1] - head * math.sin(angle - 0.42)),
        (end[0] - head * math.cos(angle + 0.42), end[1] - head * math.sin(angle + 0.42)),
    ]
    draw.polygon(pts, fill=(*fill, 235))


def draw_pill_label(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    *,
    color: tuple[int, int, int],
    fill: tuple[int, int, int, int] = (255, 255, 255, 242),
    font_size: int = 32,
) -> tuple[int, int, int, int]:
    x, y = xy
    fnt = full.font(full.TIMES_BOLD, font_size)
    tw, th = full.text_box(draw, text, fnt)
    box = (x, y, x + tw + 34, y + th + 18)
    draw.rounded_rectangle(box, radius=9, fill=fill, outline=(*color, 245), width=3)
    draw.text((x + 17, y + 8), text, font=fnt, fill=color)
    return box


def draw_e11_e33(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    size: int,
    *,
    fill: tuple[int, int, int] = full.INK,
) -> int:
    sub = max(18, round(size * 0.62))
    return full.draw_formula_run(
        draw,
        xy,
        [
            ("E", size, 0, full.TIMES_BOLD),
            ("11", sub, round(size * 0.34), full.TIMES_BOLD),
            (" / E", size, 0, full.TIMES_BOLD),
            ("33", sub, round(size * 0.34), full.TIMES_BOLD),
        ],
        fill=fill,
    )


def formula_run_width(draw: ImageDraw.ImageDraw, parts: list[tuple[str, int, int, str]]) -> int:
    width = 0
    for text, size, _offset, font_path in parts:
        width += full.text_box(draw, text, full.font(font_path, size))[0]
    return width


def e11_e33_width(draw: ImageDraw.ImageDraw, size: int) -> int:
    sub = max(18, round(size * 0.62))
    return formula_run_width(
        draw,
        [
            ("E", size, 0, full.TIMES_BOLD),
            ("11", sub, round(size * 0.34), full.TIMES_BOLD),
            (" / E", size, 0, full.TIMES_BOLD),
            ("33", sub, round(size * 0.34), full.TIMES_BOLD),
        ],
    )


def parenthesized_e11_e33_width(draw: ImageDraw.ImageDraw, size: int) -> int:
    f = full.font(full.TIMES_BOLD, size)
    return full.text_box(draw, "(", f)[0] + e11_e33_width(draw, size) + full.text_box(draw, ")", f)[0]


def draw_parenthesized_e11_e33(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    size: int,
    *,
    fill: tuple[int, int, int] = full.INK,
) -> int:
    f = full.font(full.TIMES_BOLD, size)
    x, y = xy
    draw.text((x, y), "(", font=f, fill=fill)
    x += full.text_box(draw, "(", f)[0]
    x = draw_e11_e33(draw, (x, y), size, fill=fill)
    draw.text((x, y), ")", font=f, fill=fill)
    return x + full.text_box(draw, ")", f)[0]


def draw_e_symbol(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    subscript: str,
    size: int,
    *,
    fill: tuple[int, int, int] = full.INK,
) -> int:
    sub = max(18, round(size * 0.62))
    return full.draw_formula_run(
        draw,
        xy,
        [
            ("E", size, 0, full.TIMES_BOLD),
            (subscript, sub, round(size * 0.34), full.TIMES_BOLD),
        ],
        fill=fill,
    )


def e_symbol_width(draw: ImageDraw.ImageDraw, subscript: str, size: int) -> int:
    sub = max(18, round(size * 0.62))
    return formula_run_width(
        draw,
        [
            ("E", size, 0, full.TIMES_BOLD),
            (subscript, sub, round(size * 0.34), full.TIMES_BOLD),
        ],
    )


def wrap_lines(draw: ImageDraw.ImageDraw, text: str, max_width: float, fnt) -> list[str]:
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


def draw_base() -> Image.Image:
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*full.SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*full.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*full.PHOTON, 255))
    title = "Reading photon-like shower shapes"
    draw.text((132, 76), title, font=full.font(full.TIMES_BOLD, 86), fill=full.INK)
    draw.line((132, 232, W - 132, 232), fill=(221, 226, 232, 255), width=3)
    full.add_top_right_sphenix_logo_like_slide2(img)
    return img


def draw_context_inset(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    title = "CEMC η–φ tower image"
    tf = full.font(full.TIMES_BOLD, 40)
    tw, _ = full.text_box(draw, title, tf)
    draw.text((x0 + (x1 - x0 - tw) / 2, y0), title, font=tf, fill=full.BLUE)

    info = "Δη × Δφ = 0.025 × 0.025     24,576 towers     96 η × 256 φ"
    inf = full.font(full.TIMES_BOLD, 29)
    iw, _ = full.text_box(draw, info, inf)
    draw.text((x0 + (x1 - x0 - iw) / 2, y0 + 54), info, font=inf, fill=full.INK)

    gx, gy = x0 + 72, y0 + 110
    cell_w = 44
    cell_h = 31
    rows, cols = 7, 12
    for r in range(rows):
        for c in range(cols):
            fill = (250, 244, 221, 255)
            if 2 <= r <= 4 and 5 <= c <= 7:
                fill = (241, 193, 116, 255)
            draw.rectangle(
                (gx + c * cell_w, gy + r * cell_h, gx + (c + 1) * cell_w, gy + (r + 1) * cell_h),
                fill=fill,
                outline=(204, 214, 224, 255),
                width=2,
            )
    draw.rectangle(
        (gx + 5 * cell_w, gy + 2 * cell_h, gx + 8 * cell_w, gy + 5 * cell_h),
        outline=(*full.SPHENIX_BLUE, 255),
        width=6,
    )
    draw.ellipse(
        (gx + 6.5 * cell_w - 7, gy + 3.5 * cell_h - 7, gx + 6.5 * cell_w + 7, gy + 3.5 * cell_h + 7),
        fill=(0, 0, 0, 255),
    )

    # Axes with enough separation to read as detector segmentation, not data.
    axis = (102, 119, 132)
    draw.line((gx, gy + rows * cell_h + 34, gx + cols * cell_w, gy + rows * cell_h + 34), fill=(*axis, 255), width=3)
    draw.line((gx - 34, gy + rows * cell_h, gx - 34, gy), fill=(*axis, 255), width=3)
    draw_arrow(draw, (gx + cols * cell_w - 22, gy + rows * cell_h + 34), (gx + cols * cell_w + 34, gy + rows * cell_h + 34), fill=axis, width=3, head=12)
    draw_arrow(draw, (gx - 34, gy + 22), (gx - 34, gy - 34), fill=axis, width=3, head=12)
    draw.text((gx + cols * cell_w + 44, gy + rows * cell_h + 17), "φ", font=full.font(full.TIMES_BOLD, 36), fill=axis)
    draw.text((gx - 52, gy - 78), "η", font=full.font(full.TIMES_BOLD, 36), fill=axis)
    label = "local seed-centered window"
    lf = full.font(full.TIMES_ITALIC, 29)
    lw, _ = full.text_box(draw, label, lf)
    draw.text((gx + (cols * cell_w - lw) / 2, gy + rows * cell_h + 58), label, font=lf, fill=full.MUTED)


def draw_rhs_context_block(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    row_x = x0 - 46
    row_w = x1 - row_x
    context_f = full.font(full.TIMES, 43)
    context = "Tower size: Δη × Δφ = 0.025 × 0.025"
    cw, _ = full.text_box(draw, context, context_f)
    draw.text((row_x + (row_w - cw) / 2, y0 + 28), context, font=context_f, fill=full.INK)

    callout = (row_x, y0 + 112, x1, y0 + 404)
    draw.rounded_rectangle(callout, radius=14, fill=(248, 251, 253, 255), outline=(218, 230, 239, 255), width=3)

    # Large conceptual ratio definition:
    # E11/E33 = center tower energy / local 3x3 energy
    eq_size = 48
    frac_f = full.font(full.TIMES_BOLD, 45)
    small_f = frac_f
    left_w = e11_e33_width(draw, eq_size)
    eq_f = full.font(full.TIMES_BOLD, eq_size)
    eq_w = full.text_box(draw, " = ", eq_f)[0]
    num = "center tower energy"
    den_left = "local "
    den_mid = "3×3"
    den_right = " energy"
    num_w = full.text_box(draw, num, frac_f)[0]
    den_w = (
        full.text_box(draw, den_left, frac_f)[0]
        + full.text_box(draw, den_mid, small_f)[0]
        + full.text_box(draw, den_right, frac_f)[0]
    )
    frac_w = max(num_w, den_w) + 28
    total_w = left_w + eq_w + frac_w
    start_x = callout[0] + (callout[2] - callout[0] - total_w) / 2
    y_base = callout[1] + 72
    x = draw_e11_e33(draw, (round(start_x), y_base + 32), eq_size, fill=full.BLUE)
    draw.text((x, y_base + 32), " = ", font=eq_f, fill=full.INK)
    x += eq_w
    frac_x = x + 14
    num_x = frac_x + (frac_w - num_w) / 2
    den_x = frac_x + (frac_w - den_w) / 2
    draw.text((num_x, y_base), num, font=frac_f, fill=full.INK)
    draw.line((frac_x, y_base + 70, frac_x + frac_w, y_base + 70), fill=(*full.INK, 255), width=4)
    dx = den_x
    draw.text((dx, y_base + 88), den_left, font=frac_f, fill=full.INK)
    dx += full.text_box(draw, den_left, frac_f)[0]
    draw.text((dx, y_base + 88), den_mid, font=small_f, fill=full.INK)
    dx += full.text_box(draw, den_mid, small_f)[0]
    draw.text((dx, y_base + 88), den_right, font=frac_f, fill=full.INK)

    def_f = full.font(full.TIMES, 44)
    def_y = callout[3] + 46
    definitions = [
        ("11", "center tower at CoG"),
        ("33", "local 3×3 energy sum"),
    ]
    for subscript, body in definitions:
        bullet_r = 5
        draw.ellipse(
            (row_x + 62, def_y + 22 - bullet_r, row_x + 62 + 2 * bullet_r, def_y + 22 + bullet_r),
            fill=(*full.BLUE, 255),
        )
        x = row_x + 92
        x = draw_e_symbol(draw, (x, def_y - 2), subscript, 43, fill=full.BLUE)
        colon_f = full.font(full.TIMES_BOLD, 43)
        draw.text((x + 4, def_y), ":", font=colon_f, fill=full.INK)
        x += full.text_box(draw, ":", colon_f)[0] + 18
        draw.text((x, def_y - 1), body, font=def_f, fill=full.INK)
        def_y += 66


def draw_local_maps(draw: ImageDraw.ImageDraw, x0: int, y0: int) -> None:
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
    cell = 104
    gap = 98
    label1_prefix = "Local EMCal tower slice:"
    label2 = "Does the center tower dominate the nearby 3×3 neighborhood?"
    lf1 = full.font(full.TIMES_BOLD, 48)
    lf2 = full.font(full.TIMES_ITALIC, 43)
    prefix_w, _ = full.text_box(draw, label1_prefix, lf1)
    formula_w = parenthesized_e11_e33_width(draw, 48)
    ltw2, _ = full.text_box(draw, label2, lf2)
    pair_w = 2 * 5 * cell + gap
    header_w = prefix_w + 22 + formula_w
    hx = x0 + (pair_w - header_w) / 2
    hy = y0 - 148
    draw.text((hx, hy), label1_prefix, font=lf1, fill=full.BLUE)
    draw_parenthesized_e11_e33(draw, (round(hx + prefix_w + 22), hy), 48, fill=full.BLUE)
    draw.text((x0 + (pair_w - ltw2) / 2, y0 - 82), label2, font=lf2, fill=full.MUTED)

    left = full.draw_small_tower_grid(
        draw,
        (x0, y0),
        cell,
        values,
        highlight=(2, 2, 3, 3),
        highlight_color=GEOM_TOWER_GREEN,
        secondary=(1, 1, 4, 4),
        secondary_color=GEOM_CORE_PURPLE,
    )
    right = full.draw_small_tower_grid(
        draw,
        (x0 + 5 * cell + gap, y0),
        cell,
        diffuse,
        highlight=(2, 2, 3, 3),
        highlight_color=GEOM_TOWER_GREEN,
        secondary=(1, 1, 4, 4),
        secondary_color=GEOM_CORE_PURPLE,
    )
    for grid in (left, right):
        gx0, gy0, _, _ = grid
        draw.rectangle(
            (gx0 + cell, gy0 + cell, gx0 + 4 * cell, gy0 + 4 * cell),
            outline=(*GEOM_CORE_PURPLE, 255),
            width=9,
        )
        draw.rectangle(
            (gx0 + 2 * cell, gy0 + 2 * cell, gx0 + 3 * cell, gy0 + 3 * cell),
            outline=(*GEOM_TOWER_GREEN, 255),
            width=9,
        )

        core_f = full.font(full.TIMES_BOLD, 37)
        core_label = "3×3 local core"
        ctw, cth = full.text_box(draw, core_label, core_f)
        draw.text(
            (gx0 + 2.5 * cell - ctw / 2, gy0 + cell - cth - 6),
            core_label,
            font=core_f,
            fill=GEOM_CORE_PURPLE,
        )

        tower_f = full.font(full.TIMES_BOLD, 37)
        tower_lines = ["1×1 center", "tower"]
        line_gap = 1
        line_heights = [full.text_box(draw, line, tower_f)[1] for line in tower_lines]
        total_h = sum(line_heights) + line_gap * (len(tower_lines) - 1)
        ty = gy0 + 3 * cell + 12
        for line, lh in zip(tower_lines, line_heights):
            ttw, _ = full.text_box(draw, line, tower_f)
            draw.text(
                (gx0 + 2.5 * cell - ttw / 2, ty),
                line,
                font=tower_f,
                fill=GEOM_TOWER_GREEN,
            )
            ty += lh + line_gap

        cog_f = full.font(full.TIMES_BOLD, 32)
        cog = "CoG"
        cog_tw, cog_th = full.text_box(draw, cog, cog_f)
        draw.text(
            (gx0 + 2.5 * cell - cog_tw / 2, gy0 + 2.5 * cell - 42),
            cog,
            font=cog_f,
            fill=full.INK,
        )

    for caption, grid, color in (("photon-like", left, PHOTON_LABEL_RED), ("background-like", right, full.SPHENIX_BLUE)):
        gx0, _, gx1, gy1 = grid
        cf = full.font(full.TIMES_BOLD, 42)
        tw, _ = full.text_box(draw, caption, cf)
        draw.text((gx0 + (gx1 - gx0 - tw) / 2, gy1 + 34), caption, font=cf, fill=color)


def build_slide() -> Image.Image:
    img = draw_base()
    draw = ImageDraw.Draw(img, "RGBA")
    card = (132, 258, W - 132, 1282)
    draw.rounded_rectangle(card, radius=12, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((card[0], card[1], card[0] + 12, card[3]), radius=6, fill=(*full.PHOTON_DARK, 225))

    # Keep the sequence cue without adding definitions.
    badge = (card[0] + 60, card[1] + 48, card[0] + 132, card[1] + 120)
    draw.ellipse(badge, fill=(255, 255, 255, 255), outline=(*full.PHOTON_DARK, 235), width=4)
    draw_centered(draw, "1", ((badge[0] + badge[2]) / 2, (badge[1] + badge[3]) / 2), full.font(full.TIMES_BOLD, 40), full.PHOTON_DARK)

    # Main local maps.
    draw_local_maps(draw, card[0] + 196, card[1] + 264)

    # Right-side context and interpretation, text-only by design.
    rhs = (card[0] + 1450, card[1] + 138, card[2] - 88, card[1] + 802)
    draw_rhs_context_block(draw, rhs)

    body_w = card[2] - card[0]
    strip_w = int(body_w * 0.84)
    strip_x0 = card[0] + (body_w - strip_w) // 2
    takeaway = (strip_x0, card[3] - 134, strip_x0 + strip_w, card[3] - 34)
    draw.rounded_rectangle(takeaway, radius=10, fill=(253, 251, 246, 255), outline=(218, 207, 185, 255), width=2)
    prefix = "Photon-like:"
    pf = full.font(full.TIMES_BOLD, 46)
    bf = full.font(full.TIMES, 46)
    body = " energy concentrated in the center of the local 3×3 tower neighborhood."
    total_w = full.text_box(draw, prefix, pf)[0] + full.text_box(draw, body, bf)[0]
    y = takeaway[1] + ((takeaway[3] - takeaway[1]) - full.text_box(draw, prefix, pf)[1]) / 2 - 1
    x = takeaway[0] + ((takeaway[2] - takeaway[0]) - total_w) / 2
    draw.text((x, y), prefix, font=pf, fill=full.INK)
    x += full.text_box(draw, prefix, pf)[0]
    draw.text((x, y), body, font=bf, fill=full.INK)

    full.draw_hp2026_identity_footer(img)
    return img.convert("RGB")


def write_script() -> None:
    SCRIPT_PATH.write_text(
        """# Slide 8 script - compact local core

Here I only want to introduce one representative shower-shape handle. In the local EMCal tower image, the question is whether the center tower dominates the nearby 3×3 neighborhood.

The photon-like example keeps most of the energy concentrated at the center of the local image, while the background-like example spreads more energy across nearby towers. That is what E11 over E33 captures: center-tower dominance inside the local 3×3 region.
""",
        encoding="utf-8",
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    img = build_slide()
    img.save(PNG_PATH, "PNG")
    write_script()
    MANIFEST_PATH.write_text(
        json.dumps(
            {
                "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
                "google_slides_mutation": False,
                "hp2026_main_header": {
                    "deck": "hp2026_main_talk",
                    "title_font_size": 86,
                    "subtitle_font_size": None,
                    "title_xy": [132, 76],
                    "subtitle_xy": None,
                    "divider_y": 232,
                    "exceptions": ["subtitle removed by explicit user request for this local candidate"],
                },
                "source_baseline": "scripts/slides/hp2026/fulltalk/make_hp2026_shower_no_blur_variants.py",
                "output_png": str(PNG_PATH.relative_to(ROOT)),
                "companion_script": str(SCRIPT_PATH.relative_to(ROOT)),
                "slide_intent": "Simplified Slide 8 local PNG candidate using E11/E33 as the only representative compact-core handle.",
                "e11_e33_definition_source": {
                    "code": "src/PhotonClusterBuilder.cc centers a local E77 array on the floored shower-shape CoG, sets e11 = E77[3][3], and sums e33 over di<=1 and dj<=1.",
                    "paper": "PPG12 current draft describes E11/E33 as central tower energy relative to the surrounding 3x3 tower energy.",
                },
                "removed": [
                    "Core fraction = E2x2 / Ecluster",
                    "2x2 core outline and label",
                    "right-side definition table",
                    "prompt-photon-like cluster prose box",
                    "E2x2 mention on the main slide",
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(PNG_PATH)
    print(SCRIPT_PATH)
    print(MANIFEST_PATH)


if __name__ == "__main__":
    main()

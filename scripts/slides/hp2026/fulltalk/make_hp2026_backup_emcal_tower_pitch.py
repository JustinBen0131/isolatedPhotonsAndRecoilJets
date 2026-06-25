#!/usr/bin/env python3
"""Generate a full-slide PNG backup explaining EMCal tower pitch.

This is a local slide candidate only. It does not mutate Google Slides.
"""

from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw, ImageEnhance, ImageFilter, ImageFont

import make_hp2026_fulltalk_candidates as full


ROOT = full.ROOT
OUT_DIR = ROOT / "outputs/manual-20260621-emcal-tower-pitch-backup"
OUT_PNG = OUT_DIR / "backup_emcal_tower_pitch_audience.png"
OUT_SCRIPT = OUT_DIR / "backup_emcal_tower_pitch_audience_script.md"
OUT_META = OUT_DIR / "backup_emcal_tower_pitch_audience_manifest.json"
OUT_HEADER = OUT_DIR / "backup_emcal_tower_pitch_audience.header.json"

W, H = full.W, full.H
MIN_READABLE_FONT_PX = 38
AI_CLEAN_FIGURE = ROOT / "outputs/hp2026_backup_geometry/emcal_barrel_eta_phi_ai_clean_fullbarrel.png"
SOURCE_FIGURE = ROOT / "outputs/hp2026_backup_geometry/emcal_barrel_eta_phi_original_source.png"
FALLBACK_FIGURE = ROOT / "outputs/hp2026_backup_geometry/emcal_barrel_eta_phi_crop_clean.png"
FIGURE = AI_CLEAN_FIGURE if AI_CLEAN_FIGURE.exists() else SOURCE_FIGURE if SOURCE_FIGURE.exists() else FALLBACK_FIGURE


def f(path: Path, size: int) -> ImageFont.FreeTypeFont:
    return full.font(path, max(size, MIN_READABLE_FONT_PX))


def text_size(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.ImageFont) -> tuple[int, int]:
    return full.text_box(draw, text, font)


def center_text(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    text: str,
    font: ImageFont.ImageFont,
    fill: tuple[int, int, int] = full.INK,
    *,
    line_gap: int = 8,
) -> None:
    lines = text.split("\n")
    heights = [text_size(draw, line, font)[1] for line in lines]
    total_h = sum(heights) + line_gap * max(0, len(lines) - 1)
    y = (box[1] + box[3] - total_h) / 2
    for line, height in zip(lines, heights):
        tw, _ = text_size(draw, line, font)
        draw.text(((box[0] + box[2] - tw) / 2, y), line, font=font, fill=fill)
        y += height + line_gap


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    text: str,
    xy: tuple[int, int],
    max_width: int,
    font: ImageFont.ImageFont,
    fill: tuple[int, int, int],
    *,
    line_gap: int = 8,
) -> int:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        trial = word if not current else f"{current} {word}"
        if text_size(draw, trial, font)[0] <= max_width:
            current = trial
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)

    x, y = xy
    for line in lines:
        draw.text((x, y), line, font=font, fill=fill)
        y += text_size(draw, line, font)[1] + line_gap
    return y


def erase_red_annotations(fig: Image.Image) -> Image.Image:
    """Remove red annotation arrows/text baked into the legacy source crop."""
    px = []
    for r, g, b, a in fig.getdata():
        is_red_annotation = r > 125 and r > 1.35 * max(g, 1) and r > 1.35 * max(b, 1)
        if is_red_annotation:
            px.append((255, 255, 255, 0))
        else:
            px.append((r, g, b, a))
    cleaned = Image.new("RGBA", fig.size)
    cleaned.putdata(px)
    return cleaned


def draw_shadowed_card(
    base: Image.Image,
    box: tuple[int, int, int, int],
    *,
    radius: int = 14,
    fill: tuple[int, int, int, int] = (255, 255, 255, 255),
    outline: tuple[int, int, int, int] = (*full.PANEL_EDGE, 255),
    shadow_alpha: int = 18,
) -> None:
    shadow = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow, "RGBA")
    x0, y0, x1, y1 = box
    sd.rounded_rectangle((x0 + 8, y0 + 10, x1 + 8, y1 + 10), radius=radius, fill=(20, 40, 60, shadow_alpha))
    shadow = shadow.filter(ImageFilter.GaussianBlur(9))
    base.alpha_composite(shadow)
    draw = ImageDraw.Draw(base, "RGBA")
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=3)


def draw_arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], color: tuple[int, int, int], width: int = 5) -> None:
    draw.line((*start, *end), fill=(*color, 235), width=width)
    angle = math.atan2(end[1] - start[1], end[0] - start[0])
    size = 18
    pts = [
        end,
        (int(end[0] - size * math.cos(angle - 0.48)), int(end[1] - size * math.sin(angle - 0.48))),
        (int(end[0] - size * math.cos(angle + 0.48)), int(end[1] - size * math.sin(angle + 0.48))),
    ]
    draw.polygon(pts, fill=(*color, 245))


def draw_curve_arrow(
    draw: ImageDraw.ImageDraw,
    bbox: tuple[int, int, int, int],
    start: int,
    end: int,
    color: tuple[int, int, int],
    *,
    width: int = 7,
) -> None:
    draw.arc(bbox, start=start, end=end, fill=(*color, 235), width=width)
    angle = math.radians(end)
    cx = (bbox[0] + bbox[2]) / 2
    cy = (bbox[1] + bbox[3]) / 2
    rx = (bbox[2] - bbox[0]) / 2
    ry = (bbox[3] - bbox[1]) / 2
    end_pt = (int(cx + rx * math.cos(angle)), int(cy + ry * math.sin(angle)))
    tangent = angle + math.pi / 2
    size = 16
    pts = [
        end_pt,
        (int(end_pt[0] - size * math.cos(tangent - 0.45)), int(end_pt[1] - size * math.sin(tangent - 0.45))),
        (int(end_pt[0] - size * math.cos(tangent + 0.45)), int(end_pt[1] - size * math.sin(tangent + 0.45))),
    ]
    draw.polygon(pts, fill=(*color, 245))


def draw_formula_callout(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    title: str,
    formula: str,
    note: str,
    color: tuple[int, int, int],
) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 250), outline=(213, 224, 234, 255), width=2)
    draw.rounded_rectangle((x0, y0, x0 + 12, y1), radius=7, fill=(*color, 255))
    draw.text((x0 + 30, y0 + 18), title, font=f(full.TIMES_BOLD, 25), fill=full.INK)
    draw.text((x0 + 30, y0 + 56), formula, font=f(full.TIMES_BOLD, 34), fill=color)
    draw.text((x0 + 30, y0 + 104), note, font=f(full.TIMES, 22), fill=full.MUTED)


def draw_coordinate_schematic(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    """Draw a crisp schematic of a cylindrical EMCal tower map."""
    x0, y0, x1, y1 = box

    sketch = (x0 + 36, y0 + 34, x0 + 720, y1 - 28)
    sx0, sy0, sx1, sy1 = sketch
    cy = (sy0 + sy1) // 2 + 4
    front_cx = sx0 + 162
    back_cx = sx0 + 562
    outer_rx, outer_ry = 116, 154
    inner_rx, inner_ry = 72, 96

    # Back and front barrel outlines.
    outline = (53, 62, 72)
    light = (178, 190, 202)
    draw.ellipse((back_cx - outer_rx, cy - outer_ry, back_cx + outer_rx, cy + outer_ry), outline=(*light, 190), width=4)
    draw.ellipse((back_cx - inner_rx, cy - inner_ry, back_cx + inner_rx, cy + inner_ry), outline=(*light, 130), width=3)
    top_front = (front_cx, cy - outer_ry)
    bot_front = (front_cx, cy + outer_ry)
    top_back = (back_cx, cy - outer_ry)
    bot_back = (back_cx, cy + outer_ry)
    draw.line((top_front, top_back), fill=(*outline, 210), width=5)
    draw.line((bot_front, bot_back), fill=(*outline, 210), width=5)
    for offset in [-72, -36, 0, 36, 72]:
        y_front = cy + offset
        y_back = cy + int(offset * 0.82)
        draw.line((front_cx - 2, y_front, back_cx, y_back), fill=(152, 165, 176, 170), width=2)

    draw.ellipse((front_cx - outer_rx, cy - outer_ry, front_cx + outer_rx, cy + outer_ry), outline=(*outline, 235), width=5)
    draw.ellipse((front_cx - inner_rx, cy - inner_ry, front_cx + inner_rx, cy + inner_ry), outline=(144, 157, 168, 220), width=4)

    # Highlight one eta ring and one phi wedge.
    ring_x0 = back_cx - 78
    ring_x1 = back_cx + 4
    ring_poly = [
        (ring_x0, cy - outer_ry + 12),
        (ring_x1, cy - outer_ry + 22),
        (ring_x1 + 16, cy + outer_ry - 22),
        (ring_x0 + 10, cy + outer_ry - 12),
    ]
    draw.polygon(ring_poly, fill=(*full.PHOTON, 82))
    draw.line(ring_poly + [ring_poly[0]], fill=(*full.PHOTON_DARK, 210), width=4)
    draw_curve_arrow(
        draw,
        (front_cx - outer_rx + 14, cy - outer_ry + 18, front_cx + outer_rx - 14, cy + outer_ry - 18),
        start=-35,
        end=45,
        color=full.TEAL,
        width=8,
    )

    center = (front_cx, cy)
    tower = (ring_x0 + 58, cy - 30)
    draw.ellipse((tower[0] - 9, tower[1] - 9, tower[0] + 9, tower[1] + 9), fill=(0, 0, 0, 255))
    draw.line((center, tower), fill=(48, 58, 69, 220), width=4)

    # Axes and geometric labels.
    draw_arrow(draw, (front_cx - 32, cy + outer_ry + 34), (back_cx + 120, cy + outer_ry + 34), full.SPHENIX_BLUE, width=5)
    draw.text((back_cx + 130, cy + outer_ry + 13), "z", font=f(full.TIMES_ITALIC, 31), fill=full.SPHENIX_BLUE)
    draw.text((front_cx + 16, cy + outer_ry + 47), "barrel axis", font=f(full.TIMES, 20), fill=full.MUTED)
    draw.text((center[0] - 24, center[1] - 30), "R", font=f(full.TIMES_ITALIC, 31), fill=full.MUTED)
    draw.text((ring_x0 - 2, cy - outer_ry - 36), "η ring", font=f(full.TIMES_BOLD, 25), fill=full.PHOTON_DARK)
    draw.text((front_cx - 126, cy - outer_ry - 2), "φ wraps around", font=f(full.TIMES_BOLD, 24), fill=full.TEAL)
    draw.text((tower[0] + 16, tower[1] - 4), "(x, y, z)", font=f(full.TIMES_ITALIC, 25), fill=full.INK)

    # The formulas are text, not baked into a low-resolution screenshot.
    callout_x = x0 + 760
    draw_formula_callout(
        draw,
        (callout_x, y0 + 44, x1 - 34, y0 + 180),
        "eta coordinate",
        "η = asinh(z/R)",
        "R = √(x²+y²); η is dimensionless",
        full.SPHENIX_BLUE,
    )
    draw_formula_callout(
        draw,
        (callout_x, y0 + 214, x1 - 34, y0 + 350),
        "azimuthal coordinate",
        "φ = atan2(y,x)",
        "angle around the detector, in radians",
        full.TEAL,
    )

    draw_arrow(draw, (tower[0] + 24, tower[1] - 24), (callout_x - 18, y0 + 112), full.SPHENIX_BLUE, width=4)
    draw_arrow(draw, (front_cx - 56, cy - 122), (callout_x - 18, y0 + 282), full.TEAL, width=4)


def draw_local_grid(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(250, 253, 255, 255), outline=(207, 223, 232, 255), width=2)

    grid_x0 = x0 + 92
    grid_y0 = y0 + 54
    cell = 28
    rows = 7
    cols = 11
    grid_w = cols * cell
    grid_h = rows * cell
    for i in range(cols + 1):
        x = grid_x0 + i * cell
        draw.line((x, grid_y0, x, grid_y0 + grid_h), fill=(199, 219, 231, 255), width=3 if i in (4, 7) else 2)
    for j in range(rows + 1):
        y = grid_y0 + j * cell
        draw.line((grid_x0, y, grid_x0 + grid_w, y), fill=(199, 219, 231, 255), width=3 if j in (2, 5) else 2)

    # Local 3x3 patch around the seed.
    patch = (grid_x0 + 4 * cell, grid_y0 + 2 * cell, grid_x0 + 7 * cell, grid_y0 + 5 * cell)
    draw.rectangle(patch, fill=(247, 180, 35, 70), outline=(*full.SPHENIX_BLUE, 255), width=5)
    seed = (grid_x0 + 5.5 * cell, grid_y0 + 3.5 * cell)
    draw.ellipse((seed[0] - 10, seed[1] - 10, seed[0] + 10, seed[1] + 10), fill=(0, 0, 0, 255))
    draw.text((patch[2] + 22, patch[1] + 42), "seed", font=f(full.TIMES_ITALIC, 27), fill=full.MUTED)

    axis_color = full.TEAL
    axis_y = grid_y0 + grid_h + 22
    draw_arrow(draw, (grid_x0 - 38, axis_y), (grid_x0 + grid_w + 18, axis_y), axis_color, width=5)
    draw.text((grid_x0 + grid_w + 28, axis_y - 20), "η", font=f(full.TIMES_ITALIC, 36), fill=axis_color)
    draw_arrow(draw, (grid_x0 - 38, axis_y - 2), (grid_x0 - 38, grid_y0 - 18), axis_color, width=5)
    draw.text((grid_x0 - 72, grid_y0 - 40), "φ", font=f(full.TIMES_ITALIC, 36), fill=axis_color)

    caption = "local 3×3 tower patch"
    center_text(draw, (x0 + 28, y1 - 42, x1 - 28, y1 - 8), caption, f(full.TIMES_BOLD, 25), full.BLUE, line_gap=4)


def draw_pitch_row(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    label: str,
    expression: str,
    result: str,
    color: tuple[int, int, int],
) -> None:
    draw.rounded_rectangle((x, y, x + 864, y + 96), radius=12, fill=(255, 255, 255, 255), outline=(216, 226, 235, 255), width=2)
    draw.rounded_rectangle((x, y, x + 14, y + 96), radius=8, fill=(*color, 255))
    draw.text((x + 36, y + 20), label, font=f(full.TIMES_BOLD, 29), fill=full.INK)
    draw.text((x + 228, y + 20), expression, font=f(full.TIMES, 30), fill=full.MUTED)
    draw.text((x + 604, y + 20), result, font=f(full.TIMES_BOLD, 31), fill=color)


def load_enhanced_geometry_figure() -> Image.Image:
    fig = Image.open(FIGURE).convert("RGBA")
    if FIGURE != AI_CLEAN_FIGURE:
        # Legacy fallback: use the original barrel drawing as the detector
        # schematic, but drop the old low-resolution formula text. The formulas
        # are redrawn as sharp slide text in the generator.
        if fig.width >= 900 and fig.height >= 500:
            fig = fig.crop((0, 0, 560, fig.height))
        else:
            fig = fig.crop((0, 0, int(fig.width * 0.70), fig.height))
    fig = full.crop_visible(fig, white_threshold=253, pad=8)
    # The source from the older deck is the right image; upscale and sharpen it
    # before placing it as the visual anchor of the generated slide.
    fig = fig.resize((fig.width * 4, fig.height * 4), Image.Resampling.LANCZOS)
    fig = ImageEnhance.Sharpness(fig).enhance(1.35)
    fig = ImageEnhance.Contrast(fig).enhance(1.04)
    return fig


def draw_shower_shape_inset(draw: ImageDraw.ImageDraw, image_box: tuple[int, int, int, int]) -> None:
    fx0, fy0, fx1, fy1 = image_box
    fw = fx1 - fx0
    fh = fy1 - fy0
    target = (int(fx0 + 0.60 * fw), int(fy0 + 0.44 * fh))
    inset = (fx0 + 560, fy1 - 214, fx0 + 982, fy1 - 40)

    draw.rounded_rectangle(inset, radius=12, fill=(255, 255, 255, 236), outline=(202, 219, 231, 255), width=2)
    draw.text((inset[0] + 18, inset[1] + 14), "local shower-shape window", font=f(full.TIMES_BOLD, 24), fill=full.BLUE)
    draw.text((inset[0] + 18, inset[1] + 44), "tower energies around the seed", font=f(full.TIMES_ITALIC, 20), fill=full.MUTED)

    cell = 21
    gx = inset[0] + 26
    gy = inset[1] + 78
    rows = cols = 5
    for i in range(cols + 1):
        x = gx + i * cell
        draw.line((x, gy, x, gy + rows * cell), fill=(196, 216, 228, 255), width=2)
    for j in range(rows + 1):
        y = gy + j * cell
        draw.line((gx, y, gx + cols * cell, y), fill=(196, 216, 228, 255), width=2)
    patch = (gx + cell, gy + cell, gx + 4 * cell, gy + 4 * cell)
    draw.rectangle(patch, fill=(*full.PHOTON, 72), outline=(*full.SPHENIX_BLUE, 255), width=4)
    seed = (gx + int(2.5 * cell), gy + int(2.5 * cell))
    draw.ellipse((seed[0] - 7, seed[1] - 7, seed[0] + 7, seed[1] + 7), fill=(0, 0, 0, 255))
    draw.text((gx + cols * cell + 22, gy + 30), "seed", font=f(full.TIMES_ITALIC, 24), fill=full.MUTED)
    draw.text((gx + cols * cell + 22, gy + 72), "3×3", font=f(full.TIMES_BOLD, 26), fill=full.SPHENIX_BLUE)

    draw_arrow(draw, (inset[2] - 20, inset[1] + 84), target, full.PHOTON_DARK, width=4)


def draw_bottom_chip(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    eyebrow: str,
    main: str,
    sub: str,
    color: tuple[int, int, int],
) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=14, fill=(255, 255, 255, 255), outline=(215, 226, 235, 255), width=2)
    draw.rounded_rectangle((x0, y0, x0 + 14, y1), radius=8, fill=(*color, 255))
    draw.text((x0 + 34, y0 + 18), eyebrow, font=f(full.TIMES_BOLD, 25), fill=color)
    draw_wrapped(draw, main, (x0 + 34, y0 + 56), x1 - x0 - 68, f(full.TIMES_BOLD, 30), full.INK, line_gap=4)
    draw_wrapped(draw, sub, (x0 + 34, y0 + 104), x1 - x0 - 68, f(full.TIMES, 22), full.MUTED, line_gap=3)


def draw_formula_stack(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    draw.text((x0, y0), "Coordinate definitions", font=f(full.TIMES_BOLD, 38), fill=full.INK)
    draw.text((x0, y0 + 50), "tower center at (x, y, z)", font=f(full.TIMES_ITALIC, 27), fill=full.MUTED)

    def row(y: int, color: tuple[int, int, int], lhs: str, rhs: str, note: str) -> None:
        row_h = 112
        draw.rounded_rectangle((x0, y, x1, y + row_h), radius=10, fill=(255, 255, 255, 255), outline=(218, 229, 238, 255), width=2)
        draw.rectangle((x0, y, x0 + 12, y + row_h), fill=(*color, 255))
        draw.text((x0 + 34, y + 13), lhs, font=f(full.TIMES_BOLD, 27), fill=color)
        draw.text((x0 + 34, y + 58), rhs, font=f(full.TIMES_BOLD, 35), fill=full.INK)
        nw = text_size(draw, note, f(full.TIMES, 20))[0]
        draw.text((x1 - nw - 28, y + 62), note, font=f(full.TIMES, 20), fill=full.MUTED)

    row(y0 + 100, full.SPHENIX_BLUE, "eta coordinate", "η = asinh(z/R)", "dimensionless")
    row(y0 + 232, full.TEAL, "azimuthal coordinate", "φ = atan2(y,x)", "radians")

    grid_y = y0 + 390
    draw.text((x0, grid_y), "Local shower-shape window", font=f(full.TIMES_BOLD, 33), fill=full.BLUE)
    draw.text((x0 + 26, grid_y + 42), "tower-energy pattern around the seed", font=f(full.TIMES_ITALIC, 24), fill=full.MUTED)
    gx = x0 + 32
    gy = grid_y + 92
    cell = 32
    rows = cols = 7
    for i in range(cols + 1):
        x = gx + i * cell
        draw.line((x, gy, x, gy + rows * cell), fill=(196, 216, 229, 255), width=2)
    for j in range(rows + 1):
        y = gy + j * cell
        draw.line((gx, y, gx + cols * cell, y), fill=(196, 216, 229, 255), width=2)
    patch = (gx + 2 * cell, gy + 2 * cell, gx + 5 * cell, gy + 5 * cell)
    draw.rectangle(patch, fill=(*full.PHOTON, 70), outline=(*full.SPHENIX_BLUE, 255), width=5)
    seed = (gx + int(3.5 * cell), gy + int(3.5 * cell))
    draw.ellipse((seed[0] - 8, seed[1] - 8, seed[0] + 8, seed[1] + 8), fill=(0, 0, 0, 255))
    draw.text((gx + cols * cell + 28, gy + 92), "seed", font=f(full.TIMES_ITALIC, 27), fill=full.MUTED)
    draw.text((gx + cols * cell + 28, gy + 134), "3×3", font=f(full.TIMES_BOLD, 31), fill=full.SPHENIX_BLUE)
    axis_color = full.TEAL
    draw_arrow(draw, (gx, gy + rows * cell + 34), (gx + cols * cell + 52, gy + rows * cell + 34), axis_color, width=4)
    draw.text((gx + cols * cell + 62, gy + rows * cell + 18), "η", font=f(full.TIMES_ITALIC, 34), fill=axis_color)
    draw_arrow(draw, (gx - 28, gy + rows * cell), (gx - 28, gy + 4), axis_color, width=4)
    draw.text((gx - 66, gy + 4), "φ", font=f(full.TIMES_ITALIC, 34), fill=axis_color)


def draw_clean_footer(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    illinois_path = full.TITLE_ASSET_DIR / "illinois_logo_fullcolor_rgb.png"
    hp_path = full.TITLE_ASSET_DIR / "hp2026_indico_logo.png"
    cy = 1384
    if illinois_path.exists():
        illinois = full.crop_visible(Image.open(illinois_path).convert("RGBA"), white_threshold=252)
        illinois = full.fit(illinois, 54, 62)
        base.alpha_composite(illinois, (30, cy - illinois.height // 2))
    draw.text((104, cy - 17), "Justin Bennett", font=f(full.TIMES, 31), fill=(43, 49, 57))

    center = "Hard Probes 2026 / June 24, 2026"
    center_font = f(full.TIMES_BOLD, 31)
    cw, ch = full.text_box(draw, center, center_font)
    hp = None
    if hp_path.exists():
        hp = full.white_to_transparent(Image.open(hp_path).convert("RGBA"), threshold=246)
        hp = full.fit(hp, 116, 60)
    group_w = (hp.width + 20 if hp is not None else 0) + cw
    gx = (W - group_w) // 2
    if hp is not None:
        base.alpha_composite(hp, (gx, cy - hp.height // 2))
        gx += hp.width + 20
    draw.text((gx, cy - ch // 2 - 1), center, font=center_font, fill=(43, 49, 57))


def draw_bottom_fact(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    w: int,
    title: str,
    main: str,
    detail: str,
    color: tuple[int, int, int],
) -> None:
    draw.rectangle((x, y, x + 10, y + 124), fill=(*color, 255))
    draw.text((x + 34, y + 0), title, font=f(full.TIMES_BOLD, 25), fill=color)
    draw.text((x + 34, y + 43), main, font=f(full.TIMES_BOLD, 30), fill=full.INK)
    draw.text((x + 34, y + 86), detail, font=f(full.TIMES, 22), fill=full.MUTED)


def draw_slide() -> Image.Image:
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")

    # HPslides_v1-like white canvas, large serif title, logos/footer.
    title = "Backup: EMCal geometry and tower segmentation"
    draw.text((132, 72), title, font=f(full.TIMES_BOLD, 78), fill=full.INK)
    full.add_top_right_sphenix_logo_like_slide2(img)

    # Main coordinate figure.
    if FIGURE.exists():
        fig = load_enhanced_geometry_figure()
        image_box = full.paste_fit(img, fig, (110, 180, 1490, 1006))
        # A clean connector from the tower center in the schematic to the local
        # shower-shape map in the right column.
        if FIGURE == AI_CLEAN_FIGURE:
            tower_dot = (
                image_box[0] + int(0.545 * (image_box[2] - image_box[0])),
                image_box[1] + int(0.397 * (image_box[3] - image_box[1])),
            )
        else:
            tower_dot = (
                image_box[0] + int(0.72 * (image_box[2] - image_box[0])),
                image_box[1] + int(0.42 * (image_box[3] - image_box[1])),
            )
        draw_arrow(draw, tower_dot, (1512, 770), full.PHOTON_DARK, width=4)
    draw_formula_stack(draw, (1510, 205, W - 132, 1000))

    # Bottom readout row.
    y0 = 1084
    x = 132
    fact_w = 690
    draw_bottom_fact(
        draw,
        x,
        y0,
        fact_w,
        "Angular tower pitch",
        "Δη × Δφ, not cm × cm",
        "η dimensionless; φ in radians",
        full.PHOTON_DARK,
    )
    x += 760
    draw_bottom_fact(
        draw,
        x,
        y0,
        fact_w,
        "η segmentation",
        "Δη ≈ 2.2 / 96 = 0.0229",
        "quoted as ~0.025 granularity",
        full.SPHENIX_BLUE,
    )
    x += 760
    draw_bottom_fact(
        draw,
        x,
        y0,
        fact_w,
        "φ segmentation",
        "Δφ ≈ 2π / 256 = 0.0245 rad",
        "256 bins around full azimuth",
        full.TEAL,
    )

    draw_clean_footer(img)
    return img.convert("RGB")


def write_script() -> None:
    OUT_SCRIPT.write_text(
        """# Backup - EMCal geometry and tower segmentation

This backup slide gives the geometry behind the EMCal tower notation and the local shower-shape windows.

The EMCal is a cylindrical detector layer, so the readout is naturally described in eta and phi. Eta is dimensionless, and phi is an angle measured in radians. The shower-shape variables use small local tower windows around the cluster seed in this eta-phi map.

The CEMC is segmented into 96 eta bins across about |eta| less than 1.1, so Delta eta is about 2.2 divided by 96, or 0.0229. Around the full azimuth there are 256 phi bins, so Delta phi is 2 pi divided by 256, or 0.0245 radians. This is why the paper summarizes the EMCal granularity as approximately 0.025 by 0.025 in eta-phi.
""",
        encoding="utf-8",
    )


def write_metadata() -> None:
    OUT_HEADER.write_text(
        json.dumps(
            {
                "hp2026_png_candidate": {
                    "size_px": [W, H],
                    "style": "HPslides_v1 white canvas with footer logos",
                    "slide_number_baked_in": False,
                }
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    OUT_META.write_text(
        json.dumps(
            {
                "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
                "google_slides_mutation": False,
                "output_png": str(OUT_PNG.relative_to(ROOT)),
                "speaker_script": str(OUT_SCRIPT.relative_to(ROOT)),
                "source_script": str(Path(__file__).relative_to(ROOT)),
                "source_figure": str(FIGURE.relative_to(ROOT)),
                "physics_content": {
                    "eta_bins": 96,
                    "phi_bins": 256,
                    "eta_range_used": "approximately |eta| <= 1.1, total span about 2.2",
                    "delta_eta": "2.2 / 96 = 0.0229",
                    "delta_phi": "2*pi / 256 = 0.0245 rad",
                    "paper_rounding": "approximately 0.025 x 0.025 eta-phi segmentation",
                    "units_answer": "eta is dimensionless; phi is in radians; this is not cm x cm",
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    img = draw_slide()
    img.save(OUT_PNG, quality=95)
    write_script()
    write_metadata()
    print(OUT_PNG)
    print(OUT_SCRIPT)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Render the HP2026 Slide 2 opening/motivation PNG candidate.

The slide is deterministic: text and schematic geometry are drawn with Pillow,
with no AI-generated logos or copied physics figures. It is intended as a local
candidate only; Google Slides insertion is a separate approval step.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont


W, H = 2560, 1440

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "AGENTS.md").exists())
DEFAULT_WORKSPACE = ROOT / "outputs/manual-20260601-hp2026-opening-slide/presentations/hp2026-opening-slide"
DEFAULT_OUTPUT = DEFAULT_WORKSPACE / "output"
TITLE_ASSET_DIR = ROOT / "outputs/manual-20260601-hp2026-title/presentations/hp2026-title-slide/assets"

FONT_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES = FONT_DIR / "Times New Roman.ttf"
TIMES_BOLD = FONT_DIR / "Times New Roman Bold.ttf"
TIMES_ITALIC = FONT_DIR / "Times New Roman Italic.ttf"
ARIAL = FONT_DIR / "Arial.ttf"
ARIAL_BOLD = FONT_DIR / "Arial Bold.ttf"

INK = (18, 22, 28)
MUTED = (71, 79, 91)
LIGHT_MUTED = (117, 126, 138)
BLUE = (19, 41, 75)
SPHENIX_BLUE = (30, 143, 214)
TEAL = (18, 112, 142)
TEAL_SOFT = (102, 173, 190)
PHOTON = (245, 181, 34)
PHOTON_DARK = (196, 124, 10)
SOFT_BG = (252, 253, 254)
PANEL = (246, 249, 252)
PANEL_EDGE = (218, 226, 235)
CARD = (255, 255, 255)
CARD_EDGE = (223, 229, 236)

TITLE = "Why isolated prompt photons?"
SUBTITLE = "A color-neutral hard probe for the p+p baseline at RHIC"
BRIDGE = "Next: how sPHENIX turns EMCal clusters into isolated prompt-photon candidates."

REASONS = [
    (
        "Calibrated hard scale",
        "The photon tags the hard scattering with minimal final-state interaction.",
    ),
    (
        "Isolation sharpens the sample",
        "A quiet cone suppresses decay and fragmentation photons before purity correction.",
    ),
    (
        "p+p baseline for future hard probes",
        "The cross section anchors RHIC comparisons and later gamma-jet measurements.",
    ),
]

CONTEXT_TITLE = "A p+p photon baseline for the sPHENIX hard-probes program"
CONTEXT_SUBTITLE = "RHIC p+p collisions at 200 GeV are the reference system for calibrated hard scattering."
CONTEXT_BRIDGE = "Next: define the prompt-photon object and the isolation requirement."

CONTEXT_REASONS = [
    (
        "Hard probes need a reference",
        "p+p establishes the baseline before interpreting nuclear modifications in heavy-ion data.",
    ),
    (
        "Photons preserve the hard scale",
        "Prompt photons escape color final-state interactions; isolation suppresses decay and fragmentation backgrounds.",
    ),
    (
        "sPHENIX makes it measurable",
        "EMCal shower shape, isolation, and purity plus correction steps turn candidates into a cross section.",
    ),
]


def font(path: Path, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(path), size)


def text_box(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.ImageFont) -> tuple[int, int]:
    box = draw.textbbox((0, 0), text, font=fnt)
    return box[2] - box[0], box[3] - box[1]


def open_rgba(path: Path) -> Image.Image:
    return Image.open(path).convert("RGBA")


def crop_visible(img: Image.Image, white_threshold: int = 248) -> Image.Image:
    rgba = img.convert("RGBA")
    px = rgba.load()
    min_x, min_y, max_x, max_y = rgba.width, rgba.height, -1, -1
    for y in range(rgba.height):
        for x in range(rgba.width):
            r, g, b, a = px[x, y]
            visible = a > 8 and not (r >= white_threshold and g >= white_threshold and b >= white_threshold)
            if visible:
                min_x = min(min_x, x)
                min_y = min(min_y, y)
                max_x = max(max_x, x)
                max_y = max(max_y, y)
    if max_x < min_x:
        return rgba
    pad = 8
    return rgba.crop(
        (
            max(0, min_x - pad),
            max(0, min_y - pad),
            min(rgba.width, max_x + pad + 1),
            min(rgba.height, max_y + pad + 1),
        )
    )


def fit(img: Image.Image, max_w: int, max_h: int) -> Image.Image:
    scale = min(max_w / img.width, max_h / img.height)
    size = (max(1, int(img.width * scale)), max(1, int(img.height * scale)))
    return img.resize(size, Image.Resampling.LANCZOS)


def paste_fit(base: Image.Image, img: Image.Image, box: tuple[int, int, int, int], anchor: str = "center") -> None:
    x0, y0, x1, y1 = box
    fitted = fit(img, x1 - x0, y1 - y0)
    if anchor == "left":
        x = x0
    elif anchor == "right":
        x = x1 - fitted.width
    else:
        x = x0 + ((x1 - x0) - fitted.width) // 2
    y = y0 + ((y1 - y0) - fitted.height) // 2
    base.alpha_composite(fitted, (x, y))


def load_sphenix_logo() -> Image.Image | None:
    logo = TITLE_ASSET_DIR / "sphenix-logo-white-bg_0.png"
    if not logo.exists():
        return None
    return crop_visible(open_rgba(logo), white_threshold=252)


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    text: str,
    xy: tuple[int, int],
    max_width: int,
    fnt: ImageFont.ImageFont,
    fill: tuple[int, int, int],
    line_gap: int,
) -> int:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        trial = word if not current else f"{current} {word}"
        if text_box(draw, trial, fnt)[0] <= max_width:
            current = trial
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)

    x, y = xy
    for line in lines:
        draw.text((x, y), line, font=fnt, fill=fill)
        y += text_box(draw, line, fnt)[1] + line_gap
    return y


def alpha_layer(base: Image.Image, callback) -> None:
    layer = Image.new("RGBA", base.size, (0, 0, 0, 0))
    callback(ImageDraw.Draw(layer, "RGBA"))
    base.alpha_composite(layer)


def feynman_points(
    start: tuple[float, float],
    end: tuple[float, float],
    amplitude: float,
    cycles: float,
    steps: int,
) -> list[tuple[float, float]]:
    sx, sy = start
    ex, ey = end
    dx, dy = ex - sx, ey - sy
    length = math.hypot(dx, dy)
    if length == 0:
        return [start]
    nx, ny = -dy / length, dx / length
    pts = []
    for i in range(steps + 1):
        t = i / steps
        wave = math.sin(t * cycles * math.tau) * amplitude
        pts.append((sx + dx * t + nx * wave, sy + dy * t + ny * wave))
    return pts


def draw_polyline(
    draw: ImageDraw.ImageDraw,
    pts: list[tuple[float, float]],
    fill: tuple[int, int, int, int] | tuple[int, int, int],
    width: int,
) -> None:
    draw.line([(round(x), round(y)) for x, y in pts], fill=fill, width=width, joint="curve")


def interp(a: tuple[float, float], b: tuple[float, float], t: float) -> tuple[float, float]:
    return (a[0] * (1 - t) + b[0] * t, a[1] * (1 - t) + b[1] * t)


def draw_detector_face(base: Image.Image, impact: tuple[int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    tl = (1148, 392)
    tr = (1358, 454)
    br = (1304, 664)
    bl = (1094, 604)
    draw.polygon([tl, tr, br, bl], fill=(237, 244, 249, 255), outline=(172, 194, 212, 255))

    for i in range(1, 5):
        t = i / 5
        p1 = interp(tl, bl, t)
        p2 = interp(tr, br, t)
        draw.line((p1, p2), fill=(197, 211, 223, 255), width=3)
    for i in range(1, 6):
        t = i / 6
        p1 = interp(tl, tr, t)
        p2 = interp(bl, br, t)
        draw.line((p1, p2), fill=(197, 211, 223, 255), width=3)

    glow = Image.new("RGBA", base.size, (0, 0, 0, 0))
    gdraw = ImageDraw.Draw(glow, "RGBA")
    x, y = impact
    for r, alpha in ((74, 42), (50, 72), (28, 110)):
        gdraw.ellipse((x - r, y - r, x + r, y + r), fill=(*PHOTON, alpha))
    glow = glow.filter(ImageFilter.GaussianBlur(8))
    base.alpha_composite(glow)
    draw = ImageDraw.Draw(base, "RGBA")
    draw.ellipse((x - 22, y - 22, x + 22, y + 22), fill=(*PHOTON, 230), outline=(*PHOTON_DARK, 230), width=3)
    draw.text((1138, 686), "EMCal", font=font(TIMES_ITALIC, 31), fill=LIGHT_MUTED)


def draw_icon(draw: ImageDraw.ImageDraw, kind: int, center: tuple[int, int]) -> None:
    cx, cy = center
    if kind == 0:
        pts = feynman_points((cx - 34, cy + 12), (cx + 34, cy - 12), 7, 3.6, 80)
        draw_polyline(draw, pts, (*PHOTON_DARK, 255), 4)
        draw.ellipse((cx - 44, cy - 44, cx + 44, cy + 44), outline=(195, 207, 220, 255), width=3)
        draw.line((cx - 52, cy, cx - 34, cy), fill=(195, 207, 220, 255), width=3)
        draw.line((cx + 34, cy, cx + 52, cy), fill=(195, 207, 220, 255), width=3)
    elif kind == 1:
        draw.pieslice((cx - 47, cy - 47, cx + 47, cy + 47), start=320, end=40, fill=(250, 223, 154, 120), outline=(*PHOTON_DARK, 200), width=3)
        draw.line((cx, cy, cx + 46, cy - 17), fill=(*PHOTON_DARK, 220), width=3)
        draw.line((cx, cy, cx + 46, cy + 17), fill=(*PHOTON_DARK, 220), width=3)
        draw.ellipse((cx - 10, cy - 10, cx + 10, cy + 10), fill=(*PHOTON, 230))
    else:
        draw.line((cx - 42, cy + 34, cx + 44, cy + 34), fill=(73, 104, 128, 255), width=4)
        draw.line((cx - 42, cy + 34, cx - 42, cy - 34), fill=(73, 104, 128, 255), width=4)
        draw.line((cx - 30, cy + 22, cx + 32, cy - 18), fill=(*TEAL, 230), width=5)
        draw.ellipse((cx + 22, cy - 26, cx + 42, cy - 6), fill=(*SPHENIX_BLUE, 220))


def draw_context_icon(draw: ImageDraw.ImageDraw, kind: int, center: tuple[int, int]) -> None:
    cx, cy = center
    if kind == 0:
        draw.line((cx - 46, cy + 36, cx + 48, cy + 36), fill=(72, 101, 122, 255), width=4)
        draw.line((cx - 46, cy + 36, cx - 46, cy - 36), fill=(72, 101, 122, 255), width=4)
        draw.line((cx - 35, cy + 18, cx - 5, cy - 8, cx + 38, cy - 24), fill=(*SPHENIX_BLUE, 230), width=5)
        draw.ellipse((cx + 28, cy - 34, cx + 50, cy - 12), fill=(*SPHENIX_BLUE, 230))
    elif kind == 1:
        pts = feynman_points((cx - 50, cy + 4), (cx + 50, cy - 4), 8, 5.0, 100)
        draw_polyline(draw, pts, (*PHOTON_DARK, 255), 5)
        draw.pieslice((cx - 56, cy - 56, cx + 56, cy + 56), start=330, end=30, outline=(*PHOTON_DARK, 190), width=4)
        draw.line((cx, cy, cx + 54, cy - 17), fill=(*PHOTON_DARK, 175), width=3)
        draw.line((cx, cy, cx + 54, cy + 17), fill=(*PHOTON_DARK, 175), width=3)
    else:
        for r, color in ((49, TEAL), (34, PHOTON), (18, SPHENIX_BLUE)):
            draw.ellipse((cx - r, cy - r, cx + r, cy + r), outline=(*color, 230), width=5)
        draw.line((cx, cy, cx + 48, cy - 22), fill=(*PHOTON_DARK, 240), width=5)
        draw.ellipse((cx + 40, cy - 30, cx + 58, cy - 12), fill=(*PHOTON, 230))


def draw_segmented_ring(
    draw: ImageDraw.ImageDraw,
    center: tuple[int, int],
    r_inner: int,
    r_outer: int,
    fill: tuple[int, int, int],
    outline: tuple[int, int, int],
    start_offset: int = 0,
    highlight: tuple[int, int] | None = None,
) -> None:
    cx, cy = center
    for i in range(16):
        start = start_offset + i * 22
        end = start + 16
        color = fill
        alpha = 72
        if highlight is not None and highlight[0] <= i <= highlight[1]:
            color = PHOTON
            alpha = 150
        box_outer = (cx - r_outer, cy - r_outer, cx + r_outer, cy + r_outer)
        draw.pieslice(box_outer, start=start, end=end, fill=(*color, alpha))
    draw.ellipse((cx - r_outer, cy - r_outer, cx + r_outer, cy + r_outer), outline=(*outline, 210), width=4)
    draw.ellipse((cx - r_inner, cy - r_inner, cx + r_inner, cy + r_inner), outline=(*outline, 155), width=3)


def draw_context_visual(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    panel = (126, 332, 1434, 1128)
    draw.rounded_rectangle(panel, radius=10, fill=(*PANEL, 255), outline=(*PANEL_EDGE, 255), width=2)

    # Subtle RHIC program cue.
    draw.ellipse((230, 430, 1240, 1014), outline=(200, 214, 226, 150), width=6)
    draw.arc((230, 430, 1240, 1014), start=202, end=338, fill=(*SPHENIX_BLUE, 155), width=8)
    draw.arc((230, 430, 1240, 1014), start=24, end=158, fill=(*PHOTON_DARK, 140), width=8)
    draw.text((240, 388), "RHIC p+p at 200 GeV", font=font(TIMES_ITALIC, 34), fill=MUTED)
    draw.text((1018, 976), "hard-probes reference", font=font(TIMES_ITALIC, 31), fill=LIGHT_MUTED)

    center = (704, 690)
    draw.line((312, center[1], 1096, center[1]), fill=(166, 191, 211, 255), width=6)
    draw.ellipse((center[0] - 19, center[1] - 19, center[0] + 19, center[1] + 19), fill=(*SPHENIX_BLUE, 230))
    draw.text((336, center[1] - 52), "p", font=font(TIMES_BOLD, 38), fill=BLUE)
    draw.text((1036, center[1] - 52), "p", font=font(TIMES_BOLD, 38), fill=BLUE)

    draw_segmented_ring(draw, center, 82, 120, (126, 169, 198), (83, 124, 153), start_offset=10)
    draw_segmented_ring(draw, center, 148, 210, (245, 181, 34), (183, 137, 39), start_offset=2, highlight=(1, 3))
    draw_segmented_ring(draw, center, 230, 302, (18, 112, 142), (18, 96, 124), start_offset=12)
    draw.ellipse((center[0] - 330, center[1] - 330, center[0] + 330, center[1] + 330), outline=(57, 83, 112, 130), width=4)
    draw.text((center[0] - 98, center[1] - 22), "sPHENIX", font=font(TIMES_BOLD, 39), fill=BLUE)

    photon_start = (center[0] + 28, center[1] - 22)
    photon_end = (center[0] + 205, center[1] - 112)
    cone_left = (center[0] + 297, center[1] - 190)
    cone_right = (center[0] + 294, center[1] - 42)
    alpha_layer(base, lambda d: d.polygon([photon_start, cone_left, cone_right], fill=(245, 181, 34, 26)))
    draw.line((photon_start, cone_left), fill=(219, 157, 32, 90), width=3)
    draw.line((photon_start, cone_right), fill=(219, 157, 32, 90), width=3)
    wave = feynman_points(photon_start, photon_end, 13, 4.8, 180)
    glow = Image.new("RGBA", base.size, (0, 0, 0, 0))
    gdraw = ImageDraw.Draw(glow, "RGBA")
    draw_polyline(gdraw, wave, (*PHOTON, 115), 16)
    glow = glow.filter(ImageFilter.GaussianBlur(6))
    base.alpha_composite(glow)
    draw = ImageDraw.Draw(base, "RGBA")
    draw_polyline(draw, wave, (*PHOTON, 255), 8)
    draw_polyline(draw, wave, (*PHOTON_DARK, 210), 3)
    draw.ellipse((photon_end[0] - 19, photon_end[1] - 19, photon_end[0] + 19, photon_end[1] + 19), fill=(*PHOTON, 240), outline=(*PHOTON_DARK, 230), width=3)

    for angle, length, width, alpha in ((38, 166, 7, 105), (318, 142, 6, 85), (285, 160, 5, 75)):
        ex = center[0] + math.cos(math.radians(angle)) * length
        ey = center[1] + math.sin(math.radians(angle)) * length
        draw.line((center[0], center[1], ex, ey), fill=(*TEAL_SOFT, alpha), width=width)
        draw.ellipse((ex - 9, ey - 9, ex + 9, ey + 9), fill=(*TEAL_SOFT, alpha + 25))

    legend_x, legend_y = 1048, 500
    legend = [
        ("Tracking", (126, 169, 198)),
        ("EMCal", PHOTON),
        ("HCal", TEAL),
    ]
    for i, (label, color) in enumerate(legend):
        yy = legend_y + i * 58
        draw.rounded_rectangle((legend_x, yy, legend_x + 34, yy + 22), radius=4, fill=(*color, 190))
        draw.text((legend_x + 50, yy - 6), label, font=font(TIMES, 34), fill=MUTED)
    draw.text((898, 454), "isolated photon candidate", font=font(TIMES_ITALIC, 34), fill=BLUE)

    # Measurement chain: the whole talk in one quiet line.
    chain_y = 1032
    steps = [
        ("p+p", "reference"),
        ("EMCal", "clusters"),
        ("isolated", "photons"),
        ("corrected", "cross section"),
    ]
    x = 202
    for i, (top, bottom) in enumerate(steps):
        w = (190, 220, 230, 252)[i]
        draw.rounded_rectangle((x, chain_y, x + w, chain_y + 72), radius=8, fill=(255, 255, 255, 255), outline=(216, 225, 234, 255), width=2)
        top_w, _ = text_box(draw, top, font(TIMES_BOLD, 29))
        bottom_w, _ = text_box(draw, bottom, font(TIMES, 25))
        draw.text((x + (w - top_w) / 2, chain_y + 9), top, font=font(TIMES_BOLD, 29), fill=INK)
        draw.text((x + (w - bottom_w) / 2, chain_y + 42), bottom, font=font(TIMES, 25), fill=MUTED)
        if i < len(steps) - 1:
            draw.line((x + w + 18, chain_y + 36, x + w + 72, chain_y + 36), fill=(158, 178, 194, 255), width=3)
            draw.polygon([(x + w + 72, chain_y + 36), (x + w + 58, chain_y + 27), (x + w + 58, chain_y + 45)], fill=(158, 178, 194, 255))
        x += w + 78


def draw_context_cards(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x, y = 1510, 338
    w, h, gap = 870, 224, 32
    for idx, (heading, body) in enumerate(CONTEXT_REASONS):
        top = y + idx * (h + gap)
        stripe = SPHENIX_BLUE if idx == 0 else PHOTON if idx == 1 else TEAL
        draw.rounded_rectangle((x, top, x + w, top + h), radius=10, fill=(*CARD, 255), outline=(*CARD_EDGE, 255), width=2)
        draw.rounded_rectangle((x, top, x + 14, top + h), radius=6, fill=(*stripe, 255))
        draw_context_icon(draw, idx, (x + 86, top + 112))
        draw.text((x + 162, top + 38), heading, font=font(TIMES_BOLD, 41), fill=INK)
        draw_wrapped(draw, body, (x + 164, top + 98), 632, font(TIMES, 32), fill=MUTED, line_gap=10)


def draw_conceptual_visual(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    panel = (126, 333, 1422, 1112)
    draw.rounded_rectangle(panel, radius=10, fill=(*PANEL, 255), outline=(*PANEL_EDGE, 255), width=2)

    # A quiet RHIC-beam cue leading into one hard scattering.
    y0 = 704
    draw.line((212, y0, 720, y0), fill=(178, 199, 216, 255), width=6)
    draw.line((326, y0 - 62, 435, y0 - 14), fill=(178, 199, 216, 150), width=4)
    draw.line((326, y0 + 62, 435, y0 + 14), fill=(178, 199, 216, 150), width=4)
    draw.ellipse((268, y0 - 48, 364, y0 + 48), fill=(255, 255, 255, 255), outline=(98, 139, 169, 255), width=4)
    draw.ellipse((468, y0 - 48, 564, y0 + 48), fill=(255, 255, 255, 255), outline=(98, 139, 169, 255), width=4)
    draw.text((290, y0 - 19), "p", font=font(TIMES_BOLD, 42), fill=BLUE)
    draw.text((491, y0 - 19), "p", font=font(TIMES_BOLD, 42), fill=BLUE)

    collision = (640, 704)
    for r, alpha in ((70, 36), (46, 58), (24, 96)):
        draw.ellipse((collision[0] - r, collision[1] - r, collision[0] + r, collision[1] + r), fill=(*SPHENIX_BLUE, alpha))
    draw.ellipse((collision[0] - 12, collision[1] - 12, collision[0] + 12, collision[1] + 12), fill=(*SPHENIX_BLUE, 210))
    draw.text((520, 806), "hard scattering", font=font(TIMES_ITALIC, 34), fill=MUTED)

    # Faint colored recoil: present enough for the hard process, not enough to
    # make the opening slide feel like a photon+jet talk.
    for angle, length, width, alpha in ((54, 230, 9, 125), (77, 185, 6, 90), (34, 170, 5, 80)):
        ex = collision[0] + math.cos(math.radians(angle)) * length
        ey = collision[1] + math.sin(math.radians(angle)) * length
        draw.line((collision[0], collision[1], ex, ey), fill=(*TEAL_SOFT, alpha), width=width)
        draw.ellipse((ex - 12, ey - 12, ex + 12, ey + 12), fill=(*TEAL_SOFT, alpha + 35))
    draw.text((771, 914), "recoil activity", font=font(TIMES_ITALIC, 30), fill=(102, 132, 148))

    photon_start = (682, 668)
    impact = (1210, 520)
    dx, dy = impact[0] - photon_start[0], impact[1] - photon_start[1]
    length = math.hypot(dx, dy)
    nx, ny = -dy / length, dx / length
    cone_left = (impact[0] + nx * 145, impact[1] + ny * 145)
    cone_right = (impact[0] - nx * 145, impact[1] - ny * 145)
    alpha_layer(
        base,
        lambda d: d.polygon([photon_start, cone_left, cone_right], fill=(245, 181, 34, 30)),
    )
    draw.line((photon_start, cone_left), fill=(219, 157, 32, 90), width=3)
    draw.line((photon_start, cone_right), fill=(219, 157, 32, 90), width=3)

    wave = feynman_points(photon_start, impact, 17, 8.8, 260)
    glow = Image.new("RGBA", base.size, (0, 0, 0, 0))
    gdraw = ImageDraw.Draw(glow, "RGBA")
    draw_polyline(gdraw, wave, (245, 181, 34, 120), 17)
    glow = glow.filter(ImageFilter.GaussianBlur(7))
    base.alpha_composite(glow)
    draw = ImageDraw.Draw(base, "RGBA")
    draw_polyline(draw, wave, (*PHOTON, 255), 9)
    draw_polyline(draw, wave, (*PHOTON_DARK, 210), 3)

    draw_detector_face(base, impact)
    draw = ImageDraw.Draw(base, "RGBA")
    draw.text((826, 475), "isolated prompt photon", font=font(TIMES_ITALIC, 36), fill=BLUE)
    draw.text((938, 591), "quiet cone", font=font(TIMES_ITALIC, 29), fill=(145, 111, 43))


def draw_reason_cards(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x, y = 1510, 345
    w, h, gap = 870, 216, 34
    for idx, (heading, body) in enumerate(REASONS):
        top = y + idx * (h + gap)
        stripe = SPHENIX_BLUE if idx == 0 else TEAL if idx == 2 else PHOTON
        draw.rounded_rectangle((x, top, x + w, top + h), radius=10, fill=(*CARD, 255), outline=(*CARD_EDGE, 255), width=2)
        draw.rounded_rectangle((x, top, x + 14, top + h), radius=6, fill=(*stripe, 255))
        draw_icon(draw, idx, (x + 86, top + 108))
        draw.text((x + 162, top + 41), heading, font=font(TIMES_BOLD, 43), fill=INK)
        draw_wrapped(draw, body, (x + 164, top + 101), 620, font(TIMES, 34), fill=MUTED, line_gap=11)


def render_why_photons(output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    img = Image.new("RGBA", (W, H), (*SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")

    draw.rectangle((0, 0, W, H), fill=(*SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*PHOTON, 255))

    draw.text((132, 98), TITLE, font=font(TIMES_BOLD, 90), fill=INK)
    draw.text((136, 205), SUBTITLE, font=font(TIMES_ITALIC, 43), fill=BLUE)
    draw.line((132, 292, W - 132, 292), fill=(221, 226, 232), width=3)

    draw_conceptual_visual(img)
    draw_reason_cards(img)

    draw.rounded_rectangle((132, 1192, W - 132, 1282), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw.text((174, 1216), BRIDGE, font=font(TIMES_ITALIC, 40), fill=BLUE)

    png = output_dir / "hp2026_slide02_why_isolated_prompt_photons.png"
    img.convert("RGB").save(png, "PNG")

    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "output": str(png.relative_to(ROOT)),
        "size": [W, H],
        "mode": "RGB",
        "title": TITLE,
        "slide_role": "first content slide after the title slide",
        "core_claim": "Isolated prompt photons provide a clean p+p hard-scattering baseline at RHIC.",
        "source_basis": [
            "HP2026 title-slide visual language: Times New Roman, restrained institutional palette",
            "HP2026 photon reference deck linearization: photon-only opening arc",
            "Conceptual motivation pattern from Yeonju/Hanpu/Stefan references without copied figures",
        ],
        "google_slides_mutation": False,
        "notes": [
            "No slide number or provenance footer.",
            "No copied public/internal plots or figures.",
            "Photon-only motivation; BDT, purity, unfolding, and detector details are left for later slides.",
        ],
    }
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return png


def render_sphenix_hard_probe_context(output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    img = Image.new("RGBA", (W, H), (*SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")

    draw.rectangle((0, 0, W, H), fill=(*SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*PHOTON, 255))

    logo = load_sphenix_logo()
    if logo is not None:
        paste_fit(img, logo, (2062, 82, 2380, 168), anchor="right")

    draw.text((132, 88), CONTEXT_TITLE, font=font(TIMES_BOLD, 76), fill=INK)
    draw.text((136, 192), CONTEXT_SUBTITLE, font=font(TIMES_ITALIC, 40), fill=BLUE)
    draw.line((132, 292, W - 132, 292), fill=(221, 226, 232), width=3)

    draw_context_visual(img)
    draw_context_cards(img)

    draw.rounded_rectangle((132, 1190, W - 132, 1284), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw.text((174, 1217), CONTEXT_BRIDGE, font=font(TIMES_ITALIC, 39), fill=BLUE)

    png = output_dir / "hp2026_slide02_sphenix_hard_probes_photon_baseline.png"
    img.convert("RGB").save(png, "PNG")

    manifest_path = output_dir / "manifest.json"
    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "primary_output": str(png.relative_to(ROOT)),
        "size": [W, H],
        "mode": "RGB",
        "title": CONTEXT_TITLE,
        "slide_role": "first content slide after the title slide",
        "core_claim": "The p+p isolated prompt-photon result is the clean baseline entry point for the sPHENIX hard-probes program.",
        "source_basis": [
            "Reference inventory: HPslides_v1 and current PPG12 material are truth base; Yeonju/Hanpu/Stefan guide pacing.",
            "Deck linearization: establish sPHENIX/data context before photon definitions and analysis details.",
            "Conceptual sPHENIX/RHIC schematic drawn locally; no copied detector figure or result plot.",
        ],
        "google_slides_mutation": False,
        "candidate_files": [
            str(png.relative_to(ROOT)),
            str((output_dir / "hp2026_slide02_why_isolated_prompt_photons.png").relative_to(ROOT)),
        ],
        "notes": [
            "No slide number or provenance footer.",
            "Small sPHENIX logo uses cached official blue logo from the title-slide asset set.",
            "Detector/program visual is a conceptual schematic, not an official detector rendering.",
            "BDT, purity, unfolding, and detailed detector claims are reserved for later slides.",
        ],
    }
    with manifest_path.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return png


def render(output_dir: Path, variant: str) -> Path:
    if variant == "why-photons":
        return render_why_photons(output_dir)
    if variant == "sphenix-hard-probes":
        return render_sphenix_hard_probe_context(output_dir)
    raise ValueError(f"unknown variant: {variant}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--variant",
        choices=("sphenix-hard-probes", "why-photons"),
        default="sphenix-hard-probes",
        help="Slide 2 concept to render.",
    )
    args = parser.parse_args()
    png = render(args.output_dir, args.variant)
    print(png)
    print(args.output_dir / "manifest.json")


if __name__ == "__main__":
    main()

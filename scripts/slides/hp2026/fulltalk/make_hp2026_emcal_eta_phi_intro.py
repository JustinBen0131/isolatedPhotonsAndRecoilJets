#!/usr/bin/env python3
"""Local EMCal eta-phi intro slide for HP2026.

This is a local PNG prototype only. It does not mutate Google Slides.
"""

from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw, ImageEnhance, ImageFilter, ImageFont, ImageOps

import make_hp2026_fulltalk_candidates as full


ROOT = full.ROOT
OUT_DIR = ROOT / "outputs/manual-20260610-emcal-eta-phi-intro"
OUT_PNG = OUT_DIR / "backup_emcal_geometry_tower_map_candidate.png"
OUT_META = OUT_DIR / "backup_emcal_geometry_tower_map_candidate_metadata.json"
OUT_SCRIPT = OUT_DIR / "backup_emcal_geometry_tower_map_candidate_script.md"
OUT_HEADER = OUT_DIR / "backup_emcal_geometry_tower_map_candidate.header.json"

W, H = full.W, full.H

DETECTOR_ASSET = (
    ROOT
    / "outputs/manual-20260601-hp2026-opening-slide/presentations/hp2026-opening-slide/assets/"
    / "sphenix_detector_clean_reference_user_20260602.png"
)
DETECTOR_RENDER = (
    ROOT
    / "outputs/manual-20260601-hp2026-opening-slide/presentations/hp2026-opening-slide/assets/"
    / "bnl_sphenix_detector_rendering.jpg"
)
SEGMENTATION_REF = OUT_DIR / "reference/referenced_emcal_unwrap_slide.png"


def f(path: Path, size: int) -> ImageFont.FreeTypeFont:
    return full.font(path, size)


def center_text(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    text: str,
    font: ImageFont.ImageFont,
    fill: tuple[int, int, int] = full.INK,
    *,
    stroke_width: int = 0,
    stroke_fill: tuple[int, int, int, int] | None = None,
) -> None:
    tw, th = full.text_box(draw, text, font)
    x = (box[0] + box[2] - tw) / 2
    y = (box[1] + box[3] - th) / 2
    draw.text((x, y), text, font=font, fill=fill, stroke_width=stroke_width, stroke_fill=stroke_fill)


def draw_wrapped_centered(
    draw: ImageDraw.ImageDraw,
    text: str,
    box: tuple[int, int, int, int],
    font: ImageFont.ImageFont,
    fill: tuple[int, int, int],
    *,
    line_gap: int = 8,
) -> None:
    words = text.split()
    lines: list[str] = []
    current = ""
    max_w = box[2] - box[0]
    for word in words:
        trial = word if not current else f"{current} {word}"
        if full.text_box(draw, trial, font)[0] <= max_w:
            current = trial
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)

    heights = [full.text_box(draw, line, font)[1] for line in lines]
    total_h = sum(heights) + line_gap * max(0, len(lines) - 1)
    y = (box[1] + box[3] - total_h) / 2
    for line, line_h in zip(lines, heights):
        tw, _ = full.text_box(draw, line, font)
        draw.text(((box[0] + box[2] - tw) / 2, y), line, font=font, fill=fill)
        y += line_h + line_gap


def rounded_card(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], accent: tuple[int, int, int]) -> None:
    shadow = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow, "RGBA")
    sx0, sy0, sx1, sy1 = box
    sd.rounded_rectangle((sx0 + 7, sy0 + 8, sx1 + 7, sy1 + 8), radius=12, fill=(20, 38, 55, 22))
    shadow = shadow.filter(ImageFilter.GaussianBlur(8))
    draw.bitmap((0, 0), shadow, fill=None)
    draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=3)
    draw.rounded_rectangle((box[0], box[1], box[0] + 13, box[3]), radius=7, fill=(*accent, 235))


def draw_stage_header(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    number: str,
    title: str,
    subtitle: str,
    accent: tuple[int, int, int],
) -> None:
    x0, y0, x1, _ = box
    cy = y0 + 68
    badge = (x0 + 38, cy - 30, x0 + 98, cy + 30)
    draw.ellipse(badge, fill=(255, 255, 255, 255), outline=(*accent, 230), width=4)
    center_text(draw, badge, number, f(full.TIMES_BOLD, 33), accent)
    draw.text((x0 + 122, y0 + 36), title, font=f(full.TIMES_BOLD, 42), fill=full.INK)
    draw.text((x0 + 122, y0 + 90), subtitle, font=f(full.TIMES_ITALIC, 27), fill=full.MUTED)
    draw.line((x0 + 38, y0 + 132, x1 - 38, y0 + 132), fill=(224, 231, 238, 255), width=2)


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], color: tuple[int, int, int], width: int = 6) -> None:
    draw.line((start[0], start[1], end[0], end[1]), fill=(*color, 220), width=width)
    ang = math.atan2(end[1] - start[1], end[0] - start[0])
    size = 22
    pts = [
        end,
        (int(end[0] - size * math.cos(ang - 0.42)), int(end[1] - size * math.sin(ang - 0.42))),
        (int(end[0] - size * math.cos(ang + 0.42)), int(end[1] - size * math.sin(ang + 0.42))),
    ]
    draw.polygon(pts, fill=(*color, 235))


def small_arrow(
    draw: ImageDraw.ImageDraw,
    start: tuple[int, int],
    end: tuple[int, int],
    color: tuple[int, int, int],
    width: int = 4,
    *,
    head: int = 13,
) -> None:
    draw.line((start[0], start[1], end[0], end[1]), fill=(*color, 230), width=width)
    ang = math.atan2(end[1] - start[1], end[0] - start[0])
    pts = [
        end,
        (int(end[0] - head * math.cos(ang - 0.48)), int(end[1] - head * math.sin(ang - 0.48))),
        (int(end[0] - head * math.cos(ang + 0.48)), int(end[1] - head * math.sin(ang + 0.48))),
    ]
    draw.polygon(pts, fill=(*color, 240))


def draw_detector_card(base: Image.Image, draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    draw_stage_header(draw, box, "1", "Cylindrical EMCal", "segmented calorimeter layer", full.PHOTON_DARK)
    x0, y0, x1, y1 = box
    if DETECTOR_ASSET.exists():
        det = Image.open(DETECTOR_ASSET).convert("RGBA")
        det = full.crop_visible(det, white_threshold=252, pad=6)
        det = full.fit(det, x1 - x0 - 130, y1 - y0 - 260)
        px = x0 + (x1 - x0 - det.width) // 2 + 20
        py = y0 + 174
        base.alpha_composite(det, (px, py))

        # Highlight the approximate barrel EMCal layer without relying on the
        # detector image labels.
        overlay = Image.new("RGBA", (W, H), (0, 0, 0, 0))
        od = ImageDraw.Draw(overlay, "RGBA")
        band = (
            px + int(det.width * 0.53),
            py + int(det.height * 0.40),
            px + int(det.width * 0.88),
            py + int(det.height * 0.73),
        )
        od.ellipse(band, outline=(*full.PHOTON_DARK, 230), width=8)
        od.ellipse(
            (band[0] + 18, band[1] + 18, band[2] - 18, band[3] - 18),
            outline=(*full.PHOTON, 180),
            width=8,
        )
        od.rectangle((band[0], band[1], band[2], band[1] + 72), fill=(*full.PHOTON, 34))
        base.alpha_composite(overlay)

        emcal_tag = (x1 - 238, y0 + 358, x1 - 94, y0 + 408)
        draw.rounded_rectangle(emcal_tag, radius=8, fill=(255, 248, 229, 238), outline=(*full.PHOTON_DARK, 210), width=2)
        center_text(draw, emcal_tag, "EMCal", f(full.TIMES_BOLD, 29), full.PHOTON_DARK)
        arrow(draw, (emcal_tag[0] + 12, emcal_tag[3] + 8), (px + int(det.width * 0.66), py + int(det.height * 0.58)), full.PHOTON_DARK, width=4)

        spec = (x0 + 86, y1 - 172, x1 - 86, y1 - 104)
        draw.rounded_rectangle(spec, radius=9, fill=(247, 251, 253, 255), outline=(210, 224, 234, 255), width=2)
        draw_wrapped_centered(
            draw,
            "tungsten / scintillating-fiber sampling EMCal",
            spec,
            f(full.TIMES_BOLD, 27),
            full.BLUE,
            line_gap=4,
        )

        label = (x0 + 86, y1 - 92, x1 - 86, y1 - 34)
        draw.rounded_rectangle(label, radius=9, fill=(255, 248, 229, 255), outline=(*full.PHOTON_DARK, 190), width=2)
        center_text(draw, label, "96 η × 256 φ towers", f(full.TIMES_BOLD, 32), full.BLUE)


def draw_unwrap_card(base: Image.Image, draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    draw_stage_header(draw, box, "2", "Unwrap the barrel", "same towers, flattened view", full.SPHENIX_BLUE)
    x0, y0, x1, y1 = box
    cx = (x0 + x1) // 2
    top = y0 + 185
    bottom = y0 + 495

    # Curved barrel strip.
    barrel = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    bd = ImageDraw.Draw(barrel, "RGBA")
    for i in range(17):
        t0 = i / 17
        t1 = (i + 1) / 17
        xl0 = x0 + 104 + int(426 * t0)
        xl1 = x0 + 104 + int(426 * t1)
        offset0 = int(70 * math.sin((t0 - 0.5) * math.pi))
        offset1 = int(70 * math.sin((t1 - 0.5) * math.pi))
        poly = [
            (xl0, top + offset0),
            (xl1, top + offset1),
            (xl1 + 18, bottom + offset1 // 2),
            (xl0 + 18, bottom + offset0 // 2),
        ]
        fill = (224, 244, 250, 255) if i % 2 == 0 else (238, 248, 251, 255)
        bd.polygon(poly, fill=fill)
        bd.line(poly + [poly[0]], fill=(*full.TEAL, 160), width=2)
    bd.arc((x0 + 96, top - 62, x1 - 96, top + 112), start=8, end=172, fill=(*full.TEAL, 230), width=7)
    bd.arc((x0 + 112, bottom - 68, x1 - 78, bottom + 70), start=8, end=172, fill=(*full.TEAL, 165), width=5)
    base.alpha_composite(barrel)
    cyl_label = (x0 + 155, y0 + 540, x1 - 155, y0 + 590)
    draw.rounded_rectangle(cyl_label, radius=8, fill=(255, 255, 255, 220), outline=(207, 223, 232, 255), width=2)
    center_text(draw, cyl_label, "cylindrical tower layer", f(full.TIMES_BOLD, 30), full.TEAL)

    # Unwrap direction.
    arrow(draw, (cx - 92, y0 + 660), (cx + 92, y0 + 660), full.PHOTON_DARK, width=7)
    center_text(draw, (cx - 86, y0 + 598, cx + 86, y0 + 646), "unwrap", f(full.TIMES_BOLD, 32), full.PHOTON_DARK)

    # Flattened strip.
    gx0, gy0 = x0 + 104, y0 + 720
    gw, gh = x1 - x0 - 208, 146
    draw.rounded_rectangle((gx0, gy0, gx0 + gw, gy0 + gh), radius=8, fill=(247, 251, 253, 255), outline=(*full.SPHENIX_BLUE, 180), width=3)
    cols = 18
    rows = 4
    for c in range(1, cols):
        x = gx0 + c * gw / cols
        draw.line((x, gy0, x, gy0 + gh), fill=(210, 223, 232, 255), width=2)
    for r in range(1, rows):
        y = gy0 + r * gh / rows
        draw.line((gx0, y, gx0 + gw, y), fill=(210, 223, 232, 255), width=2)
    strip_label = (gx0 + 80, gy0 + 48, gx0 + gw - 80, gy0 + 98)
    draw.rounded_rectangle(strip_label, radius=8, fill=(255, 255, 255, 225), outline=(205, 222, 232, 255), width=2)
    center_text(draw, strip_label, "flat indexed tower plane", f(full.TIMES_BOLD, 30), full.BLUE)


def tower_color(val: float) -> tuple[int, int, int, int]:
    cold = (244, 248, 250)
    warm = full.PHOTON
    t = max(0.0, min(1.0, val))
    r = int(cold[0] * (1 - t) + warm[0] * t)
    g = int(cold[1] * (1 - t) + warm[1] * t)
    b = int(cold[2] * (1 - t) + warm[2] * t)
    return (r, g, b, 255)


def draw_eta_phi_card(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    draw_stage_header(draw, box, "3", "Local η–φ tower map", "zoom from the 96 × 256 tower plane", full.TEAL)
    x0, y0, x1, y1 = box

    gx0, gy0 = x0 + 88, y0 + 246
    cell = 34
    cols, rows = 18, 13
    gx1, gy1 = gx0 + cols * cell, gy0 + rows * cell
    draw.rounded_rectangle((gx0 - 26, gy0 - 26, gx1 + 44, gy1 + 64), radius=12, fill=(249, 252, 253, 255), outline=(217, 226, 235, 255), width=3)

    seed_c, seed_r = 9, 6
    for r in range(rows):
        for c in range(cols):
            dist2 = ((c - seed_c) / 1.35) ** 2 + ((r - seed_r) / 1.55) ** 2
            val = math.exp(-0.55 * dist2)
            x = gx0 + c * cell
            y = gy0 + r * cell
            draw.rectangle((x, y, x + cell, y + cell), fill=tower_color(val * 0.95), outline=(213, 224, 232, 255), width=2)

    # Local window and seed.
    win = (gx0 + 6 * cell, gy0 + 3 * cell, gx0 + 13 * cell, gy0 + 10 * cell)
    draw.rectangle(win, outline=(*full.SPHENIX_BLUE, 255), width=6)
    draw.ellipse(
        (
            gx0 + seed_c * cell + cell / 2 - 10,
            gy0 + seed_r * cell + cell / 2 - 10,
            gx0 + seed_c * cell + cell / 2 + 10,
            gy0 + seed_r * cell + cell / 2 + 10,
        ),
        fill=(0, 0, 0, 255),
    )
    seed_label = (gx0 + 13 * cell + 12, gy0 + 6 * cell - 20, gx0 + 17 * cell + 20, gy0 + 6 * cell + 30)
    draw.rounded_rectangle(seed_label, radius=7, fill=(255, 255, 255, 225), outline=(203, 216, 227, 255), width=2)
    center_text(draw, seed_label, "seed", f(full.TIMES_ITALIC, 25), full.MUTED)

    # Axes.
    arrow(draw, (gx0 - 8, gy1 + 36), (gx1 + 2, gy1 + 36), full.SPHENIX_BLUE, width=5)
    arrow(draw, (gx0 - 36, gy1 + 8), (gx0 - 36, gy0 + 4), full.TEAL, width=5)
    draw.text((gx1 + 20, gy1 + 15), "η", font=f(full.TIMES_ITALIC, 42), fill=full.SPHENIX_BLUE)
    draw.text((gx0 - 74, gy0 - 44), "φ", font=f(full.TIMES_ITALIC, 42), fill=full.TEAL)

    local_label = (gx0 + 6 * cell - 12, gy0 + 2 * cell - 54, gx0 + 13 * cell + 12, gy0 + 2 * cell - 8)
    draw.rounded_rectangle(local_label, radius=7, fill=(255, 255, 255, 235), outline=(*full.SPHENIX_BLUE, 190), width=2)
    center_text(draw, local_label, "local cluster window", f(full.TIMES_BOLD, 28), full.SPHENIX_BLUE)

    caption = (
        x0 + 78,
        y1 - 126,
        x1 - 78,
        y1 - 34,
    )
    draw.rounded_rectangle(caption, radius=8, fill=(239, 246, 250, 255), outline=(211, 225, 235, 255), width=2)
    draw_wrapped_centered(
        draw,
        "Shower-shape variables summarize how the energy is arranged in this local tower window.",
        caption,
        f(full.TIMES_BOLD, 30),
        full.BLUE,
        line_gap=6,
    )


def draw_fact_chip(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    label: str,
    value: str,
    accent: tuple[int, int, int],
) -> None:
    draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 245), outline=(*accent, 150), width=2)
    draw.text((box[0] + 22, box[1] + 15), label, font=f(full.TIMES_BOLD, 25), fill=accent)
    draw.text((box[0] + 22, box[1] + 47), value, font=f(full.TIMES_BOLD, 29), fill=full.BLUE)


def draw_wrapped_left(
    draw: ImageDraw.ImageDraw,
    text: str,
    xy: tuple[int, int],
    max_w: int,
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
        if full.text_box(draw, trial, font)[0] <= max_w:
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
        y += full.text_box(draw, line, font)[1] + line_gap
    return y


def paste_rounded_image(
    base: Image.Image,
    image: Image.Image,
    box: tuple[int, int, int, int],
    *,
    radius: int = 14,
    border: tuple[int, int, int] = full.PANEL_EDGE,
    border_width: int = 3,
    fit_mode: str = "contain",
    fill: tuple[int, int, int] = (255, 255, 255),
) -> tuple[int, int, int, int]:
    x0, y0, x1, y1 = box
    w, h = x1 - x0, y1 - y0
    image = image.convert("RGBA")
    if fit_mode == "cover":
        fitted = ImageOps.fit(image, (w, h), method=Image.Resampling.LANCZOS, centering=(0.5, 0.5))
    else:
        fitted = Image.new("RGBA", (w, h), (*fill, 255))
        tmp = full.fit(image, w, h)
        fitted.alpha_composite(tmp, ((w - tmp.width) // 2, (h - tmp.height) // 2))

    mask = Image.new("L", (w, h), 0)
    md = ImageDraw.Draw(mask)
    md.rounded_rectangle((0, 0, w, h), radius=radius, fill=255)
    base.paste(fitted, (x0, y0), mask)
    d = ImageDraw.Draw(base, "RGBA")
    d.rounded_rectangle(box, radius=radius, outline=(*border, 255), width=border_width)
    return box


def paste_detector_with_emcal_highlight(
    base: Image.Image,
    draw: ImageDraw.ImageDraw,
    image: Image.Image,
    box: tuple[int, int, int, int],
) -> None:
    paste_rounded_image(base, image, box, radius=12, border=(210, 224, 234), fit_mode="cover")
    x0, y0, x1, y1 = box

    overlay = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    od = ImageDraw.Draw(overlay, "RGBA")
    # Dim the full detector just enough that the EMCal layer becomes the first
    # read without destroying the source rendering.
    od.rounded_rectangle((x0, y0, x1, y1), radius=12, fill=(255, 255, 255, 58))

    # The EMCal is the projective barrel layer just outside the tracking
    # volume. Highlight the orange barrel skin and its visible annular face,
    # instead of drawing an oversized detector-wide oval.
    ring = (
        x0 + int((x1 - x0) * 0.515),
        y0 + int((y1 - y0) * 0.455),
        x0 + int((x1 - x0) * 0.735),
        y0 + int((y1 - y0) * 0.785),
    )
    for width, alpha in [(24, 58), (15, 122)]:
        od.ellipse(ring, outline=(*full.PHOTON, alpha), width=width)
    od.ellipse(ring, outline=(*full.PHOTON_DARK, 245), width=7)
    inner = (ring[0] + 20, ring[1] + 20, ring[2] - 20, ring[3] - 20)
    od.ellipse(inner, outline=(*full.PHOTON_DARK, 170), width=4)
    # A longitudinal band points to the same barrel skin, not only the face.
    band = [
        (x0 + int((x1 - x0) * 0.315), y0 + int((y1 - y0) * 0.455)),
        (x0 + int((x1 - x0) * 0.565), y0 + int((y1 - y0) * 0.505)),
        (x0 + int((x1 - x0) * 0.565), y0 + int((y1 - y0) * 0.575)),
        (x0 + int((x1 - x0) * 0.315), y0 + int((y1 - y0) * 0.525)),
    ]
    od.polygon(band, fill=(*full.PHOTON, 76))
    od.line(band + [band[0]], fill=(*full.PHOTON_DARK, 220), width=5, joint="curve")

    label = (x0 + 342, y0 + 22, x0 + 592, y0 + 80)
    od.rounded_rectangle(label, radius=9, fill=(255, 248, 229, 245), outline=(*full.PHOTON_DARK, 215), width=2)
    base.alpha_composite(overlay)
    center_text(draw, label, "EMCal layer", f(full.TIMES_BOLD, 31), full.PHOTON_DARK)
    arrow(draw, (label[0] + 34, label[3] + 8), (ring[0] + 34, ring[1] + 96), full.PHOTON_DARK, width=4)


def draw_real_card(
    base: Image.Image,
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    accent: tuple[int, int, int],
    title: str,
    subtitle: str,
) -> None:
    rounded_card(draw, box, accent)
    x0, y0, x1, _ = box
    draw.text((x0 + 44, y0 + 27), title, font=f(full.TIMES_BOLD, 47), fill=full.INK)
    draw_wrapped_left(
        draw,
        subtitle,
        (x0 + 44, y0 + 82),
        x1 - x0 - 88,
        f(full.TIMES_ITALIC, 30),
        full.MUTED,
        line_gap=4,
    )
    draw.line((x0 + 44, y0 + 126, x1 - 44, y0 + 126), fill=(218, 229, 236, 255), width=2)


def load_detector_anchor() -> Image.Image:
    source = DETECTOR_RENDER if DETECTOR_RENDER.exists() else DETECTOR_ASSET
    img = Image.open(source).convert("RGBA")
    img = ImageEnhance.Contrast(img).enhance(1.04)
    img = ImageEnhance.Sharpness(img).enhance(1.15)
    return img


def load_segmentation_crop() -> Image.Image:
    img = Image.open(SEGMENTATION_REF).convert("RGBA")
    # Crop the visual barrel -> sector -> module/tower strip, not the dense
    # explanatory text above it.
    crop = img.crop((0, 465, 1548, 850))
    crop = ImageEnhance.Contrast(crop).enhance(1.05)
    crop = ImageEnhance.Sharpness(crop).enhance(1.2)
    return crop


def load_tower_block_crop() -> Image.Image:
    img = Image.open(SEGMENTATION_REF).convert("RGBA")
    crop = img.crop((1115, 55, 1575, 390))
    crop = ImageEnhance.Contrast(crop).enhance(1.04)
    crop = ImageEnhance.Sharpness(crop).enhance(1.15)
    return crop


def draw_sector_module_schematic(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    """Code-native sector sketch: one barrel sector subdivided into modules."""
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(
        box,
        radius=12,
        fill=(248, 252, 254, 255),
        outline=(207, 222, 232, 255),
        width=2,
    )
    center_text(
        draw,
        (x0 + 18, y0 + 12, x1 - 18, y0 + 56),
        "barrel sector → 24 modules along φ",
        f(full.TIMES_BOLD, 28),
        full.INK,
    )

    # Annular detector cross-section cue with EMCal as a highlighted outer
    # calorimeter segment, rather than a generic pie slice.
    cx, cy = x0 + 108, y0 + 144
    outer = (cx - 78, cy - 78, cx + 78, cy + 78)
    mid = (cx - 58, cy - 58, cx + 58, cy + 58)
    inner = (cx - 36, cy - 36, cx + 36, cy + 36)
    draw.ellipse(outer, fill=(246, 250, 252, 255), outline=(56, 76, 87, 255), width=4)
    draw.ellipse(mid, outline=(164, 183, 194, 255), width=4)
    draw.ellipse(inner, fill=(255, 255, 255, 255), outline=(190, 202, 210, 255), width=3)
    for deg in range(0, 360, 15):
        draw.arc(outer, deg + 2, deg + 9, fill=(*full.TEAL, 190), width=5)
    for deg in range(-20, 22, 7):
        draw.line(
            (
                cx + 42 * math.cos(math.radians(deg)),
                cy + 42 * math.sin(math.radians(deg)),
                cx + 78 * math.cos(math.radians(deg)),
                cy + 78 * math.sin(math.radians(deg)),
            ),
            fill=(*full.PHOTON_DARK, 160),
            width=2,
        )
    draw.arc(outer, -20, 20, fill=(*full.PHOTON_DARK, 255), width=13)
    draw.arc((cx - 67, cy - 67, cx + 67, cy + 67), -20, 20, fill=(*full.PHOTON, 230), width=9)
    center_text(draw, (cx - 58, cy + 84, cx + 58, cy + 120), "CEMC", f(full.TIMES_BOLD, 23), full.TEAL)

    arrow(draw, (x0 + 200, y0 + 146), (x0 + 286, y0 + 146), full.PHOTON_DARK, width=5)

    # Unrolled projective sector strip: 24 tapered module cells with a visible
    # end-cap and one selected module. This is a schematic, not a data plot.
    sx0, sy0 = x0 + 314, y0 + 76
    sx1, sy1 = x1 - 38, y0 + 216
    sector = [
        (sx0, sy0 + 26),
        (sx1, sy0 + 2),
        (sx1 - 30, sy1 - 8),
        (sx0 + 24, sy1 - 28),
    ]
    shadow = [(p[0] + 4, p[1] + 5) for p in sector]
    draw.polygon(shadow, fill=(25, 52, 68, 18))
    draw.polygon(sector, fill=(238, 250, 253, 255))
    draw.line(sector + [sector[0]], fill=(*full.TEAL, 245), width=4, joint="curve")
    n_modules = 24
    top0, top1, bot0, bot1 = sector[0], sector[1], sector[3], sector[2]
    for i in range(1, n_modules):
        t = i / n_modules
        xa = top0[0] * (1 - t) + top1[0] * t
        ya = top0[1] * (1 - t) + top1[1] * t
        xb = bot0[0] * (1 - t) + bot1[0] * t
        yb = bot0[1] * (1 - t) + bot1[1] * t
        seam = (111, 180, 199, 215) if i % 4 == 0 else (139, 199, 214, 185)
        draw.line((xa, ya, xb, yb), fill=seam, width=2)
    for j in range(1, 4):
        t = j / 4
        xa = top0[0] * (1 - t) + bot0[0] * t
        ya = top0[1] * (1 - t) + bot0[1] * t
        xb = top1[0] * (1 - t) + bot1[0] * t
        yb = top1[1] * (1 - t) + bot1[1] * t
        draw.line((xa, ya, xb, yb), fill=(200, 224, 231, 190), width=2)

    hi0 = 9 / n_modules
    hi1 = 10 / n_modules
    highlight = [
        (
            top0[0] * (1 - hi0) + top1[0] * hi0,
            top0[1] * (1 - hi0) + top1[1] * hi0,
        ),
        (
            top0[0] * (1 - hi1) + top1[0] * hi1,
            top0[1] * (1 - hi1) + top1[1] * hi1,
        ),
        (
            bot0[0] * (1 - hi1) + bot1[0] * hi1,
            bot0[1] * (1 - hi1) + bot1[1] * hi1,
        ),
        (
            bot0[0] * (1 - hi0) + bot1[0] * hi0,
            bot0[1] * (1 - hi0) + bot1[1] * hi0,
        ),
    ]
    draw.polygon(highlight, fill=(*full.PHOTON, 185), outline=(*full.PHOTON_DARK, 255))
    draw.line((sx0 + 16, sy1 - 20, sx1 - 18, sy1 - 2), fill=(255, 255, 255, 155), width=2)
    draw.line((sx1, sy0 + 2, sx1 - 30, sy1 - 8), fill=(42, 83, 96, 210), width=5)
    module_label = (sx0 + 112, sy1 + 16, sx0 + 308, sy1 + 58)
    draw.rounded_rectangle(module_label, radius=8, fill=(255, 248, 229, 255), outline=(*full.PHOTON_DARK, 185), width=2)
    center_text(draw, module_label, "one module", f(full.TIMES_BOLD, 24), full.PHOTON_DARK)
    arrow(draw, (sx0 + 210, sy1 + 14), (sx0 + 184, sy1 - 30), full.PHOTON_DARK, width=4)


def draw_module_tower_diagram(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    """Clean room-readable redraw of the module/block/tower hierarchy."""
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(
        box,
        radius=12,
        fill=(255, 255, 255, 250),
        outline=(205, 222, 232, 255),
        width=2,
    )

    center_text(
        draw,
        (x0 + 16, y0 + 12, x1 - 16, y0 + 54),
        "one module = four 2×2 tower blocks",
        f(full.TIMES_BOLD, 28),
        full.INK,
    )

    block_gap = 14
    block_w = 84
    block_h = 92
    total_w = 4 * block_w + 3 * block_gap
    bx0 = x0 + (x1 - x0 - total_w) // 2
    by0 = y0 + 74

    for b in range(4):
        lx = bx0 + b * (block_w + block_gap)
        block = (lx, by0, lx + block_w, by0 + block_h)
        draw.rounded_rectangle(
            block,
            radius=7,
            fill=(250, 253, 255, 255),
            outline=(*full.SPHENIX_BLUE, 210),
            width=3,
        )
        cell_w = block_w // 2
        cell_h = block_h // 2
        for r in range(2):
            for c in range(2):
                xx = lx + c * cell_w
                yy = by0 + r * cell_h
                draw.rectangle(
                    (xx + 4, yy + 4, xx + cell_w - 4, yy + cell_h - 4),
                    fill=(235, 247, 252, 255),
                    outline=(190, 211, 223, 255),
                    width=2,
                )
        center_text(
            draw,
            (lx, by0 + block_h + 10, lx + block_w, by0 + block_h + 43),
            "2×2",
            f(full.TIMES_BOLD, 25),
            full.SPHENIX_BLUE,
        )

    brace_y = by0 + block_h + 56
    draw.line((bx0, brace_y, bx0 + total_w, brace_y), fill=(*full.SPHENIX_BLUE, 230), width=4)
    draw.line((bx0, brace_y - 10, bx0, brace_y + 10), fill=(*full.SPHENIX_BLUE, 230), width=4)
    draw.line((bx0 + total_w, brace_y - 10, bx0 + total_w, brace_y + 10), fill=(*full.SPHENIX_BLUE, 230), width=4)
    center_text(
        draw,
        (bx0 + 80, brace_y + 10, bx0 + total_w - 80, brace_y + 48),
        "16 towers",
        f(full.TIMES_BOLD, 27),
        full.BLUE,
    )

    arrow_y = by0 + block_h + 78
    arrow(draw, (bx0, arrow_y), (bx0 + total_w, arrow_y), full.SPHENIX_BLUE, width=4)
    draw.text((bx0 + total_w + 12, arrow_y - 21), "+φ", font=f(full.TIMES_ITALIC, 30), fill=full.SPHENIX_BLUE)
    arrow(draw, (bx0 - 30, by0 + block_h + 14), (bx0 - 30, by0 + 8), full.TEAL, width=4)
    draw.text((bx0 - 78, by0 + 12), "+η", font=f(full.TIMES_ITALIC, 30), fill=full.TEAL)


def draw_global_tower_plane(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    """Room-readable schematic of the full EMCal tower readout as an eta-phi plane."""
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(
        box,
        radius=12,
        fill=(248, 252, 254, 255),
        outline=(207, 222, 232, 255),
        width=2,
    )
    center_text(
        draw,
        (x0 + 18, y0 + 14, x1 - 18, y0 + 58),
        "cylindrical barrel → flat tower image",
        f(full.TIMES_BOLD, 29),
        full.INK,
    )

    # Left: cylindrical EMCal skin, drawn as a segmented annulus.
    cx, cy = x0 + 118, y0 + 160
    outer = (cx - 80, cy - 80, cx + 80, cy + 80)
    inner = (cx - 52, cy - 52, cx + 52, cy + 52)
    draw.ellipse(outer, fill=(246, 250, 252, 255), outline=(48, 70, 82, 255), width=4)
    draw.ellipse(inner, fill=(255, 255, 255, 255), outline=(169, 188, 198, 255), width=4)
    for deg in range(0, 360, 12):
        draw.arc(outer, deg + 1, deg + 7, fill=(*full.TEAL, 205), width=5)
    draw.arc(outer, 315, 45, fill=(*full.PHOTON_DARK, 255), width=13)
    draw.arc((cx - 68, cy - 68, cx + 68, cy + 68), 315, 45, fill=(*full.PHOTON, 230), width=9)
    center_text(draw, (cx - 70, cy + 88, cx + 70, cy + 126), "CEMC barrel", f(full.TIMES_BOLD, 23), full.TEAL)

    arrow(draw, (x0 + 214, y0 + 160), (x0 + 296, y0 + 160), full.PHOTON_DARK, width=5)

    # Right: global readout image. We draw a coarse version of the 96 x 256
    # tower plane, explicitly labeling which direction is eta and phi.
    gx0, gy0 = x0 + 324, y0 + 86
    gx1, gy1 = x1 - 44, y0 + 234
    draw.rounded_rectangle((gx0, gy0, gx1, gy1), radius=8, fill=(255, 255, 255, 255), outline=(*full.SPHENIX_BLUE, 190), width=3)
    n_eta, n_phi = 16, 8
    for i in range(1, n_eta):
        xx = gx0 + (gx1 - gx0) * i / n_eta
        draw.line((xx, gy0, xx, gy1), fill=(204, 220, 230, 255), width=2)
    for j in range(1, n_phi):
        yy = gy0 + (gy1 - gy0) * j / n_phi
        draw.line((gx0, yy, gx1, yy), fill=(204, 220, 230, 255), width=2)
    # A small highlighted analysis window foreshadows the right card.
    wx0 = gx0 + int((gx1 - gx0) * 0.58)
    wy0 = gy0 + int((gy1 - gy0) * 0.34)
    wx1 = gx0 + int((gx1 - gx0) * 0.78)
    wy1 = gy0 + int((gy1 - gy0) * 0.66)
    draw.rectangle((wx0, wy0, wx1, wy1), fill=(*full.PHOTON, 55), outline=(*full.SPHENIX_BLUE, 255), width=4)
    draw.ellipse(((wx0 + wx1) / 2 - 7, (wy0 + wy1) / 2 - 7, (wx0 + wx1) / 2 + 7, (wy0 + wy1) / 2 + 7), fill=(0, 0, 0, 255))
    draw.text((gx0 + 8, gy0 - 38), "φ index", font=f(full.TIMES_ITALIC, 27), fill=full.TEAL)
    draw.text((gx1 - 80, gy1 + 10), "η index", font=f(full.TIMES_ITALIC, 27), fill=full.SPHENIX_BLUE)
    draw.line((gx0, gy1 + 36, gx1, gy1 + 36), fill=(*full.SPHENIX_BLUE, 230), width=3)
    draw.line((gx0 - 30, gy0, gx0 - 30, gy1), fill=(*full.TEAL, 230), width=3)

    note = (x0 + 42, y0 + 286, x1 - 42, y0 + 372)
    draw.rounded_rectangle(note, radius=9, fill=(239, 247, 251, 255), outline=(207, 222, 232, 255), width=2)
    center_text(
        draw,
        (note[0] + 14, note[1] + 10, note[2] - 14, note[3] - 10),
        "96 η bins × 256 φ bins",
        f(full.TIMES_BOLD, 36),
        full.BLUE,
    )

    count = (x0 + 42, y0 + 404, x1 - 42, y0 + 514)
    draw.rounded_rectangle(count, radius=9, fill=(255, 248, 229, 255), outline=(*full.PHOTON_DARK, 150), width=2)
    center_text(
        draw,
        (count[0] + 14, count[1] + 10, count[2] - 14, count[1] + 56),
        "24,576 EMCal towers",
        f(full.TIMES_BOLD, 35),
        full.PHOTON_DARK,
    )
    center_text(
        draw,
        (count[0] + 14, count[1] + 56, count[2] - 14, count[3] - 10),
        "each tower: Δη × Δφ = 0.025 × 0.025",
        f(full.TIMES_BOLD, 26),
        full.BLUE,
    )

    bridge = (x0 + 42, y0 + 558, x1 - 42, y1 - 32)
    draw.rounded_rectangle(bridge, radius=9, fill=(255, 255, 255, 245), outline=(207, 222, 232, 255), width=2)
    draw_wrapped_centered(
        draw,
        "analysis zooms from the global tower plane to a local window around each photon candidate",
        bridge,
        f(full.TIMES_BOLD, 28),
        full.BLUE,
        line_gap=6,
    )


def draw_clean_tower_window(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    *,
    cell: int = 24,
    cols: int = 11,
    rows: int = 8,
) -> None:
    x0, y0, x1, y1 = box
    grid_w, grid_h = cols * cell, rows * cell
    gx0 = x0 + (x1 - x0 - grid_w) // 2
    gy0 = y0 + 12
    seed_c, seed_r = cols // 2, rows // 2
    draw.rounded_rectangle((x0, y0, x1, y1), radius=12, fill=(248, 252, 254, 255), outline=(211, 225, 235, 255), width=2)
    for r in range(rows):
        for c in range(cols):
            dist2 = ((c - seed_c) / 1.35) ** 2 + ((r - seed_r) / 1.45) ** 2
            val = math.exp(-0.58 * dist2)
            x = gx0 + c * cell
            y = gy0 + r * cell
            draw.rectangle((x, y, x + cell, y + cell), fill=tower_color(val * 0.98), outline=(208, 222, 232, 255), width=2)
    local = (gx0 + 3 * cell, gy0 + 2 * cell, gx0 + 8 * cell, gy0 + 6 * cell)
    draw.rectangle(local, outline=(*full.SPHENIX_BLUE, 255), width=6)
    draw.ellipse(
        (
            gx0 + seed_c * cell + cell / 2 - 10,
            gy0 + seed_r * cell + cell / 2 - 10,
            gx0 + seed_c * cell + cell / 2 + 10,
            gy0 + seed_r * cell + cell / 2 + 10,
        ),
        fill=(0, 0, 0, 255),
    )
    draw.text((gx0 + 8 * cell + 14, gy0 + 4 * cell - 15), "seed", font=f(full.TIMES_ITALIC, 24), fill=full.MUTED)
    arrow(draw, (gx0 - 6, gy0 + grid_h + 34), (gx0 + grid_w + 4, gy0 + grid_h + 34), full.SPHENIX_BLUE, width=5)
    arrow(draw, (gx0 - 34, gy0 + grid_h + 8), (gx0 - 34, gy0 + 6), full.TEAL, width=5)
    draw.text((gx0 + grid_w + 22, gy0 + grid_h + 12), "η", font=f(full.TIMES_ITALIC, 42), fill=full.SPHENIX_BLUE)
    draw.text((gx0 - 78, gy0 - 40), "φ", font=f(full.TIMES_ITALIC, 42), fill=full.TEAL)


def draw_stage_label(
    draw: ImageDraw.ImageDraw,
    center: tuple[int, int],
    number: str,
    text: str,
    accent: tuple[int, int, int],
    *,
    width: int = 390,
) -> None:
    x, y = center
    box = (x - width // 2, y - 36, x + width // 2, y + 36)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 248), outline=(*accent, 185), width=2)
    badge = (box[0] + 16, y - 24, box[0] + 64, y + 24)
    draw.ellipse(badge, fill=(255, 255, 255, 255), outline=(*accent, 220), width=3)
    center_text(draw, badge, number, f(full.TIMES_BOLD, 25), accent)
    draw.text((box[0] + 78, y - 20), text, font=f(full.TIMES_BOLD, 31), fill=full.INK)


def draw_vector_barrel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> tuple[int, int]:
    """Stylized EMCal barrel with visible segmentation and highlighted skin."""
    x0, y0, x1, y1 = box
    cx0, cx1 = x0 + 170, x1 - 115
    cy = (y0 + y1) // 2 + 12
    ry = 205
    body_top = cy - ry
    body_bot = cy + ry

    # Main cylinder body.
    body = (cx0, body_top, cx1, body_bot)
    draw.rounded_rectangle(body, radius=22, fill=(42, 55, 62, 255), outline=(13, 22, 28, 255), width=4)
    for i in range(9):
        t = i / 8
        x = cx0 + int((cx1 - cx0) * t)
        draw.line((x, body_top + 8, x, body_bot - 8), fill=(98, 116, 123, 170), width=2)
    for j in range(7):
        y = body_top + 28 + j * (body_bot - body_top - 56) / 6
        draw.arc((cx0 - 50, y - 26, cx1 + 52, y + 26), 176, 358, fill=(94, 113, 120, 130), width=2)

    # Back/front ellipses.
    draw.ellipse((cx0 - 102, cy - ry, cx0 + 102, cy + ry), fill=(64, 76, 82, 255), outline=(13, 22, 28, 255), width=4)
    draw.ellipse((cx1 - 104, cy - ry, cx1 + 104, cy + ry), fill=(36, 47, 54, 255), outline=(13, 22, 28, 255), width=5)
    for rr, col, w in [(172, full.TEAL, 6), (136, full.SPHENIX_BLUE, 4), (78, full.PHOTON_DARK, 4)]:
        draw.ellipse((cx1 - rr, cy - rr, cx1 + rr, cy + rr), outline=(*col, 210), width=w)

    # EMCal highlighted skin on the front and the top surface.
    draw.arc((cx1 - 210, cy - 210, cx1 + 210, cy + 210), 322, 42, fill=(*full.PHOTON_DARK, 255), width=13)
    draw.arc((cx1 - 177, cy - 177, cx1 + 177, cy + 177), 322, 42, fill=(*full.PHOTON, 230), width=10)
    skin = [
        (cx0 + 20, body_top + 24),
        (cx1 - 56, body_top + 42),
        (cx1 - 42, body_top + 112),
        (cx0 + 34, body_top + 92),
    ]
    draw.polygon(skin, fill=(*full.PHOTON, 74))
    draw.line(skin + [skin[0]], fill=(*full.PHOTON_DARK, 210), width=4, joint="curve")
    for i in range(12):
        t = i / 11
        xa = skin[0][0] * (1 - t) + skin[1][0] * t
        xb = skin[3][0] * (1 - t) + skin[2][0] * t
        ya = skin[0][1] * (1 - t) + skin[1][1] * t
        yb = skin[3][1] * (1 - t) + skin[2][1] * t
        draw.line((xa, ya, xb, yb), fill=(*full.PHOTON_DARK, 120), width=2)

    # Beam and tower segmentation cue.
    draw.line((cx0 - 120, cy, cx1 + 120, cy), fill=(*full.SPHENIX_BLUE, 230), width=5)
    draw.ellipse((cx1 - 18, cy - 18, cx1 + 18, cy + 18), fill=(255, 255, 255, 255), outline=(*full.SPHENIX_BLUE, 230), width=4)

    tag = (x0 + 194, body_top - 70, x0 + 542, body_top - 18)
    draw.rounded_rectangle(tag, radius=9, fill=(255, 248, 229, 248), outline=(*full.PHOTON_DARK, 180), width=2)
    center_text(draw, tag, "segmented EMCal tower skin", f(full.TIMES_BOLD, 27), full.PHOTON_DARK)
    return (cx1 + 108, cy - 4)


def draw_unrolled_sector(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> tuple[int, int]:
    x0, y0, x1, y1 = box
    top = y0 + 88
    left = x0 + 28
    right = x1 - 26
    height = y1 - y0 - 170
    strip = [
        (left, top + 34),
        (right, top - 36),
        (right - 26, top + height + 30),
        (left + 22, top + height - 8),
    ]
    draw.polygon(strip, fill=(239, 249, 252, 255))
    draw.line(strip + [strip[0]], fill=(*full.TEAL, 210), width=4, joint="curve")

    # Tower segmentation on the unrolled strip.
    for i in range(18):
        t = i / 17
        xa = strip[0][0] * (1 - t) + strip[1][0] * t
        ya = strip[0][1] * (1 - t) + strip[1][1] * t
        xb = strip[3][0] * (1 - t) + strip[2][0] * t
        yb = strip[3][1] * (1 - t) + strip[2][1] * t
        draw.line((xa, ya, xb, yb), fill=(106, 178, 195, 175), width=2)
    for j in range(7):
        t = j / 6
        xa = strip[0][0] * (1 - t) + strip[3][0] * t
        ya = strip[0][1] * (1 - t) + strip[3][1] * t
        xb = strip[1][0] * (1 - t) + strip[2][0] * t
        yb = strip[1][1] * (1 - t) + strip[2][1] * t
        draw.line((xa, ya, xb, yb), fill=(192, 217, 226, 210), width=2)

    # Peel curl.
    draw.arc((left - 26, top - 2, right + 30, top + 142), 180, 355, fill=(*full.TEAL, 230), width=7)
    draw.arc((left - 18, top + height - 78, right + 40, top + height + 72), 184, 350, fill=(*full.TEAL, 145), width=5)

    mini = (x0 + 60, y1 - 98, x1 - 60, y1 - 24)
    draw.rounded_rectangle(mini, radius=10, fill=(255, 255, 255, 245), outline=(205, 222, 232, 255), width=2)
    center_text(draw, mini, "same towers, now indexed as rows and columns", f(full.TIMES_BOLD, 28), full.BLUE)
    return (right + 22, top + height // 2)


def draw_flat_eta_phi_plane(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    gx0, gy0 = x0 + 18, y0 + 66
    cell = 39
    cols, rows = 20, 14
    gx1, gy1 = gx0 + cols * cell, gy0 + rows * cell
    draw.rounded_rectangle((gx0 - 26, gy0 - 28, gx1 + 54, gy1 + 82), radius=14, fill=(249, 252, 253, 255), outline=(214, 226, 235, 255), width=3)

    seed_c, seed_r = 10, 7
    for r in range(rows):
        for c in range(cols):
            dist2 = ((c - seed_c) / 1.45) ** 2 + ((r - seed_r) / 1.55) ** 2
            val = math.exp(-0.55 * dist2)
            x = gx0 + c * cell
            y = gy0 + r * cell
            draw.rectangle((x, y, x + cell, y + cell), fill=tower_color(val * 0.98), outline=(211, 224, 232, 255), width=2)

    win = (gx0 + 6 * cell, gy0 + 3 * cell, gx0 + 14 * cell, gy0 + 11 * cell)
    draw.rectangle(win, outline=(*full.SPHENIX_BLUE, 255), width=6)
    draw.ellipse(
        (
            gx0 + seed_c * cell + cell / 2 - 10,
            gy0 + seed_r * cell + cell / 2 - 10,
            gx0 + seed_c * cell + cell / 2 + 10,
            gy0 + seed_r * cell + cell / 2 + 10,
        ),
        fill=(0, 0, 0, 255),
    )
    local_label = (gx0 + 6 * cell - 12, gy0 + 2 * cell - 48, gx0 + 14 * cell + 12, gy0 + 2 * cell - 4)
    draw.rounded_rectangle(local_label, radius=8, fill=(255, 255, 255, 238), outline=(*full.SPHENIX_BLUE, 190), width=2)
    center_text(draw, local_label, "local cluster window", f(full.TIMES_BOLD, 27), full.SPHENIX_BLUE)

    seed_label = (gx0 + 14 * cell + 18, gy0 + 7 * cell - 22, gx1 - 18, gy0 + 7 * cell + 28)
    draw.rounded_rectangle(seed_label, radius=8, fill=(255, 255, 255, 235), outline=(203, 216, 227, 255), width=2)
    center_text(draw, seed_label, "seed", f(full.TIMES_ITALIC, 25), full.MUTED)

    arrow(draw, (gx0 - 8, gy1 + 38), (gx1 + 10, gy1 + 38), full.SPHENIX_BLUE, width=5)
    arrow(draw, (gx0 - 38, gy1 + 8), (gx0 - 38, gy0 + 4), full.TEAL, width=5)
    draw.text((gx1 + 26, gy1 + 17), "η", font=f(full.TIMES_ITALIC, 43), fill=full.SPHENIX_BLUE)
    draw.text((gx0 - 80, gy0 - 45), "φ", font=f(full.TIMES_ITALIC, 43), fill=full.TEAL)


def draw_tower_grid_setup_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    """Single large geometry panel: barrel readout -> global grid -> local seed window."""
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(
        box,
        radius=12,
        fill=(249, 252, 253, 255),
        outline=(218, 229, 237, 255),
        width=2,
    )

    center_text(
        draw,
        (x0 + 24, y0 + 8, x1 - 24, y0 + 58),
        "CEMC barrel readout → flat η–φ tower image",
        f(full.TIMES_BOLD, 39),
        full.INK,
    )

    # Barrel annulus on the left.
    cx, cy = x0 + 152, y0 + 185
    outer = (cx - 86, cy - 86, cx + 86, cy + 86)
    inner = (cx - 54, cy - 54, cx + 54, cy + 54)
    draw.ellipse(outer, fill=(246, 250, 252, 255), outline=(49, 71, 84, 255), width=5)
    draw.ellipse(inner, fill=(255, 255, 255, 255), outline=(169, 188, 198, 255), width=4)
    for deg in range(0, 360, 10):
        draw.arc(outer, deg + 1, deg + 6, fill=(*full.TEAL, 205), width=5)
    draw.arc(outer, 314, 46, fill=(*full.PHOTON_DARK, 255), width=15)
    draw.arc((cx - 78, cy - 78, cx + 78, cy + 78), 314, 46, fill=(*full.PHOTON, 230), width=10)
    center_text(draw, (cx - 98, cy + 96, cx + 98, cy + 140), "CEMC barrel", f(full.TIMES_BOLD, 28), full.TEAL)

    arrow(draw, (x0 + 270, y0 + 185), (x0 + 360, y0 + 185), full.PHOTON_DARK, width=6)

    # Global tower image.
    gx0, gy0 = x0 + 392, y0 + 94
    gx1, gy1 = x1 - 64, y0 + 292
    draw.rounded_rectangle((gx0, gy0, gx1, gy1), radius=8, fill=(255, 255, 255, 255), outline=(*full.SPHENIX_BLUE, 230), width=4)
    n_eta, n_phi = 24, 10
    for i in range(1, n_eta):
        xx = gx0 + (gx1 - gx0) * i / n_eta
        draw.line((xx, gy0, xx, gy1), fill=(203, 220, 230, 255), width=2)
    for j in range(1, n_phi):
        yy = gy0 + (gy1 - gy0) * j / n_phi
        draw.line((gx0, yy, gx1, yy), fill=(203, 220, 230, 255), width=2)

    wx0 = gx0 + int((gx1 - gx0) * 0.60)
    wy0 = gy0 + int((gy1 - gy0) * 0.34)
    wx1 = gx0 + int((gx1 - gx0) * 0.78)
    wy1 = gy0 + int((gy1 - gy0) * 0.68)
    draw.rectangle((wx0, wy0, wx1, wy1), fill=(*full.PHOTON, 64), outline=(*full.SPHENIX_BLUE, 255), width=5)
    draw.ellipse(((wx0 + wx1) / 2 - 8, (wy0 + wy1) / 2 - 8, (wx0 + wx1) / 2 + 8, (wy0 + wy1) / 2 + 8), fill=(0, 0, 0, 255))

    center_text(draw, (gx0, gy0 - 44, gx1, gy0 - 9), "global η–φ tower plane", f(full.TIMES_BOLD, 33), full.BLUE)
    small_arrow(draw, (gx0 - 18, gy1 + 28), (gx1 + 10, gy1 + 28), full.SPHENIX_BLUE, width=4, head=13)
    small_arrow(draw, (gx0 - 32, gy1 + 2), (gx0 - 32, gy0 - 2), full.TEAL, width=4, head=13)
    global_phi_label = (gx0 - 74, gy0 - 8, gx0 - 35, gy0 + 43)
    global_eta_label = (gx1 + 1, gy1 - 25, gx1 + 43, gy1 + 15)
    draw.rounded_rectangle(global_phi_label, radius=7, fill=(249, 252, 253, 232), outline=(207, 222, 232, 150), width=1)
    draw.rounded_rectangle(global_eta_label, radius=7, fill=(249, 252, 253, 232), outline=(207, 222, 232, 150), width=1)
    center_text(draw, global_phi_label, "φ", f(full.TIMES_ITALIC, 39), full.TEAL)
    center_text(draw, global_eta_label, "η", f(full.TIMES_ITALIC, 39), full.SPHENIX_BLUE)

    # Tower specification strip.
    spec_y0 = y0 + 338
    spec = (x0 + 48, spec_y0, x1 - 48, spec_y0 + 78)
    draw.rounded_rectangle(spec, radius=10, fill=(255, 252, 244, 255), outline=(*full.PHOTON_DARK, 185), width=2)
    spec_cols = [
        (spec[0] + 12, spec[1] + 8, spec[0] + (spec[2] - spec[0]) // 3 - 6, spec[3] - 8),
        (spec[0] + (spec[2] - spec[0]) // 3 + 6, spec[1] + 8, spec[0] + 2 * (spec[2] - spec[0]) // 3 - 6, spec[3] - 8),
        (spec[0] + 2 * (spec[2] - spec[0]) // 3 + 6, spec[1] + 8, spec[2] - 12, spec[3] - 8),
    ]
    for xx in [spec_cols[1][0] - 12, spec_cols[2][0] - 12]:
        draw.line((xx, spec[1] + 12, xx, spec[3] - 12), fill=(230, 215, 180, 255), width=2)
    center_text(draw, spec_cols[0], "96 η × 256 φ bins", f(full.TIMES_BOLD, 31), full.BLUE)
    center_text(draw, spec_cols[1], "24,576 towers", f(full.TIMES_BOLD, 31), full.PHOTON_DARK)
    center_text(draw, spec_cols[2], "Δη × Δφ = 0.025 × 0.025", f(full.TIMES_BOLD, 29), full.SPHENIX_BLUE)

    # Zoom into the local seed window. The highlighted region in the global map
    # and the local-window label provide the connection without crossing the
    # tower-spec strip.
    lx0, ly0 = x0 + 78, y0 + 470
    cell = 23
    cols, rows = 13, 8
    lx1, ly1 = lx0 + cols * cell, ly0 + rows * cell
    zoom_frame = (lx0 - 26, ly0 - 34, lx1 + 42, ly1 + 52)
    draw.rounded_rectangle(zoom_frame, radius=12, fill=(255, 255, 255, 255), outline=(210, 224, 234, 255), width=2)
    seed_c, seed_r = cols // 2, rows // 2
    for r in range(rows):
        for c in range(cols):
            dist2 = ((c - seed_c) / 1.35) ** 2 + ((r - seed_r) / 1.45) ** 2
            val = math.exp(-0.58 * dist2)
            xx = lx0 + c * cell
            yy = ly0 + r * cell
            draw.rectangle((xx, yy, xx + cell, yy + cell), fill=tower_color(val * 0.98), outline=(211, 224, 232, 255), width=2)
    win = (lx0 + 4 * cell, ly0 + 2 * cell, lx0 + 9 * cell, ly0 + 7 * cell)
    draw.rectangle(win, outline=(*full.SPHENIX_BLUE, 255), width=5)
    seed_x = lx0 + seed_c * cell + cell / 2
    seed_y = ly0 + seed_r * cell + cell / 2
    draw.ellipse(
        (
            seed_x - 9,
            seed_y - 9,
            seed_x + 9,
            seed_y + 9,
        ),
        fill=(0, 0, 0, 255),
    )
    center_text(draw, (zoom_frame[0] + 8, zoom_frame[1] + 4, zoom_frame[2] - 8, ly0 - 5), "local seed-centered window", f(full.TIMES_BOLD, 29), full.BLUE)
    # Local coordinate axes and seed callout.
    phi_x = lx0 - 22
    eta_y = ly1 + 26
    small_arrow(draw, (phi_x, ly1 + 18), (phi_x, ly0 + 4), full.TEAL, width=4, head=12)
    small_arrow(draw, (lx0 - 2, eta_y), (lx1 + 18, eta_y), full.SPHENIX_BLUE, width=4, head=12)
    local_phi_label = (phi_x - 35, ly0 + 62, phi_x - 3, ly0 + 106)
    local_eta_label = (lx1 - 4, eta_y - 29, lx1 + 38, eta_y + 10)
    draw.rounded_rectangle(local_phi_label, radius=6, fill=(255, 255, 255, 232), outline=(204, 219, 229, 150), width=1)
    draw.rounded_rectangle(local_eta_label, radius=6, fill=(255, 255, 255, 232), outline=(204, 219, 229, 150), width=1)
    center_text(draw, local_phi_label, "φ", f(full.TIMES_ITALIC, 33), full.TEAL)
    center_text(draw, local_eta_label, "η", f(full.TIMES_ITALIC, 33), full.SPHENIX_BLUE)

    seed_box = (lx1 + 10, int(seed_y - 25), lx1 + 82, int(seed_y + 17))
    draw.rounded_rectangle(seed_box, radius=7, fill=(255, 255, 255, 238), outline=(195, 208, 218, 255), width=2)
    center_text(draw, seed_box, "seed", f(full.TIMES_ITALIC, 24), full.MUTED)
    draw.line((seed_box[0] - 4, (seed_box[1] + seed_box[3]) // 2, seed_x + 12, seed_y - 2), fill=(95, 108, 120, 210), width=3)

    explain = (x0 + 475, y0 + 450, x1 - 56, y0 + 696)
    draw.rounded_rectangle(explain, radius=10, fill=(239, 247, 251, 255), outline=(207, 222, 232, 255), width=2)
    draw_wrapped_centered(
        draw,
        "Shower-shape variables summarize the local tower-energy pattern around the cluster seed.",
        explain,
        f(full.TIMES_BOLD, 35),
        full.BLUE,
        line_gap=10,
    )

def draw_slide() -> Image.Image:
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*full.SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*full.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*full.PHOTON, 255))
    title = "Backup: EMCal tower-map geometry"
    subtitle = "CEMC tower readout becomes the η–φ image used for local shower-shape inputs."
    draw.text((132, 76), title, font=f(full.TIMES_BOLD, 86), fill=full.INK)
    draw.text(
        (136, 190),
        subtitle,
        font=f(full.TIMES_ITALIC, 56),
        fill=full.MUTED,
    )
    draw.line((132, 286, W - 132, 286), fill=(221, 226, 232, 255), width=3)
    full.add_top_right_sphenix_logo_like_slide2(img)

    content_top = 324
    content_bottom = 1236
    gap = 40
    left_w = 760
    left = (132, content_top, 132 + left_w, content_bottom)
    right = (left[2] + gap, content_top, W - 132, content_bottom)

    draw_real_card(
        img,
        draw,
        left,
        full.PHOTON_DARK,
        "Detector geometry",
        "EMCal barrel layer around the tracking volume",
    )
    det_area = (left[0] + 42, left[1] + 148, left[2] - 42, left[1] + 672)
    paste_detector_with_emcal_highlight(img, draw, load_detector_anchor(), det_area)
    tag = (left[0] + 92, left[1] + 700, left[2] - 92, left[1] + 764)
    draw.rounded_rectangle(tag, radius=10, fill=(255, 248, 229, 238), outline=(*full.PHOTON_DARK, 210), width=2)
    center_text(draw, tag, "CEMC / EMCal", f(full.TIMES_BOLD, 34), full.PHOTON_DARK)

    spec = (left[0] + 56, left[1] + 792, left[2] - 56, left[3] - 52)
    draw.rounded_rectangle(spec, radius=9, fill=(247, 251, 253, 244), outline=(207, 223, 232, 255), width=2)
    draw_wrapped_centered(
        draw,
        "tungsten / scintillating-fiber\nsampling calorimeter",
        spec,
        f(full.TIMES_BOLD, 31),
        full.BLUE,
        line_gap=5,
    )

    draw_real_card(
        img,
        draw,
        right,
        full.SPHENIX_BLUE,
        "Tower-grid readout",
        "global η–φ tower image, then local seed window",
    )
    grid_panel = (right[0] + 42, right[1] + 148, right[2] - 42, right[3] - 52)
    draw_tower_grid_setup_panel(draw, grid_panel)

    full.draw_hp2026_identity_footer(img)
    return img.convert("RGB")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    img = draw_slide()
    img.save(OUT_PNG, quality=95)
    OUT_SCRIPT.write_text(
        """# Backup - EMCal tower-map geometry

This backup slide gives the detector geometry behind the shower-shape maps.

The sPHENIX CEMC is the electromagnetic calorimeter layer around the tracking volume. It is a tungsten and scintillating-fiber sampling calorimeter, segmented into a tower readout.

The important point for the shower-shape slides is not the construction hierarchy. The important point is that the cylindrical barrel becomes a 96 by 256 eta-phi tower image, with 24,576 EMCal towers total and tower granularity of 0.025 by 0.025.

For a photon candidate, the shower-shape variables use a local window around the cluster seed, so the physics object is a small eta-phi tower-energy image.
""",
        encoding="utf-8",
    )
    OUT_HEADER.write_text(
        json.dumps(
            {
                "hp2026_main_header": {
                    "deck": "hp2026_main_talk",
                    "title_font_size": 86,
                    "subtitle_font_size": 56,
                    "title_xy": [132, 76],
                    "subtitle_xy": [136, 190],
                    "divider_y": 286,
                }
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    meta = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "google_slides_mutation": False,
        "output_png": str(OUT_PNG.relative_to(ROOT)),
        "speaker_script": str(OUT_SCRIPT.relative_to(ROOT)),
        "source_script": str(Path(__file__).relative_to(ROOT)),
        "canvas": [W, H],
        "source_inputs": [
            str(DETECTOR_RENDER.relative_to(ROOT)) if DETECTOR_RENDER.exists() else str(DETECTOR_ASSET.relative_to(ROOT)),
            "ppg12codeGit/wiki/physics/detector/emcal-performance.md",
            "ppg12codeGit/wiki/physics/detector/sphenix-overview.md",
        ],
        "design_intent": "Backup geometry slide explaining how the CEMC barrel maps to a global eta-phi tower plane and then to the local tower image used for shower-shape inputs.",
    }
    OUT_META.write_text(json.dumps(meta, indent=2) + "\n")
    print(OUT_PNG)


if __name__ == "__main__":
    main()

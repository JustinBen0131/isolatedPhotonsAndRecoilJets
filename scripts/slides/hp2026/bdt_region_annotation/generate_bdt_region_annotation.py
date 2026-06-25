#!/usr/bin/env python3
"""Annotate the PPG12 photon-ID BDT score figure for HP2026 slides."""

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


REPO = Path(__file__).resolve().parents[4]
BASE = (
    REPO
    / "outputs/manual-20260601-hp2026-fulltalk/presentations/hp2026-fulltalk"
    / "assets/paper_figures/fig2_bdt_score.png"
)
OUT = REPO / "outputs/hp2026_bdt_region_annotation/bdt_score_regions_signal_background_neither.png"


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
    ]
    for candidate in candidates:
        path = Path(candidate)
        if path.exists():
            return ImageFont.truetype(str(path), size=size)
    return ImageFont.load_default()


def score_to_x(score: float, *, scale: int) -> int:
    # Pixel coordinates of the plot frame in the extracted figure.
    x0, x1 = 186 * scale, 959 * scale
    return round(x0 + score * (x1 - x0))


def draw_double_arrow(draw: ImageDraw.ImageDraw, x0: int, x1: int, y: int, color: tuple[int, int, int], width: int) -> None:
    draw.line((x0, y, x1, y), fill=color, width=width)
    head = 16
    half = 10
    draw.polygon([(x0, y), (x0 + head, y - half), (x0 + head, y + half)], fill=color)
    draw.polygon([(x1, y), (x1 - head, y - half), (x1 - head, y + half)], fill=color)


def draw_label(
    draw: ImageDraw.ImageDraw,
    center: tuple[int, int],
    text: str,
    text_font: ImageFont.FreeTypeFont,
    color: tuple[int, int, int],
    outline: tuple[int, int, int],
    pad_x: int,
    pad_y: int,
) -> None:
    cx, cy = center
    bbox = draw.textbbox((0, 0), text, font=text_font)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    box = (cx - tw // 2 - pad_x, cy - th // 2 - pad_y, cx + tw // 2 + pad_x, cy + th // 2 + pad_y)
    draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255), outline=outline, width=3)
    draw.text((cx - tw // 2, cy - th // 2 - 1), text, font=text_font, fill=color)


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    scale = 2
    img = Image.open(BASE).convert("RGB")
    img = img.resize((img.width * scale, img.height * scale), Image.Resampling.LANCZOS)
    draw = ImageDraw.Draw(img)

    blue = (42, 143, 229)
    red = (231, 86, 73)
    gray = (110, 110, 110)
    light_gray = (210, 210, 210)

    # Approximate 14 < ET < 18 GeV annotation at ET ~= 16 GeV:
    # non-tight lower edge = 0.7333 - 0.01333 ET ~= 0.52
    # non-tight upper edge = 0.6844 + 0.00156 ET ~= 0.71
    # tight threshold = 0.8156 - 0.00156 ET ~= 0.79
    nt_lower = 0.52
    nt_upper = 0.71
    tight = 0.79
    y_arrow = 442 * scale
    y_label = 418 * scale

    draw_double_arrow(draw, score_to_x(nt_lower, scale=scale), score_to_x(nt_upper, scale=scale), y_arrow, blue, 5)
    draw_double_arrow(draw, score_to_x(tight + 0.015, scale=scale), score_to_x(0.98, scale=scale), y_arrow, red, 5)

    label_font = font(28, bold=True)
    gap_font = font(30, bold=True)
    draw_label(
        draw,
        ((score_to_x(nt_lower, scale=scale) + score_to_x(nt_upper, scale=scale)) // 2, y_label),
        "Background-enhanced",
        label_font,
        blue,
        blue,
        14,
        8,
    )
    draw_label(
        draw,
        ((score_to_x(tight + 0.015, scale=scale) + score_to_x(0.98, scale=scale)) // 2, y_label),
        "Signal-enhanced",
        label_font,
        red,
        red,
        14,
        8,
    )

    x_nt_lower = score_to_x(nt_lower, scale=scale)
    x_gap0 = score_to_x(nt_upper, scale=scale)
    x_gap1 = score_to_x(tight, scale=scale)
    y0, y1 = 382 * scale, 486 * scale
    for x in (x_nt_lower, x_gap0, x_gap1):
        for yy in range(y0, y1, 18):
            draw.line((x, yy, x, min(yy + 9, y1)), fill=gray, width=3)
    draw.line((x_gap0, y_arrow, x_gap1, y_arrow), fill=light_gray, width=5)
    draw_label(
        draw,
        ((x_gap0 + x_gap1) // 2, y0 - 28),
        "neither",
        gap_font,
        gray,
        gray,
        12,
        6,
    )

    img.save(OUT, dpi=(300, 300))
    print(OUT)


if __name__ == "__main__":
    main()

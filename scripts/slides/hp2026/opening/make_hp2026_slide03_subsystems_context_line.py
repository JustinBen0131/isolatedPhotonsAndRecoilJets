#!/usr/bin/env python3
"""Local Slide 3 refinement: add a concise sPHENIX context line.

This uses the current accepted Slide 3 PNG as a raster baseline. It preserves
the title, footer, right detector rendering, and subsystem-card internals while
moving the left card stack down to make room for a single audience-facing
context statement above it.

No Google Slides mutation is performed.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw


ROOT = next(
    p
    for p in Path(__file__).resolve().parents
    if (p / "README.md").exists() and (p / "scripts").exists() and (p / "src").exists()
)
SCRIPT_DIR = ROOT / "scripts/slides/hp2026/opening"
sys.path.insert(0, str(SCRIPT_DIR))

import make_hp2026_opening_motivation_slide as opening  # noqa: E402


BASELINE = (
    ROOT
    / "outputs/manual-20260611-slides3-7-header-normalized/"
    / "slide03_sphenix_subsystems_header_normalized.png"
)
OUT_DIR = ROOT / "outputs/manual-20260611-slide03-subsystems-context-line"
FINAL_PNG = OUT_DIR / "slide03_sphenix_subsystems_context_line.png"
MANIFEST = OUT_DIR / "slide03_sphenix_subsystems_context_line_manifest.json"


def text_width(draw: ImageDraw.ImageDraw, text: str, font) -> int:
    box = draw.textbbox((0, 0), text, font=font)
    return box[2] - box[0]


def draw_rich_context(draw: ImageDraw.ImageDraw) -> None:
    # Deck-current baseline: a two-line intro above the subsystem cards, not a
    # boxed callout or RHS caption.
    x = 132
    line1_y = 258
    line2_y = 314

    line1_bold = opening.font(opening.TIMES_BOLD, 45)
    line1_body = opening.font(opening.TIMES, 45)
    line2_body = opening.font(opening.TIMES, 45)
    line2_bold = opening.font(opening.TIMES_BOLD, 45)

    line1_parts = [
        ("sPHENIX:", line1_bold, opening.INK),
        (" full-azimuth RHIC detector at BNL", line1_body, opening.INK),
    ]
    cursor = x
    for text, font, fill in line1_parts:
        draw.text((cursor, line1_y), text, font=font, fill=fill)
        cursor += text_width(draw, text, font)

    line2_parts = [
        ("central coverage ", line2_body, opening.MUTED),
        ("|\u03b7| < 1.1", line2_bold, opening.INK),
    ]
    cursor = x
    for text, font, fill in line2_parts:
        draw.text((cursor, line2_y), text, font=font, fill=fill)
        cursor += text_width(draw, text, font)


def render() -> Path:
    if not BASELINE.exists():
        raise FileNotFoundError(f"missing current Slide 3 baseline: {BASELINE}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    img = Image.open(BASELINE).convert("RGBA")
    draw = ImageDraw.Draw(img, "RGBA")
    bg = (*opening.SOFT_BG, 255)

    # Preserve all card pixels exactly by moving one crop.
    crop_box = (52, 340, 1240, 1176)
    shift_x = 36
    shift_y = 88
    cards = img.crop(crop_box)

    # Clear old card region plus the new statement band. Keep title/footer/right
    # detector untouched.
    draw.rectangle((52, 260, 1328, 1268), fill=bg)
    draw_rich_context(draw)
    img.alpha_composite(cards, (crop_box[0] + shift_x, crop_box[1] + shift_y))

    img.convert("RGB").save(FINAL_PNG, "PNG")
    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "google_slides_mutation": False,
        "baseline_png": str(BASELINE.relative_to(ROOT)),
        "output_png": str(FINAL_PNG.relative_to(ROOT)),
        "change_scope": "Deck-current baseline: move left subsystem-card stack down/right and add the two-line sPHENIX identity/coverage text above it; title/footer/right detector preserved.",
        "left_card_crop": crop_box,
        "left_card_shift_x_px": shift_x,
        "left_card_shift_y_px": shift_y,
        "context_text": "sPHENIX: full-azimuth RHIC detector at BNL / central coverage |eta| < 1.1",
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return FINAL_PNG


def main() -> None:
    print(render())


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Refine the HP2026 Slide 3 subsystem header without changing the body.

The current deck slide uses a locally reordered subsystem PNG. This helper
uses that PNG as the baseline and redraws only the header band so the slide
matches the HP2026 main-talk title/subtitle contract.
"""

from __future__ import annotations

import json
import shutil
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


BASELINE_PNG = (
    ROOT
    / "outputs/manual-20260608-slide3-subsystems-reorder/"
    / "hp2026_slide03_sphenix_subsystems_order_calorimetry_forward_tracking.png"
)
OUTDIR = ROOT / "outputs/manual-20260610-slide03_subsystems_header_refined"
FINAL_PNG = OUTDIR / "slide03_sphenix_subsystems_header_refined.png"
SCRIPT_SRC = (
    ROOT
    / "outputs/manual-20260601-hp2026-opening-slide/presentations/hp2026-opening-slide/output/"
    / "hp2026_slide02_sphenix_subsystems_recreated_script.md"
)
SCRIPT_OUT = OUTDIR / "slide03_sphenix_subsystems_header_refined_script.md"
MANIFEST = OUTDIR / "slide03_sphenix_subsystems_header_refined_manifest.json"

TITLE = "sPHENIX subsystems"
SUBTITLE = "A full-azimuth RHIC experiment at BNL with central coverage |η| < 1.1."

HP2026_MAIN_HEADER = {
    "deck": "hp2026_main_talk",
    "title_font_size": 86,
    "subtitle_font_size": 56,
    "title_xy": [132, 76],
    "subtitle_xy": [136, 190],
    "divider_y": 286,
}


def write_header_spec(path: Path) -> None:
    path.with_suffix(".header.json").write_text(
        json.dumps({"hp2026_main_header": HP2026_MAIN_HEADER}, indent=2) + "\n",
        encoding="utf-8",
    )


def render() -> Path:
    if not BASELINE_PNG.exists():
        raise FileNotFoundError(f"missing current reordered baseline: {BASELINE_PNG}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    img = Image.open(BASELINE_PNG).convert("RGBA")
    draw = ImageDraw.Draw(img, "RGBA")

    # Preserve the top blue/yellow rule and the body below the header divider.
    draw.rectangle((0, 30, opening.W, HP2026_MAIN_HEADER["divider_y"] - 1), fill=(*opening.SOFT_BG, 255))
    draw.text(
        tuple(HP2026_MAIN_HEADER["title_xy"]),
        TITLE,
        font=opening.font(opening.TIMES_BOLD, HP2026_MAIN_HEADER["title_font_size"]),
        fill=opening.INK,
    )
    draw.text(
        tuple(HP2026_MAIN_HEADER["subtitle_xy"]),
        SUBTITLE,
        font=opening.font(opening.TIMES_ITALIC, HP2026_MAIN_HEADER["subtitle_font_size"]),
        fill=opening.MUTED,
    )
    draw.line(
        (132, HP2026_MAIN_HEADER["divider_y"], opening.W - 132, HP2026_MAIN_HEADER["divider_y"]),
        fill=(221, 226, 232, 255),
        width=3,
    )

    logo = opening.load_sphenix_logo()
    if logo is not None:
        opening.paste_fit(img, logo, (2188, 58, 2432, 164), anchor="right")

    img.convert("RGB").save(FINAL_PNG, "PNG")
    write_header_spec(FINAL_PNG)

    if SCRIPT_SRC.exists():
        shutil.copy2(SCRIPT_SRC, SCRIPT_OUT)

    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "google_slides_mutation": False,
        "baseline_png": str(BASELINE_PNG.relative_to(ROOT)),
        "output_png": str(FINAL_PNG.relative_to(ROOT)),
        "speaker_script": str(SCRIPT_OUT.relative_to(ROOT)) if SCRIPT_OUT.exists() else None,
        "change_scope": "Header band only: title, subtitle, divider, and sPHENIX logo redrawn; slide body/footer preserved from current reordered PNG.",
        "title": TITLE,
        "subtitle": SUBTITLE,
        "hp2026_main_header": HP2026_MAIN_HEADER,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return FINAL_PNG


def main() -> None:
    print(render())


if __name__ == "__main__":
    main()

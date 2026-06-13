#!/usr/bin/env python3
"""Normalize HP2026 Slides 3-7 headers to the current Slide 8-10 contract.

This is a local PNG-only pass.  It preserves each accepted slide body and footer
as raster input and redraws only the title/subtitle/divider header band:

- remove the subtitle
- redraw the title at the Slide 8/9/10 location and size
- redraw the title divider line at the Slide 8/9/10 y-position

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

import make_hp2026_opening_motivation_slide as hp  # noqa: E402


OUT_DIR = ROOT / "outputs/manual-20260611-slides3-7-header-normalized"

HEADER_SPEC = {
    "title_font_size": 86,
    "title_xy": [132, 76],
    "subtitle": None,
    "divider_y": 232,
    "divider_x0": 132,
    "divider_x1": 2428,
}

SLIDES = [
    {
        "slide": 3,
        "title": "sPHENIX subsystems",
        "source": ROOT
        / "outputs/manual-20260610-slide03_subsystems_header_refined/"
        / "slide03_sphenix_subsystems_header_refined.png",
        "output": "slide03_sphenix_subsystems_header_normalized.png",
    },
    {
        "slide": 4,
        "title": "Run 24 p+p dataset used for this measurement",
        "source": ROOT
        / "outputs/manual-20260607-run24pp-lumi-asset/"
        / "hp2026_slide04_analysis_dataset_context_v13_lint_expanded.png",
        "output": "slide04_run24pp_dataset_header_normalized.png",
    },
    {
        "slide": 5,
        "title": "Prompt photons as color-neutral hard probes",
        "source": ROOT
        / "outputs/manual-20260609-slide5_7_prompt_photons_expand_sequence/"
        / "slide05_state01_production_expanded.png",
        "output": "slide05_prompt_photons_production_header_normalized.png",
    },
    {
        "slide": 6,
        "title": "Prompt photons as color-neutral hard probes",
        "source": ROOT
        / "outputs/manual-20260609-slide5_7_prompt_photons_expand_sequence/"
        / "slide06_state02_isolation_expanded.png",
        "output": "slide06_prompt_photons_isolation_header_normalized.png",
    },
    {
        "slide": 7,
        "title": "Prompt photons as color-neutral hard probes",
        "source": ROOT
        / "outputs/manual-20260609-slide5_7_prompt_photons_expand_sequence/"
        / "slide07_state03_current_full_context.png",
        "output": "slide07_prompt_photons_context_header_normalized.png",
    },
]


def normalize_header(entry: dict[str, object]) -> Path:
    source = Path(entry["source"])
    if not source.exists():
        raise FileNotFoundError(f"Missing source PNG for slide {entry['slide']}: {source}")

    img = Image.open(source).convert("RGBA")
    draw = ImageDraw.Draw(img, "RGBA")
    bg = (252, 253, 254, 255)

    # Clear the text header band without touching the top blue/yellow rule or
    # the right sPHENIX logo. The source body/footer stay unchanged.
    draw.rectangle((88, 48, 2156, 258), fill=bg)
    # Subtitle lines can run nearly to the right margin; this band is below the
    # logo, so clearing it fully does not disturb the logo artwork.
    draw.rectangle((88, 170, 2508, 258), fill=bg)

    # Remove the older lower divider used by prior header variants.
    draw.rectangle((HEADER_SPEC["divider_x0"], 270, HEADER_SPEC["divider_x1"], 300), fill=bg)

    draw.text(
        tuple(HEADER_SPEC["title_xy"]),
        str(entry["title"]),
        font=hp.font(hp.TIMES_BOLD, HEADER_SPEC["title_font_size"]),
        fill=hp.INK,
    )
    draw.line(
        (
            HEADER_SPEC["divider_x0"],
            HEADER_SPEC["divider_y"],
            HEADER_SPEC["divider_x1"],
            HEADER_SPEC["divider_y"],
        ),
        fill=(221, 226, 232, 255),
        width=3,
    )

    out = OUT_DIR / str(entry["output"])
    img.convert("RGB").save(out, "PNG")
    return out


def make_contact_sheet(paths: list[Path]) -> Path:
    thumbs = []
    for entry, path in zip(SLIDES, paths):
        im = Image.open(path).convert("RGB").resize((640, 360))
        d = ImageDraw.Draw(im)
        d.rounded_rectangle((8, 8, 58, 42), radius=6, fill=(255, 255, 255), outline=(210, 218, 228), width=2)
        d.text((24, 11), str(entry["slide"]), font=hp.font(hp.TIMES_BOLD, 24), fill=hp.INK)
        thumbs.append(im)

    sheet = Image.new("RGB", (1280, 1080), (248, 250, 252))
    for i, thumb in enumerate(thumbs):
        sheet.paste(thumb, ((i % 2) * 640, (i // 2) * 360))
    out = OUT_DIR / "slides03_07_header_normalized_contact_sheet.png"
    sheet.save(out, "PNG")
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    paths = [normalize_header(entry) for entry in SLIDES]
    contact = make_contact_sheet(paths)
    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "google_slides_mutation": False,
        "scope": "Header-only raster regeneration for current accepted Slides 3-7.",
        "header_spec": HEADER_SPEC,
        "outputs": [str(p.relative_to(ROOT)) for p in paths],
        "contact_sheet": str(contact.relative_to(ROOT)),
    }
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    for p in paths:
        print(p)
    print(contact)


if __name__ == "__main__":
    main()

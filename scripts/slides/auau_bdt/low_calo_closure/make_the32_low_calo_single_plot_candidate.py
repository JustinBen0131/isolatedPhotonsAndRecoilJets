#!/usr/bin/env python3
"""Build a one-plot THE-32 low-calo explanation candidate."""

from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


SLIDE_W = 2560
SLIDE_H = 1440
OUTDIR = Path("dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603")
SOURCE_PNG = OUTDIR / "the32_low_calo_pathology_slide_v1.png"
SOURCE_JSON = OUTDIR / "the32_low_calo_pathology_slide_v1.json"
OUT_PNG = OUTDIR / "the32_low_calo_single_plot_explanation_v1.png"
OUT_MD = OUTDIR / "the32_low_calo_single_plot_explanation_v1.md"
OUT_MANIFEST = OUTDIR / "the32_low_calo_single_plot_explanation_v1.json"

INK = (24, 34, 48)
MUTED = (71, 84, 103)
RED = (180, 35, 24)
BLUE = (31, 111, 179)
PALE_YELLOW = (255, 247, 214)
YELLOW_BORDER = (227, 179, 65)
BORDER = (195, 205, 219)


def font(size: int, *, bold: bool = False, italic: bool = False) -> ImageFont.FreeTypeFont:
    candidates = []
    if bold and italic:
        candidates = [
            "/System/Library/Fonts/Supplemental/Times New Roman Bold Italic.ttf",
            "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf",
        ]
    elif bold:
        candidates = [
            "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf",
            "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
        ]
    elif italic:
        candidates = [
            "/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf",
            "/System/Library/Fonts/Supplemental/Arial Italic.ttf",
        ]
    else:
        candidates = [
            "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
            "/System/Library/Fonts/Supplemental/Arial.ttf",
        ]
    for path in candidates:
        if Path(path).exists():
            return ImageFont.truetype(path, size=size)
    return ImageFont.load_default(size=size)


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    *,
    max_width: int,
    fill: tuple[int, int, int],
    text_font: ImageFont.FreeTypeFont,
    line_spacing: int = 10,
) -> int:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        trial = word if not current else f"{current} {word}"
        if draw.textbbox((0, 0), trial, font=text_font)[2] <= max_width:
            current = trial
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    x, y = xy
    for line in lines:
        draw.text((x, y), line, font=text_font, fill=fill)
        y += text_font.size + line_spacing
    return y


def rounded_box(draw: ImageDraw.ImageDraw, box, fill, outline, width=3, radius=20) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def main() -> int:
    source = Image.open(SOURCE_PNG).convert("RGB")
    # Crop only the evidence plot from the validated slide: axes, red rejected tail,
    # black low-calo envelope, legend, and density colorbar. Side panels are omitted.
    plot_crop = source.crop((15, 300, 1710, 1060))

    canvas = Image.new("RGB", (SLIDE_W, SLIDE_H), "white")
    draw = ImageDraw.Draw(canvas)

    title_font = font(62, bold=True)
    sub_font = font(34)
    small_font = font(30)
    note_font = font(38)
    note_bold = font(42, bold=True)

    draw.text((105, 70), "The issue: some events have too little calorimeter energy for their centrality", font=title_font, fill=INK)
    draw.text(
        (110, 160),
        "This one plot defines the event-quality cut: reject the red tail below the black centrality-conditioned envelope.",
        font=sub_font,
        fill=MUTED,
    )

    target_w = 2100
    target_h = int(plot_crop.height * target_w / plot_crop.width)
    plot_big = plot_crop.resize((target_w, target_h), Image.Resampling.LANCZOS)
    px = (SLIDE_W - target_w) // 2
    py = 245
    draw.rounded_rectangle((px - 18, py - 18, px + target_w + 18, py + target_h + 18), radius=16, fill=(247, 249, 252), outline=BORDER, width=2)
    canvas.paste(plot_big, (px, py))

    with SOURCE_JSON.open() as f:
        pathology = json.load(f)
    event_rows = pathology["event_counts"]
    total_events = sum(int(r["event_total"]) for r in event_rows)
    rejected_events = sum(int(r["event_rejected"]) for r in event_rows)
    rejected_frac = rejected_events / total_events if total_events else 0.0

    rounded_box(draw, (120, 1225, 2440, 1350), fill=PALE_YELLOW, outline=YELLOW_BORDER, width=3, radius=18)
    draw.text((165, 1252), "Clean claim:", font=note_bold, fill=INK)
    claim = (
        f"Using only centrality and total calo energy, the diagnostic finds "
        f"{rejected_events:,}/{total_events:,} events ({rejected_frac:.1%}) in the low-calo tail. "
        "That is the pathology the upstream cut is designed to remove before BDT training."
    )
    draw_wrapped(draw, (420, 1255), claim, max_width=1950, fill=INK, text_font=small_font, line_spacing=8)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    canvas.save(OUT_PNG)

    OUT_MD.write_text(
        "\n".join(
            [
                "# THE-32 single-plot low-calo explanation",
                "",
                "This candidate removes the side-panel bar chart and source table from the original pathology slide.",
                "The spoken point is: the red event tail sits below the normal total-calo-energy envelope at the same centrality, so the cut removes a centrality-mismatched event-quality tail before BDT training.",
                "",
            ]
        )
    )
    OUT_MANIFEST.write_text(
        json.dumps(
            {
                "schema": "THE32_LOW_CALO_SINGLE_PLOT_EXPLANATION_V1",
                "source_png": str(SOURCE_PNG),
                "source_json": str(SOURCE_JSON),
                "output_png": str(OUT_PNG),
                "output_script": str(OUT_MD),
                "method": "Crop the validated main pathology panel only; no data recomputation or relabeling.",
                "crop_box_px": [15, 300, 1710, 1060],
                "event_total": total_events,
                "event_rejected": rejected_events,
                "event_rejected_fraction": rejected_frac,
                "cut_inputs": "centrality and log10(CEMC + IHCal + OHCal + 1) only",
                "excluded_from_candidate": [
                    "subsystem bar chart",
                    "source-split rejection table",
                    "post-hoc BDT closure metrics",
                ],
                "png_dimensions": [SLIDE_W, SLIDE_H],
            },
            indent=2,
            sort_keys=True,
        )
    )
    print(OUT_PNG)
    print(OUT_MD)
    print(OUT_MANIFEST)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Build alternate THE-32 low-calo explanation visuals."""

from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


W = 2560
H = 1440
OUTDIR = Path("dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603")
SOURCE_PNG = OUTDIR / "the32_low_calo_pathology_slide_v1.png"
SOURCE_JSON = OUTDIR / "the32_low_calo_pathology_slide_v1.json"

INK = (24, 34, 48)
MUTED = (71, 84, 103)
LIGHT_TEXT = (255, 255, 255)
BLUE = (32, 112, 181)
BLUE_LIGHT = (216, 235, 252)
RED = (184, 35, 24)
RED_LIGHT = (255, 230, 225)
GREEN = (20, 133, 92)
GREEN_LIGHT = (225, 245, 237)
YELLOW = (255, 247, 214)
YELLOW_BORDER = (227, 179, 65)
BORDER = (198, 208, 222)
GRAY = (246, 248, 251)


def font(size: int, *, bold: bool = False, italic: bool = False) -> ImageFont.FreeTypeFont:
    names = []
    if bold and italic:
        names = ["Times New Roman Bold Italic.ttf", "Arial Bold Italic.ttf"]
    elif bold:
        names = ["Times New Roman Bold.ttf", "Arial Bold.ttf"]
    elif italic:
        names = ["Times New Roman Italic.ttf", "Arial Italic.ttf"]
    else:
        names = ["Times New Roman.ttf", "Arial.ttf"]
    for name in names:
        path = Path("/System/Library/Fonts/Supplemental") / name
        if path.exists():
            return ImageFont.truetype(str(path), size=size)
    return ImageFont.load_default(size=size)


F_TITLE = font(60, bold=True)
F_SUB = font(34)
F_H1 = font(46, bold=True)
F_H2 = font(38, bold=True)
F_BODY = font(32)
F_BODY_BOLD = font(32, bold=True)
F_SMALL = font(27)
F_SMALL_BOLD = font(27, bold=True)
F_TINY = font(22)


def draw_wrapped(draw: ImageDraw.ImageDraw, xy, text: str, *, max_width: int, text_font, fill=INK, spacing=9) -> int:
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
        y += text_font.size + spacing
    return y


def box(draw, xy, *, fill, outline=BORDER, width=3, radius=22):
    draw.rounded_rectangle(xy, radius=radius, fill=fill, outline=outline, width=width)


def title(draw, main: str, sub: str) -> None:
    draw.text((100, 70), main, font=F_TITLE, fill=INK)
    draw_wrapped(draw, (104, 155), sub, max_width=2200, text_font=F_SUB, fill=MUTED, spacing=7)


def load_counts() -> tuple[list[dict], int, int]:
    payload = json.loads(SOURCE_JSON.read_text())
    rows = payload["event_counts"]
    total = sum(int(r["event_total"]) for r in rows)
    rejected = sum(int(r["event_rejected"]) for r in rows)
    return rows, total, rejected


def aggregate_by_cent(rows: list[dict]) -> list[tuple[str, int, int]]:
    out = []
    for cent in ["0-20%", "20-50%", "50-80%"]:
        subset = [r for r in rows if r["centrality_bin"] == cent]
        total = sum(int(r["event_total"]) for r in subset)
        rejected = sum(int(r["event_rejected"]) for r in subset)
        out.append((cent, total, rejected))
    return out


def save(canvas: Image.Image, name: str) -> Path:
    path = OUTDIR / name
    canvas.save(path)
    return path


def variant_schematic(rows: list[dict], total: int, rejected: int) -> Path:
    im = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(im)
    title(
        d,
        "Most events behave normally; the cut removes the obvious mismatch",
        "Use this as the clean audience explanation: central events should not look calorimeter-empty.",
    )

    box(d, (145, 275, 1185, 1040), fill=GREEN_LIGHT, outline=(143, 210, 181))
    box(d, (1375, 275, 2415, 1040), fill=RED_LIGHT, outline=(232, 145, 137))
    d.text((220, 330), "Normal event", font=F_H1, fill=GREEN)
    d.text((1450, 330), "Event removed by cut", font=F_H1, fill=RED)

    # Simple detector stacks.
    for x0, color in [(320, BLUE), (1550, RED)]:
        d.ellipse((x0 + 110, 470, x0 + 470, 830), outline=color, width=10)
        d.ellipse((x0 + 40, 400, x0 + 540, 900), outline=color, width=8)
        d.rectangle((x0 + 250, 380, x0 + 330, 920), fill=color)
    d.text((310, 930), "Centrality says: central event", font=F_BODY_BOLD, fill=INK)
    d.text((360, 982), "Calo energy also looks central", font=F_BODY, fill=MUTED)
    d.text((1530, 930), "Centrality says: central event", font=F_BODY_BOLD, fill=INK)
    d.text((1548, 982), "But calo energy looks too low", font=F_BODY, fill=MUTED)

    d.line((1195, 660, 1365, 660), fill=INK, width=5)
    d.polygon([(1365, 660), (1335, 642), (1335, 678)], fill=INK)
    d.text((1137, 600), "compare", font=F_SMALL_BOLD, fill=INK)

    box(d, (170, 1135, 2390, 1318), fill=YELLOW, outline=YELLOW_BORDER)
    d.text((220, 1175), "Clean claim:", font=F_H2, fill=INK)
    draw_wrapped(
        d,
        (505, 1179),
        f"The diagnostic finds {rejected:,}/{total:,} events ({rejected/total:.1%}) where centrality and total calorimeter energy disagree. The upstream cut removes those mismatched events before BDT training.",
        max_width=1850,
        text_font=F_BODY,
        fill=INK,
        spacing=8,
    )
    d.text((170, 1346), "Schematic view; count comes from the THE-32 diagnostic sample.", font=F_TINY, fill=MUTED)
    return save(im, "the32_low_calo_variant_a_schematic_mismatch_v1.png")


def variant_centrality_summary(rows: list[dict], total: int, rejected: int) -> Path:
    im = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(im)
    title(
        d,
        "Where the low-calo issue shows up",
        "This version avoids the 2D density plot: it only shows how much of each centrality region is removed.",
    )
    aggs = aggregate_by_cent(rows)
    x0, y0 = 240, 380
    bar_w, bar_h = 560, 520
    gap = 120
    scale_max = 0.10
    for i, (label, tot, rej) in enumerate(aggs):
        x = x0 + i * (bar_w + gap)
        frac = rej / tot
        box(d, (x, y0, x + bar_w, y0 + bar_h), fill=GRAY, outline=BORDER)
        d.text((x + 45, y0 + 35), label, font=F_H1, fill=INK)
        d.text((x + bar_w - 210, y0 + 48), "0-10% scale", font=F_TINY, fill=MUTED)
        d.text((x + 45, y0 + 105), "events removed", font=F_BODY, fill=MUTED)
        slot_top = y0 + 185
        red_h = int((bar_h - 210) * min(frac / scale_max, 1.0))
        base = y0 + bar_h - 80
        red_top = max(slot_top, base - red_h)
        if red_top > slot_top:
            d.rectangle((x + 90, slot_top, x + bar_w - 90, red_top), fill=BLUE_LIGHT)
        d.rectangle((x + 90, red_top, x + bar_w - 90, base), fill=RED)
        d.rectangle((x + 90, slot_top, x + bar_w - 90, base), outline=INK, width=3)
        d.text((x + 170, base - red_h - 72), f"{frac:.1%}", font=F_H1, fill=RED)
        d.text((x + 95, base + 25), f"{rej:,} of {tot:,}", font=F_SMALL, fill=INK)

    box(d, (250, 1070, 2310, 1260), fill=YELLOW, outline=YELLOW_BORDER)
    d.text((300, 1114), "Interpretation:", font=F_H2, fill=INK)
    draw_wrapped(
        d,
        (575, 1118),
        "The issue is concentrated in the more central part of the sample. In the peripheral 50-80% bin, almost nothing is removed, which is what you want from a targeted event-quality cut.",
        max_width=1700,
        text_font=F_BODY,
        fill=INK,
    )
    d.text((250, 1315), f"Total diagnostic sample: {rejected:,}/{total:,} events removed ({rejected/total:.1%}).", font=F_SMALL, fill=MUTED)
    return save(im, "the32_low_calo_variant_b_centrality_summary_v1.png")


def variant_scale(rows: list[dict], total: int, rejected: int) -> Path:
    im = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(im)
    kept = total - rejected
    title(
        d,
        "The cut is targeted, not a broad reshaping of the sample",
        "This view answers the scale question first: almost all events stay; only the low-calo event-quality tail is removed.",
    )
    box(d, (210, 405, 2350, 710), fill=GRAY, outline=BORDER)
    margin = 80
    bar_x0, bar_x1 = 320, 2240
    bar_y0, bar_y1 = 525, 615
    d.rectangle((bar_x0, bar_y0, bar_x1, bar_y1), fill=GREEN)
    red_w = int((bar_x1 - bar_x0) * rejected / total)
    d.rectangle((bar_x1 - red_w, bar_y0, bar_x1, bar_y1), fill=RED)
    d.text((bar_x0, 445), f"Kept: {kept:,} events", font=F_H2, fill=GREEN)
    d.text((bar_x1 - 610, 445), f"Removed: {rejected:,}", font=F_H2, fill=RED)
    d.text((bar_x0, 650), "96.6% kept", font=F_BODY_BOLD, fill=GREEN)
    d.text((bar_x1 - 270, 650), "3.4% removed", font=F_BODY_BOLD, fill=RED)

    box(d, (260, 835, 1135, 1115), fill=GREEN_LIGHT, outline=(143, 210, 181))
    box(d, (1395, 835, 2305, 1115), fill=RED_LIGHT, outline=(232, 145, 137))
    d.text((330, 890), "Kept events", font=F_H1, fill=GREEN)
    draw_wrapped(d, (330, 960), "Consistent centrality and total calorimeter energy.", max_width=700, text_font=F_BODY, fill=INK)
    d.text((1465, 890), "Removed events", font=F_H1, fill=RED)
    draw_wrapped(d, (1465, 960), "Too little total calorimeter energy for the assigned centrality.", max_width=730, text_font=F_BODY, fill=INK)

    d.text((285, 1228), "Best use:", font=F_H2, fill=INK)
    draw_wrapped(
        d,
        (500, 1232),
        "Use this if the audience is worried that the cut throws away too much. It shows the cut is small and targeted.",
        max_width=1700,
        text_font=F_BODY,
        fill=INK,
    )
    return save(im, "the32_low_calo_variant_c_targeted_scale_v1.png")


def variant_minimal_evidence(rows: list[dict], total: int, rejected: int) -> Path:
    source = Image.open(SOURCE_PNG).convert("RGB")
    crop = source.crop((15, 300, 1710, 1060))
    im = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(im)
    title(
        d,
        "Data view: red events are the low-calo tail",
        "This is closest to the original evidence plot, but with the side plots removed.",
    )
    target_w = 2100
    target_h = int(crop.height * target_w / crop.width)
    big = crop.resize((target_w, target_h), Image.Resampling.LANCZOS)
    px = (W - target_w) // 2
    py = 240
    box(d, (px - 18, py - 18, px + target_w + 18, py + target_h + 18), fill=GRAY, outline=BORDER)
    im.paste(big, (px, py))
    box(d, (150, 1220, 2410, 1345), fill=YELLOW, outline=YELLOW_BORDER)
    d.text((200, 1248), "One-sentence read:", font=F_H2, fill=INK)
    draw_wrapped(
        d,
        (575, 1253),
        f"At the same centrality, the red events sit below the normal total-calo band, so the cut removes {rejected:,}/{total:,} event-quality outliers before training.",
        max_width=1800,
        text_font=F_BODY,
        fill=INK,
    )
    return save(im, "the32_low_calo_variant_d_minimal_evidence_plot_v1.png")


def main() -> int:
    rows, total, rejected = load_counts()
    outputs = [
        variant_schematic(rows, total, rejected),
        variant_centrality_summary(rows, total, rejected),
        variant_scale(rows, total, rejected),
        variant_minimal_evidence(rows, total, rejected),
    ]
    note = OUTDIR / "the32_low_calo_explanation_variants_v1.md"
    note.write_text(
        "\n".join(
            [
                "# THE-32 low-calo explanation variants",
                "",
                "Variant A is the clearest no-axis schematic: centrality says central, calo energy says too low.",
                "Variant B summarizes the cut by centrality bin and shows the issue is central-region dominated.",
                "Variant C is a scale check: almost everything is kept, only the low-calo tail is removed.",
                "Variant D keeps the original evidence plot but removes the extra bar/table panels.",
                "",
            ]
        )
    )
    manifest = OUTDIR / "the32_low_calo_explanation_variants_v1.json"
    manifest.write_text(
        json.dumps(
            {
                "schema": "THE32_LOW_CALO_EXPLANATION_VARIANTS_V1",
                "source_json": str(SOURCE_JSON),
                "source_png": str(SOURCE_PNG),
                "event_total": total,
                "event_rejected": rejected,
                "event_rejected_fraction": rejected / total,
                "outputs": [str(p) for p in outputs],
                "speaker_note": str(note),
                "google_slides_mutated": False,
            },
            indent=2,
            sort_keys=True,
        )
    )
    for out in outputs:
        print(out)
    print(note)
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Overlay signal/background counts on the Stack B slide-28 PNG."""

from __future__ import annotations

import csv
import shutil
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


REPO = Path(__file__).resolve().parents[1]
SLIDE = (
    REPO
    / "dataOutput/auauMLDiagnosticRuns/ppg12_weighted_basev3e_stack_20260521_231600/"
    "slideReady/remade_slides/slide25_score_separation_stackB_replacement_v9.png"
)
VERSIONED = SLIDE.with_name("slide25_score_separation_stackB_replacement_v12_counts.png")
BACKUP = SLIDE.with_name("slide25_score_separation_stackB_replacement_v9_before_counts.png")
SUMMARY = (
    REPO
    / "dataOutput/auauMLDiagnosticRuns/ppg12_weighted_basev3e_stack_20260521_231600/"
    "slideReady/score_separation_3x3_fullstat/stack_b_score_et_cent/score_separation_3x3_summary.csv"
)

PANEL_X = [591, 995, 1398]
PANEL_Y = [236, 450, 663]
CENTERS = [(0.0, 20.0), (20.0, 50.0), (50.0, 80.0)]
MODELS = [
    "Input BDT, 8 E_T x 7 centrality",
    "Input MLP, single small NN",
    "Stack B: scores + E_T + centrality",
]
TEXT = "#4B5563"
BOX_FILL = (255, 255, 255, 236)
BOX_EDGE = "#E5E7EB"
SIGNAL = "#0B6FB3"
BACKGROUND = "#C43B54"
PINK_FILL = "#F7E7F1"
ORANGE_FILL = "#F9EFDF"
GREEN_FILL = "#DDEFE8"
PINK_EDGE = "#C04DA2"
ORANGE_EDGE = "#E6862C"
GREEN_EDGE = "#19885A"
ROW_BOXES = [
    (535, 228, 1840, 444, PINK_EDGE, PINK_FILL),
    (535, 452, 1840, 654, ORANGE_EDGE, ORANGE_FILL),
    (535, 663, 1840, 918, GREEN_EDGE, GREEN_FILL),
]
ET = "{ET}"


def font(size: int, *, bold: bool = False) -> ImageFont.FreeTypeFont:
    names = [
        "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
        "/Library/Fonts/Times New Roman Bold.ttf" if bold else "/Library/Fonts/Times New Roman.ttf",
        "/System/Library/Fonts/Times.ttc",
    ]
    for name in names:
        path = Path(name)
        if path.exists():
            return ImageFont.truetype(str(path), size=size)
    return ImageFont.load_default()


def compact_count(value: int) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f}M"
    if value >= 10_000:
        return f"{value / 1_000:.0f}k"
    return f"{value:,}"


def advance(draw: ImageDraw.ImageDraw, text: str, text_font: ImageFont.FreeTypeFont) -> int:
    if not text:
        return 0
    bbox = draw.textbbox((0, 0), text, font=text_font)
    return bbox[2] - bbox[0]


def draw_text_with_et(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    *,
    text_font: ImageFont.FreeTypeFont,
    fill: str,
    sub_font: ImageFont.FreeTypeFont | None = None,
) -> None:
    x, y = xy
    sub = sub_font or font(max(12, int(getattr(text_font, "size", 22) * 0.64)))
    parts = text.split(ET)
    for idx, part in enumerate(parts):
        if part:
            draw.text((x, y), part, font=text_font, fill=fill)
            x += advance(draw, part, text_font)
        if idx < len(parts) - 1:
            draw.text((x, y), "E", font=text_font, fill=fill)
            x += advance(draw, "E", text_font) - 1
            draw.text((x, y + int(getattr(text_font, "size", 22) * 0.34)), "T", font=sub, fill=fill)
            x += advance(draw, "T", sub) + 2


def read_counts() -> dict[tuple[str, tuple[float, float]], tuple[int, int]]:
    counts: dict[tuple[str, tuple[float, float]], tuple[int, int]] = {}
    with SUMMARY.open() as handle:
        for row in csv.DictReader(handle):
            key = (row["model"], (float(row["cent_lo"]), float(row["cent_hi"])))
            counts[key] = (int(row["signal_entries"]), int(row["background_entries"]))
    return counts


def redraw_top_left_legend(draw: ImageDraw.ImageDraw) -> None:
    label_font = font(26)
    # Replace the old pills with slightly wider, fixed-width legend chips.
    draw.rectangle((66, 122, 455, 174), fill=(255, 255, 255, 255))
    signal_box = (74, 130, 210, 166)
    bkg_box = (226, 130, 426, 166)
    draw.rounded_rectangle(signal_box, radius=12, fill=(248, 251, 255, 255), outline="#CBD5E1", width=1)
    draw.rounded_rectangle(bkg_box, radius=12, fill=(255, 247, 249, 255), outline="#E9C2CB", width=1)
    draw.line((95, 148, 137, 148), fill=SIGNAL, width=5)
    draw.text((149, 135), "signal", font=label_font, fill="#111827")
    draw.line((250, 148, 292, 148), fill=BACKGROUND, width=5)
    draw.text((304, 135), "background", font=label_font, fill="#111827")


def redraw_model_cards(draw: ImageDraw.ImageDraw) -> None:
    cards = [
        {
            "box": (72, 278, 447, 421),
            "fill": PINK_FILL,
            "edge": PINK_EDGE,
            "title": "Input BDT",
            "line1": "global BDT score",
            "line2": "baseline score shape",
            "auc": "AUC: 0.795  0.816  0.830",
            "accent": "#111827",
        },
        {
            "box": (72, 447, 447, 590),
            "fill": ORANGE_FILL,
            "edge": ORANGE_EDGE,
            "title": "Input MLP",
            "line1": "small NN score",
            "line2": "same inputs/test rows",
            "auc": "AUC: 0.803  0.827  0.842",
            "accent": "#111827",
        },
        {
            "box": (72, 615, 447, 758),
            "fill": GREEN_FILL,
            "edge": GREEN_EDGE,
            "title": "Stack B",
            "line1": f"score pair + {ET} + centrality",
            "line2": "compact stack shown below",
            "auc": "AUC: 0.809  0.832  0.846",
            "accent": GREEN_EDGE,
        },
    ]
    title_font = font(28, bold=True)
    body_font = font(25)
    small_font = font(21)
    note_font = font(21)
    for card in cards:
        x0, y0, x1, y1 = card["box"]
        draw.rounded_rectangle((x0, y0, x1, y1), radius=12, fill=card["fill"], outline=None)
        draw.rounded_rectangle((x0, y0, x0 + 8, y1), radius=4, fill=card["edge"], outline=None)
        title_color = card.get("accent", "#111827")
        draw.text((x0 + 24, y0 + 20), card["title"], font=title_font, fill=title_color)
        draw_text_with_et(
            draw,
            (x0 + 24, y0 + 54),
            card["line1"],
            text_font=body_font,
            fill="#111827",
        )
        draw.text((x0 + 24, y0 + 82), card["line2"], font=small_font, fill="#111827")
        draw.text((x0 + 24, y0 + 111), card["auc"], font=body_font, fill="#111827")


def redraw_stack_variant_choice(draw: ImageDraw.ImageDraw) -> None:
    title_font = font(25, bold=True)
    body_font = font(21)
    body_bold = font(21, bold=True)
    x0, y0, x1, y1 = (72, 792, 529, 902)
    draw.rounded_rectangle((x0, y0, x1, y1), radius=16, fill=(246, 249, 248, 255), outline="#D1D5DB", width=1)
    draw.text((92, 813), "Stack variant choice", font=title_font, fill="#111827")
    draw.text((92, 846), "A: scores + full features", font=body_font, fill="#4B5563")
    draw.text((322, 846), "0.803 / 0.829 / 0.844", font=body_font, fill="#4B5563")
    draw_text_with_et(
        draw,
        (92, 874),
        f"B: scores + {ET} + cent",
        text_font=body_bold,
        fill=GREEN_EDGE,
    )
    draw.text((322, 874), "0.809 / 0.832 / 0.846", font=body_bold, fill=GREEN_EDGE)


def redraw_plot_row_boxes(draw: ImageDraw.ImageDraw) -> None:
    # Reinforce the row grouping without drawing border strokes over the axes.
    # The base slide already has row-color bands behind the plots; these bars
    # and gap separators make the table structure clear while leaving every
    # embedded plot fully inside its white panel.
    for x0, y0, x1, y1, edge, fill in ROW_BOXES:
        draw.rounded_rectangle((x0, y0, x0 + 15, y1), radius=8, fill=edge, outline=None)
        draw.rounded_rectangle((x1 - 8, y0 + 10, x1, y1 - 10), radius=5, fill=fill, outline=None)

    separators = [
        (548, 448, 1830, 452, PINK_EDGE),
        (548, 656, 1830, 660, ORANGE_EDGE),
    ]
    for x0, y0, x1, y1, color in separators:
        draw.rounded_rectangle((x0, y0, x1, y1), radius=4, fill=color, outline=None)


def redraw_readout(draw: ImageDraw.ImageDraw) -> None:
    title_font = font(26, bold=True)
    body_font = font(25)
    body_green = font(25)
    footer_font = font(21)

    draw.rounded_rectangle((72, 930, 1850, 1020), radius=16, fill=(248, 250, 252, 255), outline=None)
    draw.text((102, 956), "Readout", font=title_font, fill="#111827")
    draw.text(
        (206, 955),
        "Stack B keeps the signal peak at high score while concentrating background near zero.",
        font=body_font,
        fill="#111827",
    )
    draw_text_with_et(
        draw,
        (206, 987),
        f"Full-stat validation: the compact stack is best in each centrality bin using scores + {ET} + centrality.",
        text_font=body_green,
        fill=GREEN_EDGE,
    )
    draw.rectangle((68, 1030, 620, 1070), fill=(255, 255, 255, 255))
    draw_text_with_et(
        draw,
        (76, 1042),
        f"Validation uses held-out test split, 15 < {ET} < 35 GeV.",
        text_font=footer_font,
        fill="#6B7280",
    )


def main() -> None:
    if not BACKUP.exists():
        shutil.copy2(SLIDE, BACKUP)

    counts = read_counts()
    image = Image.open(BACKUP).convert("RGBA")
    draw = ImageDraw.Draw(image, "RGBA")
    text_font = font(18)
    redraw_top_left_legend(draw)
    redraw_model_cards(draw)
    redraw_stack_variant_choice(draw)
    redraw_plot_row_boxes(draw)
    redraw_readout(draw)

    for row_idx, model in enumerate(MODELS):
        for col_idx, cent in enumerate(CENTERS):
            signal, background = counts[(model, cent)]
            label = f"S {compact_count(signal)}  B {compact_count(background)}"
            x = PANEL_X[col_idx] + 16
            y = PANEL_Y[row_idx] + 47
            bbox = draw.textbbox((x, y), label, font=text_font)
            pad_x = 5
            pad_y = 3
            box = (bbox[0] - pad_x, bbox[1] - pad_y, bbox[2] + pad_x, bbox[3] + pad_y)
            draw.rounded_rectangle(box, radius=4, fill=BOX_FILL, outline=BOX_EDGE, width=1)
            draw.text((x, y), label, font=text_font, fill=TEXT)

    image.convert("RGB").save(SLIDE)
    image.convert("RGB").save(VERSIONED)
    print(SLIDE)
    print(VERSIONED)


if __name__ == "__main__":
    main()

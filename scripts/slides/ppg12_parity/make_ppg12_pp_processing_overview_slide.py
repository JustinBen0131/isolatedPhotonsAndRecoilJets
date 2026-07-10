#!/usr/bin/env python3
"""Generate a full-slide PNG for the canonical pp processing overview."""

from __future__ import annotations

import json
from pathlib import Path
from textwrap import wrap

from PIL import Image, ImageDraw, ImageFont


REPO = Path(__file__).resolve().parents[3]
OUT_DIR = (
    REPO
    / "dataOutput"
    / "ppg12Parity"
    / "control_plane"
    / "slide_candidates"
    / "20260706_pp_processing_overview"
)
PNG_PATH = OUT_DIR / "ppg12_pp_processing_overview.png"
SCRIPT_PATH = OUT_DIR / "ppg12_pp_processing_overview_speaker_script.md"
MANIFEST_PATH = OUT_DIR / "ppg12_pp_processing_overview_manifest.json"

W, H = 2560, 1440

FONT_REG = "/System/Library/Fonts/Supplemental/Times New Roman.ttf"
FONT_BOLD = "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"
FONT_ITALIC = "/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf"


def font(path: str, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(path, size=size)


F_TITLE = font(FONT_BOLD, 82)
F_CARD_TITLE = font(FONT_BOLD, 46)
F_LABEL = font(FONT_BOLD, 34)
F_BODY = font(FONT_REG, 34)
F_BODY_SMALL = font(FONT_REG, 30)
F_BODY_TINY = font(FONT_REG, 28)
F_FORMULA = font(FONT_REG, 40)
F_FORMULA_BOLD = font(FONT_BOLD, 42)
F_CHIP = font(FONT_BOLD, 30)
F_NOTE = font(FONT_REG, 36)


def text_size(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.FreeTypeFont) -> tuple[int, int]:
    box = draw.textbbox((0, 0), text, font=fnt)
    return box[2] - box[0], box[3] - box[1]


def centered_text(
    draw: ImageDraw.ImageDraw,
    xy: tuple[float, float],
    text: str,
    fnt: ImageFont.FreeTypeFont,
    fill: str,
) -> None:
    tw, th = text_size(draw, text, fnt)
    draw.text((xy[0] - tw / 2, xy[1] - th / 2), text, font=fnt, fill=fill)


def wrap_lines(text: str, width: int) -> list[str]:
    return wrap(text, width=width, break_long_words=False, break_on_hyphens=False)


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    fnt: ImageFont.FreeTypeFont,
    fill: str,
    width: int,
    line_gap: int = 8,
) -> int:
    x, y = xy
    for line in wrap_lines(text, width):
        draw.text((x, y), line, font=fnt, fill=fill)
        _, lh = text_size(draw, line, fnt)
        y += lh + line_gap
    return y


def rounded_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    outline: str,
    accent: str,
    title: str,
) -> None:
    x1, y1, x2, y2 = box
    draw.rounded_rectangle(box, radius=26, fill="#ffffff", outline=outline, width=3)
    draw.rounded_rectangle((x1, y1, x1 + 18, y2), radius=12, fill=accent)
    draw.text((x1 + 54, y1 + 34), title, font=F_CARD_TITLE, fill="#111111")


def chip(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    label: str,
    fill: str,
    outline: str | None = None,
    text_fill: str = "#111111",
) -> None:
    draw.rounded_rectangle(box, radius=18, fill=fill, outline=outline or fill, width=2)
    centered_text(draw, ((box[0] + box[2]) / 2, (box[1] + box[3]) / 2), label, F_CHIP, text_fill)


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], color: str) -> None:
    draw.line((start, end), fill=color, width=6)
    ex, ey = end
    sx, sy = start
    if abs(ex - sx) >= abs(ey - sy):
        direction = 1 if ex > sx else -1
        pts = [(ex, ey), (ex - direction * 28, ey - 16), (ex - direction * 28, ey + 16)]
    else:
        direction = 1 if ey > sy else -1
        pts = [(ex, ey), (ex - 16, ey - direction * 28), (ex + 16, ey - direction * 28)]
    draw.polygon(pts, fill=color)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    img = Image.new("RGB", (W, H), "#ffffff")
    draw = ImageDraw.Draw(img)

    # Title axis and top band.
    title_x = 116
    title_y = 82
    draw.text((title_x, title_y), "Canonical pp processing overview", font=F_TITLE, fill="#090909")

    card_y1, card_y2 = 240, 1085
    gap = 48
    card_w = 725
    x1 = title_x
    x2 = x1 + card_w + gap
    x3 = x2 + card_w + gap

    card1 = (x1, card_y1, x1 + card_w, card_y2)
    card2 = (x2, card_y1, x2 + card_w, card_y2)
    card3 = (x3, card_y1, x3 + card_w, card_y2)

    rounded_card(draw, card1, "#c9d6df", "#2b7bba", "Data periods")
    rounded_card(draw, card2, "#cfd8cf", "#2f9d55", "SIM lanes")
    rounded_card(draw, card3, "#dacfe5", "#c03a91", "Event weight")

    # Data card.
    c1x, c1y = card1[0], card1[1]
    draw_wrapped(
        draw,
        (c1x + 58, c1y + 120),
        "Keep pp data separated by running period before the final weighted comparison.",
        F_BODY,
        "#202020",
        35,
        10,
    )
    chip(draw, (c1x + 92, c1y + 295, c1x + 342, c1y + 365), "0 mrad", "#dcecff", "#2b7bba")
    chip(draw, (c1x + 388, c1y + 295, c1x + 638, c1y + 365), "1.5 mrad", "#dcecff", "#2b7bba")
    arrow(draw, (c1x + 365, c1y + 395), (c1x + 365, c1y + 490), "#5f6d78")
    draw_wrapped(
        draw,
        (c1x + 72, c1y + 520),
        "The period split is the anchor for the matching SIM vertex reference and exposure bookkeeping.",
        F_BODY_SMALL,
        "#202020",
        35,
        9,
    )
    draw.rounded_rectangle((c1x + 70, c1y + 695, c1x + 655, c1y + 792), radius=18, fill="#ffffff", outline="#8fb9df", width=3)
    draw_wrapped(
        draw,
        (c1x + 104, c1y + 718),
        "Analysis output keeps this period identity until the final merge.",
        F_BODY_TINY,
        "#202020",
        42,
        6,
    )

    # SIM card.
    c2x, c2y = card2[0], card2[1]
    # SIM lane card: three intentionally spaced tiers.
    tier_x1, tier_x2 = c2x + 62, c2x + 663
    draw.rounded_rectangle((tier_x1, c2y + 112, tier_x2, c2y + 270), radius=20, fill="#ffffff", outline="#d8e6db", width=2)
    draw.text((c2x + 86, c2y + 132), "period", font=F_LABEL, fill="#111111")
    chip(draw, (c2x + 290, c2y + 126, c2x + 455, c2y + 186), "0 mrad", "#e7f4ea", "#2f9d55")
    chip(draw, (c2x + 480, c2y + 126, c2x + 638, c2y + 186), "1.5 mrad", "#e7f4ea", "#2f9d55")

    draw.rounded_rectangle((tier_x1, c2y + 302, tier_x2, c2y + 460), radius=20, fill="#ffffff", outline="#eadcf2", width=2)
    draw.text((c2x + 86, c2y + 322), "interaction", font=F_LABEL, fill="#111111")
    chip(draw, (c2x + 305, c2y + 316, c2x + 457, c2y + 376), "SI", "#f3ebf9", "#9b54bc")
    chip(draw, (c2x + 486, c2y + 316, c2x + 638, c2y + 376), "DI", "#f3ebf9", "#9b54bc")

    draw.rounded_rectangle((tier_x1, c2y + 492, tier_x2, c2y + 720), radius=20, fill="#ffffff", outline="#f1d8aa", width=2)
    draw.text((c2x + 86, c2y + 512), "hard process", font=F_LABEL, fill="#111111")
    chip(draw, (c2x + 112, c2y + 582, c2x + 626, c2y + 642), "photon5 / 10 / 20", "#fff1d8", "#d58a1f")
    chip(draw, (c2x + 112, c2y + 656, c2x + 626, c2y + 716), "jet8 / 12 / 20 / 30 / 40", "#fff1d8", "#d58a1f")

    draw.rounded_rectangle((c2x + 72, c2y + 752, c2x + 650, c2y + 816), radius=18, fill="#ffffff", outline="#8eb69a", width=3)
    centered_text(draw, (c2x + card_w / 2, c2y + 783), "processed independently; summed once", F_BODY_TINY, "#202020")

    # Weight card.
    c3x, c3y = card3[0], card3[1]
    centered_text(draw, (c3x + card_w / 2, c3y + 155), "event weight =", F_FORMULA_BOLD, "#111111")
    weight_rows = [
        ("generator cross-section", "#2b7bba"),
        ("× period exposure", "#2f9d55"),
        ("× vertex reweight", "#c03a91"),
        ("× SI/DI mixture", "#d58a1f"),
    ]
    yy = c3y + 240
    for text, color in weight_rows:
        draw.rounded_rectangle((c3x + 80, yy, c3x + 645, yy + 78), radius=18, fill="#ffffff", outline=color, width=4)
        centered_text(draw, (c3x + card_w / 2, yy + 39), text, F_FORMULA, color)
        yy += 108
    draw_wrapped(
        draw,
        (c3x + 72, c3y + 690),
        "The merge step is a histogram sum, not a later shape rescale.",
        F_BODY_SMALL,
        "#202020",
        35,
        8,
    )
    draw.rounded_rectangle((c3x + 76, c3y + 758, c3x + 642, c3y + 828), radius=18, fill="#ffffff", outline="#ccb2d7", width=3)
    centered_text(draw, (c3x + card_w / 2, c3y + 792), "weights enter before the merge", F_BODY_TINY, "#202020")

    # Arrows between cards.
    arrow(draw, (card1[2] + 10, (card_y1 + card_y2) // 2), (card2[0] - 18, (card_y1 + card_y2) // 2), "#b7b7b7")
    arrow(draw, (card2[2] + 10, (card_y1 + card_y2) // 2), (card3[0] - 18, (card_y1 + card_y2) // 2), "#b7b7b7")

    # Bottom contract band.
    band = (title_x, 1156, W - title_x, 1324)
    draw.rounded_rectangle(band, radius=24, fill="#ffffff", outline="#2f2f2f", width=3)
    draw.rectangle((band[0], band[1], band[0] + 16, band[3]), fill="#111111")
    draw.text((band[0] + 54, band[1] + 35), "Canonical merge rule", font=F_LABEL, fill="#111111")
    draw_wrapped(
        draw,
        (band[0] + 430, band[1] + 31),
        "Split first, weight at the event/component level, then merge. Any downstream pp combined result should follow this contract unless it is explicitly labeled as a diagnostic exception.",
        F_NOTE,
        "#161616",
        92,
        6,
    )

    img.save(PNG_PATH)

    SCRIPT_PATH.write_text(
        """# Speaker Script: Canonical pp Processing Overview

This slide is the processing contract we want to use for pp going forward.

The important point is that we do not collapse everything into one undifferentiated sample at the start. Data stays split by running period, and simulation is split by period, by single- versus double-interaction component, and by hard-process slice.

Each simulated event then carries the full event weight before the merge: the generator cross-section weight, the period or luminosity exposure, the period-specific vertex reweighting, and the single/double-interaction mixture weight.

After that, the merge is just a sum of already-weighted component histograms. This is the policy that lets downstream photon+jet and inclusive-jet pp products be compared to PPG12 in a controlled way, instead of relying on ad hoc rescaling later.
""",
        encoding="utf-8",
    )

    MANIFEST_PATH.write_text(
        json.dumps(
            {
                "artifact": str(PNG_PATH),
                "script": str(SCRIPT_PATH),
                "width_px": W,
                "height_px": H,
                "slide_policy": "full-slide PNG candidate, no Google Slides mutation",
                "title": "Canonical pp processing overview",
                "notes": [
                    "Opening processing overview only.",
                    "Does not include the Fig.5 source-stage diagnostic caveat by user request.",
                    "No slide number or provenance footer baked into the PNG.",
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    print(PNG_PATH)
    print(SCRIPT_PATH)
    print(MANIFEST_PATH)


if __name__ == "__main__":
    main()

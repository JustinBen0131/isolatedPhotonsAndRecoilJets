#!/usr/bin/env python3
"""Make a full-slide clean-run strategy candidate for scaled-trigger QA."""

from __future__ import annotations

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


BASE = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auau/scaledTriggerRunByRunQA/"
    "scaled_trigger_run_by_run_20260521_202504"
)
PLOT_PATH = BASE / "clean_efficiency_turnon_tables" / "combined_clean63_scaled_trigger_1x2.png"
OUT_DIR = BASE / "slide_candidates"
OUT_PNG = OUT_DIR / "clean_full_efficiency_combined_strategy_after_slide3_candidate.png"

FONT_REG = "/System/Library/Fonts/Supplemental/Times New Roman.ttf"
FONT_BOLD = "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"
FONT_ITALIC = "/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf"


def font(path: str, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(path, size=size)


W, H = 3840, 2160
INK = (18, 24, 32)
MUTED = (82, 92, 105)
LINE = (188, 198, 207)
GREEN = (22, 105, 67)
GREEN_FILL = (235, 247, 240)
BLUE_FILL = (238, 244, 252)
GRAY_FILL = (247, 249, 251)

F_TITLE = font(FONT_BOLD, 76)
F_SUB = font(FONT_REG, 38)
F_HEAD = font(FONT_BOLD, 54)
F_BODY = font(FONT_REG, 44)
F_BODY_BOLD = font(FONT_BOLD, 44)
F_SPHX = font(FONT_BOLD, 36)
F_ITAL = font(FONT_ITALIC, 36)
F_SUBSCRIPT = font(FONT_ITALIC, 24)


def text_size(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.FreeTypeFont) -> tuple[int, int]:
    box = draw.textbbox((0, 0), text, font=fnt)
    return box[2] - box[0], box[3] - box[1]


def wrap_text(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.FreeTypeFont, max_width: int) -> list[str]:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        trial = word if not current else f"{current} {word}"
        if text_size(draw, trial, fnt)[0] <= max_width:
            current = trial
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    return lines


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    text: str,
    xy: tuple[int, int],
    fnt: ImageFont.FreeTypeFont,
    fill: tuple[int, int, int],
    max_width: int,
    line_gap: int = 8,
) -> int:
    x, y = xy
    for line in wrap_text(draw, text, fnt, max_width):
        draw.text((x, y), line, font=fnt, fill=fill)
        y += text_size(draw, line, fnt)[1] + line_gap
    return y


def draw_box(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int, int, int],
    title: str,
    body_lines: list[str],
    fill: tuple[int, int, int],
) -> None:
    x0, y0, x1, y1 = xy
    draw.rounded_rectangle(xy, radius=18, fill=fill, outline=LINE, width=3)
    draw.text((x0 + 30, y0 + 18), title, font=F_HEAD, fill=GREEN)
    y = y0 + 88
    for line in body_lines:
        if line.startswith("BOLD:"):
            y = draw_wrapped(draw, line[5:], (x0 + 34, y), F_BODY_BOLD, INK, x1 - x0 - 68, 4)
        else:
            y = draw_wrapped(draw, line, (x0 + 34, y), F_BODY, INK, x1 - x0 - 68, 4)
        y += 6


def draw_collision_label(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    prefix = "Internal  Au+Au, "
    draw.text((x, y), prefix, font=F_ITAL, fill=INK)
    px, _ = text_size(draw, prefix, F_ITAL)
    sx = x + px
    draw.text((sx, y), "√s", font=F_ITAL, fill=INK)
    root_w, _ = text_size(draw, "√s", F_ITAL)
    draw.text((sx + root_w - 2, y + 21), "NN", font=F_SUBSCRIPT, fill=INK)
    draw.text((sx + root_w + 33, y), "=200 GeV", font=F_ITAL, fill=INK)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    img = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(img)

    draw.text((76, 42), "Clean full-efficiency runs: combined QA and selection logic", font=F_TITLE, fill=INK)
    draw.text(
        (80, 132),
        "63/620 runs (10.2%) selected from run-by-run Trigger/MBD turn-on metrics; clean-run histograms are summed before forming ratios.",
        font=F_SUB,
        fill=MUTED,
    )
    draw.text((80, 184), "sPHENIX", font=F_SPHX, fill=INK)
    sx, _ = text_size(draw, "sPHENIX", F_SPHX)
    draw_collision_label(draw, 88 + sx, 184)

    plot = Image.open(PLOT_PATH).convert("RGB")
    # Keep the two scientific panels and their axis labels, but remove the
    # repeated title and the bottom stats box from the source diagnostic.
    crop = plot.crop((54, 220, 2350, 1235))
    plot_w = 3300
    plot_h = round(plot_w * crop.height / crop.width)
    crop = crop.resize((plot_w, plot_h), Image.Resampling.LANCZOS)
    plot_x = (W - crop.width) // 2
    plot_y = 250
    shadow = Image.new("RGBA", (crop.width + 24, crop.height + 24), (0, 0, 0, 0))
    shadow_draw = ImageDraw.Draw(shadow)
    shadow_draw.rounded_rectangle((14, 14, crop.width + 14, crop.height + 14), radius=16, fill=(0, 0, 0, 28))
    img.paste(shadow, (plot_x - 12, plot_y - 12), shadow)
    draw.rounded_rectangle((plot_x - 8, plot_y - 8, plot_x + crop.width + 8, plot_y + crop.height + 8), radius=16, fill="white", outline=LINE, width=3)
    img.paste(crop, (plot_x, plot_y))

    card_y0 = 1744
    card_y1 = 2090
    gap = 26
    card_w = (W - 2 * 76 - 2 * gap) // 3
    card1 = (76, card_y0, 76 + card_w, card_y1)
    card2 = (76 + card_w + gap, card_y0, 76 + 2 * card_w + gap, card_y1)
    card3 = (76 + 2 * (card_w + gap), card_y0, W - 76, card_y1)
    draw_box(
        draw,
        card1,
        "How Singled Out",
        [
            "BOLD:Pass = tail unity + low-E off.",
            "Tail: MBD >=100 and P10/MBD, P12/MBD in [0.97, 1.08].",
            "LHS overlay is QA only; RHS turn-on is the selector.",
        ],
        BLUE_FILL,
    )
    draw_box(
        draw,
        card2,
        "Low-E Veto",
        [
            "P10/P12 maxima by Emax band:",
            "1-3: .05/.02; 2-4: .10/.04.",
            "3-5: .25/.10; 4-6: .50/.22.",
            "Rejects 71260-like early turn-ons.",
        ],
        GRAY_FILL,
    )
    draw_box(
        draw,
        card3,
        "Result",
        [
            "BOLD:63/620 runs = 10.2% clean.",
            "Combined tail: P10/MBD=1.026, P12/MBD=1.020.",
            "Low 1-3 GeV: 0.0052/0.0012; mid 6-9: 0.746/0.567.",
        ],
        GREEN_FILL,
    )

    img.save(OUT_PNG, quality=95)
    print(OUT_PNG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

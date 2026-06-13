#!/usr/bin/env python3
"""Local Slide 14 split: main BDT-score focus plus shower-shape backup."""

from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = next(
    p
    for p in Path(__file__).resolve().parents
    if (p / "README.md").exists() and (p / "scripts").exists() and (p / "src").exists()
)
SCRIPT_DIR = ROOT / "scripts/slides/hp2026/fulltalk"
sys.path.insert(0, str(SCRIPT_DIR))

import make_hp2026_fulltalk_candidates as hp  # noqa: E402


OUTDIR = ROOT / "outputs/manual-20260610-slide14_bdt_focus_split"
MAIN_PNG = OUTDIR / "slide14_main_bdt_score_focus.png"
BACKUP_PNG = OUTDIR / "backup_slide14_shower_shape_inputs.png"
MANIFEST = OUTDIR / "slide14_bdt_focus_split_manifest.json"

HP2026_MAIN_HEADER = {
    "deck": "hp2026_main_talk",
    "title_font_size": 86,
    "subtitle_font_size": 56,
    "title_xy": [132, 76],
    "subtitle_xy": [136, 190],
    "divider_y": 286,
}


def draw_shell(
    title: str,
    subtitle: str | None = None,
    *,
    title_size: int = 86,
    title_xy: tuple[int, int] = (132, 76),
) -> Image.Image:
    img = Image.new("RGBA", (hp.W, hp.H), (*hp.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, hp.W, hp.H), fill=(*hp.SOFT_BG, 255))
    draw.rectangle((0, 0, hp.W, 22), fill=(*hp.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, hp.W, 30), fill=(*hp.PHOTON, 255))
    draw.text(title_xy, title, font=hp.font(hp.TIMES_BOLD, title_size), fill=hp.INK)
    if subtitle:
        draw.text((136, 190), subtitle, font=hp.font(hp.TIMES_ITALIC, 56), fill=hp.MUTED)
    draw.line((132, 286, hp.W - 132, 286), fill=(221, 226, 232, 255), width=3)
    hp.add_top_right_sphenix_logo_like_slide2(img)
    return img


def save_header_spec(path: Path) -> None:
    path.with_suffix(".header.json").write_text(
        json.dumps({"hp2026_main_header": HP2026_MAIN_HEADER}, indent=2) + "\n",
        encoding="utf-8",
    )


def card(base: Image.Image, box: tuple[int, int, int, int], accent: tuple[int, int, int]) -> ImageDraw.ImageDraw:
    hp.shadow(base, box, radius=12)
    draw = ImageDraw.Draw(base, "RGBA")
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*hp.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((box[0], box[1], box[0] + 12, box[3]), radius=6, fill=(*accent, 235))
    return draw


def draw_wrapped_lines(
    draw: ImageDraw.ImageDraw,
    lines: list[tuple[str, int, tuple[int, int, int], bool]],
    x: int,
    y: int,
    width: int,
    *,
    gap: int = 8,
) -> int:
    cursor = y
    for text, size, fill, bold in lines:
        fnt = hp.font(hp.TIMES_BOLD if bold else hp.TIMES, size)
        hp.draw_wrapped(draw, text, (x, cursor), width, fnt, fill=fill, line_gap=4)
        approx_lines = max(1, int((hp.text_box(draw, text, fnt)[0] + width - 1) / width))
        cursor += approx_lines * (size + 8) + gap
    return cursor


def draw_step_card(
    base: Image.Image,
    box: tuple[int, int, int, int],
    number: str,
    title: str,
    body: str,
    accent: tuple[int, int, int],
    fill: tuple[int, int, int],
) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    draw.rounded_rectangle(box, radius=12, fill=(*fill, 255), outline=(*accent, 210), width=2)
    badge = (box[0] + 24, box[1] + 22, box[0] + 78, box[1] + 76)
    draw.ellipse(badge, fill=(255, 255, 255, 255), outline=(*accent, 235), width=3)
    nf = hp.font(hp.TIMES_BOLD, 31)
    nw, nh = hp.text_box(draw, number, nf)
    draw.text((badge[0] + (54 - nw) / 2, badge[1] + (54 - nh) / 2 - 1), number, font=nf, fill=accent)
    draw.text((box[0] + 96, box[1] + 21), title, font=hp.font(hp.TIMES_BOLD, 31), fill=accent)
    hp.draw_wrapped(draw, body, (box[0] + 96, box[1] + 62), box[2] - box[0] - 124, hp.font(hp.TIMES, 26), fill=hp.INK, line_gap=3)


def annotate_bdt_score_regions(plot: Image.Image) -> Image.Image:
    """Mark non-tight/tight score regions without shading over the data."""
    out = plot.convert("RGBA").copy()

    # Coordinates measured on the cropped PPG12 Fig. 2 raster used by this deck.
    x0, y0, x1, y1 = 149, 9, 919, 748
    x_cut = int(round(x0 + 0.55 * (x1 - x0)))
    blue = hp.SPHENIX_BLUE
    red = (222, 64, 50)
    draw = ImageDraw.Draw(out, "RGBA")

    dash_y = y0 + 16
    while dash_y < y1 - 16:
        draw.line((x_cut, dash_y, x_cut, min(dash_y + 16, y1 - 16)), fill=(*hp.INK, 150), width=3)
        dash_y += 28

    def double_arrow(start: tuple[int, int], end: tuple[int, int], color: tuple[int, int, int]) -> None:
        draw.line((start[0], start[1], end[0], end[1]), fill=(*color, 240), width=6)
        head = 16
        for x, direction in ((start[0], 1), (end[0], -1)):
            y = start[1]
            draw.polygon(
                [
                    (x, y),
                    (x + direction * head, y - 11),
                    (x + direction * head, y + 11),
                ],
                fill=(*color, 245),
            )

    arrow_y = y0 + 382
    label_y = y0 + 326
    left_start = (x0 + 36, arrow_y)
    left_end = (x_cut - 28, arrow_y)
    right_start = (x_cut + 28, arrow_y)
    right_end = (x1 - 36, arrow_y)
    double_arrow(left_start, left_end, blue)
    double_arrow(right_start, right_end, red)

    label_font = hp.font(hp.TIMES_BOLD, 35)
    labels = [
        ("non-tight", blue, left_start[0], left_end[0]),
        ("tight", red, right_start[0], right_end[0]),
    ]
    for text, fill, lx0, lx1 in labels:
        w, h = hp.text_box(draw, text, label_font)
        tx = lx0 + ((lx1 - lx0) - w) / 2
        draw.rounded_rectangle(
            (tx - 16, label_y - 8, tx + w + 16, label_y + h + 8),
            radius=10,
            fill=(255, 255, 255, 220),
            outline=(*fill, 170),
            width=2,
        )
        draw.text((tx, label_y), text, font=label_font, fill=fill)
    return out


def draw_centered_rich_line(
    draw: ImageDraw.ImageDraw,
    segments: list[tuple[str, tuple[int, int, int], ImageFont.FreeTypeFont]],
    center_x: int,
    y: int,
) -> int:
    widths = [hp.text_box(draw, text, font)[0] for text, _, font in segments]
    heights = [hp.text_box(draw, text, font)[1] for text, _, font in segments]
    x = center_x - sum(widths) / 2
    for (text, fill, font), width in zip(segments, widths):
        draw.text((x, y), text, font=font, fill=fill)
        x += width
    return max(heights)


def draw_rich_line(
    draw: ImageDraw.ImageDraw,
    segments: list[tuple[str, tuple[int, int, int], ImageFont.FreeTypeFont]],
    x: int,
    y: int,
) -> int:
    heights = [hp.text_box(draw, text, font)[1] for text, _, font in segments]
    cursor = x
    for text, fill, font in segments:
        draw.text((cursor, y), text, font=font, fill=fill)
        cursor += hp.text_box(draw, text, font)[0]
    return max(heights)


def draw_main_bdt_slide() -> Image.Image:
    img = draw_shell(
        "Photon-ID BDT score defines the ID axis",
        None,
        title_size=108,
        title_xy=(132, 104),
    )
    draw = ImageDraw.Draw(img, "RGBA")

    main = (132, 318, 2390, 1288)
    d = card(img, main, hp.SPHENIX_BLUE)

    plot = hp.crop_visible(Image.open(hp.figure_path("fig2_bdt_score")).convert("RGBA"), white_threshold=252, pad=8)
    plot = annotate_bdt_score_regions(plot)
    hp.paste_fit(img, plot, (main[0] + 1010, main[1] + 26, main[2] - 58, main[3] - 32), anchor="center")

    left = (main[0] + 82, main[1] + 74, main[0] + 1020, main[3] - 78)
    section_gap = 34
    statement_h = 208
    statement_boxes = [
        (
            (left[0], left[1], left[2], left[1] + statement_h),
            (184, 119, 18),
            (255, 250, 237),
            "NCB & preselection",
            [
                [("NCB = ", hp.INK, hp.font(hp.TIMES_BOLD, 35)), ("non-collisional background removal", hp.INK, hp.font(hp.TIMES, 35))],
                [("Preselection = ", hp.INK, hp.font(hp.TIMES_BOLD, 35)), ("clean candidate sample before BDT scoring", hp.INK, hp.font(hp.TIMES, 35))],
            ],
        ),
        (
            (left[0], left[1] + statement_h + section_gap, left[2], left[1] + 2 * statement_h + section_gap),
            (84, 96, 111),
            (248, 249, 251),
            "Tight / non-tight regions",
            [
                [("non-tight = ", (35, 86, 220), hp.font(hp.TIMES_BOLD, 35)), ("BDT score left of the threshold", hp.INK, hp.font(hp.TIMES, 35))],
                [("tight = ", (222, 64, 50), hp.font(hp.TIMES_BOLD, 35)), ("BDT score right of the threshold", hp.INK, hp.font(hp.TIMES, 35))],
            ],
        ),
    ]
    for box, accent, fill, heading, rows in statement_boxes:
        draw.rounded_rectangle(box, radius=14, fill=(*fill, 255), outline=(*accent, 115), width=2)
        draw.text((box[0] + 34, box[1] + 25), heading, font=hp.font(hp.TIMES_BOLD, 41), fill=accent)
        y_rows = box[1] + 92
        for row in rows:
            y_rows += draw_rich_line(draw, row, box[0] + 42, y_rows) + 16

    take = (left[0], left[1] + 2 * statement_h + 2 * section_gap + 8, left[2], left[3])
    draw.line((take[0] + 18, take[1] + 6, take[2] - 18, take[1] + 6), fill=(214, 225, 236, 255), width=2)
    center_x = (take[0] + take[2]) // 2
    lead_font = hp.font(hp.TIMES_BOLD, 44)
    body_font = hp.font(hp.TIMES_BOLD, 44)
    line_gap = 36
    line1 = [
        ("Signal MC", (222, 64, 50), lead_font),
        (" peaks at high BDT score", hp.INK, body_font),
    ]
    line2 = [
        ("Inclusive MC", (35, 86, 220), lead_font),
        (" peaks at low BDT score", hp.INK, body_font),
    ]
    total_h = 2 * hp.text_box(draw, "Signal MC", lead_font)[1] + line_gap
    y = take[1] + ((take[3] - take[1]) - total_h) / 2 - 2
    y += draw_centered_rich_line(draw, line1, center_x, y) + line_gap
    draw_centered_rich_line(draw, line2, center_x, y)

    hp.draw_hp2026_identity_footer(img)
    return img


def draw_backup_shower_shape_slide() -> Image.Image:
    img = draw_shell(
        "Backup: shower-shape inputs after NCB cleaning",
        "BDTs use shower-shape inputs for NCB cleaning and tight/non-tight photon ID.",
    )
    draw = ImageDraw.Draw(img, "RGBA")

    plot_card = (132, 330, 1650, 1276)
    read_card = (1700, 330, 2390, 1276)
    d = card(img, plot_card, hp.PHOTON_DARK)
    d.text((plot_card[0] + 56, plot_card[1] + 30), "After NCB + preselection: shower-shape space", font=hp.font(hp.TIMES_BOLD, 44), fill=hp.INK)
    d.text((plot_card[0] + 56, plot_card[1] + 84), "data remain background-rich, but shape variables become interpretable ID inputs", font=hp.font(hp.TIMES_ITALIC, 32), fill=hp.MUTED)
    d.line((plot_card[0] + 56, plot_card[1] + 132, plot_card[2] - 36, plot_card[1] + 132), fill=(221, 228, 236), width=2)
    plot = hp.crop_visible(Image.open(hp.figure_path("fig1_shower_shape")).convert("RGBA"), white_threshold=252, pad=8)
    hp.paste_fit(img, plot, (plot_card[0] + 64, plot_card[1] + 154, plot_card[2] - 64, plot_card[3] - 54), anchor="center")

    d = card(img, read_card, hp.TEAL)
    d.text((read_card[0] + 56, read_card[1] + 34), "What the inputs show", font=hp.font(hp.TIMES_BOLD, 45), fill=hp.INK)
    d.line((read_card[0] + 56, read_card[1] + 96, read_card[2] - 36, read_card[1] + 96), fill=(221, 228, 236), width=2)
    items = [
        ("w_eta", "Shower width", "Prompt-like showers are narrower than inclusive-jet background."),
        ("E3x2 / E3x5", "Energy sharing", "Compact deposits sit closer to the prompt-photon pattern."),
        ("BDT role", "Main-talk use", "These inputs motivate the BDT score used as the photon-ID axis."),
    ]
    y = read_card[1] + 138
    for idx, (label, headline, body) in enumerate(items, start=1):
        box = (read_card[0] + 56, y, read_card[2] - 42, y + 184)
        fill = (239, 248, 253) if idx == 1 else (255, 249, 235) if idx == 2 else (239, 249, 247)
        accent = hp.SPHENIX_BLUE if idx == 1 else hp.PHOTON_DARK if idx == 2 else hp.TEAL
        draw.rounded_rectangle(box, radius=12, fill=(*fill, 255), outline=(*accent, 190), width=2)
        draw.rounded_rectangle((box[0] + 24, box[1] + 28, box[0] + 38, box[3] - 28), radius=7, fill=(*accent, 235))
        display_label = {
            "w_eta": "wη",
            "E3x2 / E3x5": "E3x2 / E3x5",
        }.get(label, label)
        draw.text((box[0] + 62, box[1] + 22), display_label, font=hp.font(hp.TIMES_BOLD, 35), fill=accent)
        draw.text((box[0] + 62, box[1] + 66), headline, font=hp.font(hp.TIMES_BOLD, 31), fill=hp.INK)
        hp.draw_wrapped(draw, body, (box[0] + 62, box[1] + 104), box[2] - box[0] - 92, hp.font(hp.TIMES, 29), fill=hp.INK, line_gap=4)
        y += 214

    hp.draw_hp2026_identity_footer(img)
    return img


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    main_img = draw_main_bdt_slide()
    backup_img = draw_backup_shower_shape_slide()
    main_img.convert("RGB").save(MAIN_PNG, "PNG")
    backup_img.convert("RGB").save(BACKUP_PNG, "PNG")
    save_header_spec(MAIN_PNG)
    save_header_spec(BACKUP_PNG)
    MANIFEST.write_text(
        json.dumps(
            {
                "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
                "deck_mutation": "none; local PNG prototypes only",
                "generator": str(Path(__file__).resolve()),
                "outputs": {
                    "main_bdt_focus": str(MAIN_PNG),
                    "backup_shower_shape": str(BACKUP_PNG),
                },
                "source_assets": {
                    "bdt_score": str(hp.figure_path("fig2_bdt_score")),
                    "shower_shape": str(hp.figure_path("fig1_shower_shape")),
                },
                "design_intent": "Move the BDT-score distribution to the main slide as the dominant ID-axis object, and move the two shower-shape input distributions to a backup slide.",
                "public_label_caveat": "Embedded paper figures currently carry the source plot labels from the available PPG12 draft assets.",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(MAIN_PNG)
    print(BACKUP_PNG)


if __name__ == "__main__":
    main()

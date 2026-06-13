#!/usr/bin/env python3
"""Standalone generator for HP2026 deck Slide 11: "Photon identification via BDT score".

Bifurcated from make_hp2026_slide14_bdt_focus_split.py (draw_main_bdt_slide) so this
slide can be iterated independently of Codex's concurrent work on that file. It imports
ONLY the stable base module make_hp2026_fulltalk_candidates (hp); the slide composition
and header geometry are copied/owned here.

Fixes vs the deck slide, to match the HP2026 main-header contract:
  - title size 108 -> 86, title position (132,104) -> (132,76);
  - title underline y=286 -> y=232;
  - main card uses the same full vertical band as the adjacent isolation/purity slides.
Inner content (plot, statement boxes, takeaway) is positioned relative to the card and
reflows automatically.

PNG-only candidate. Does not mutate Google Slides.
"""

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


OUTDIR = ROOT / "outputs/manual-20260611-hp2026-slide12-bdt-id-axis-fix"
MAIN_PNG = OUTDIR / "hp2026_slide12_bdt_id_axis_candidate.png"
SCRIPT_PATH = OUTDIR / "hp2026_slide12_bdt_id_axis_script.md"
MANIFEST = OUTDIR / "manifest.json"

# Contract header (matches slides 9 and 13).
HP2026_MAIN_HEADER = {
    "deck": "hp2026_main_talk",
    "title_font_size": 86,
    "subtitle_font_size": None,
    "title_xy": [132, 76],
    "subtitle_xy": None,
    "divider_y": 232,
}


def draw_shell(title: str, subtitle: str | None = None) -> Image.Image:
    img = Image.new("RGBA", (hp.W, hp.H), (*hp.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, hp.W, hp.H), fill=(*hp.SOFT_BG, 255))
    draw.rectangle((0, 0, hp.W, 22), fill=(*hp.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, hp.W, 30), fill=(*hp.PHOTON, 255))
    draw.text(tuple(HP2026_MAIN_HEADER["title_xy"]), title, font=hp.font(hp.TIMES_BOLD, HP2026_MAIN_HEADER["title_font_size"]), fill=hp.INK)
    draw.line((132, HP2026_MAIN_HEADER["divider_y"], hp.W - 132, HP2026_MAIN_HEADER["divider_y"]), fill=(221, 226, 232, 255), width=3)
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


def annotate_bdt_score_regions(plot: Image.Image) -> Image.Image:
    """Mark background-like/photon-like score regions without shading over the data."""
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
        ("background-like", blue, left_start[0], left_end[0]),
        ("photon-like", red, right_start[0], right_end[0]),
    ]
    for text, fill, lx0, lx1 in labels:
        local_font = label_font if text == "photon-like" else hp.font(hp.TIMES_BOLD, 30)
        w, h = hp.text_box(draw, text, local_font)
        tx = lx0 + ((lx1 - lx0) - w) / 2
        draw.rounded_rectangle(
            (tx - 16, label_y - 8, tx + w + 16, label_y + h + 8),
            radius=10,
            fill=(255, 255, 255, 220),
            outline=(*fill, 170),
            width=2,
        )
        draw.text((tx, label_y), text, font=local_font, fill=fill)
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
    # Contract header: title 86 @ (132,76), underline @ 232.
    img = draw_shell("Photon identification via BDT score", None)
    draw = ImageDraw.Draw(img, "RGBA")

    # Contract card height: same full vertical band used by the neighboring
    # isolation and purity slides.
    main = (132, 258, 2390, 1308)
    d = card(img, main, hp.SPHENIX_BLUE)

    plot = hp.crop_visible(Image.open(hp.figure_path("fig2_bdt_score")).convert("RGBA"), white_threshold=252, pad=8)
    plot = annotate_bdt_score_regions(plot)
    hp.paste_fit(img, plot, (main[0] + 1010, main[1] + 26, main[2] - 58, main[3] - 32), anchor="center")

    left = (main[0] + 82, main[1] + 74, main[0] + 1020, main[3] - 78)
    section_gap = 38
    statement_h = 238
    row_size = 40
    heading_size = 47
    statement_boxes = [
        (
            (left[0], left[1], left[2], left[1] + statement_h),
            (184, 119, 18),
            (255, 250, 237),
            "NCB & preselection",
            [
                [("NCB = ", hp.INK, hp.font(hp.TIMES_BOLD, row_size)), ("non-collisional background removal", hp.INK, hp.font(hp.TIMES, row_size))],
                [("Preselection = ", hp.INK, hp.font(hp.TIMES_BOLD, row_size)), ("clean sample before BDT scoring", hp.INK, hp.font(hp.TIMES, row_size))],
            ],
        ),
        (
            (left[0], left[1] + statement_h + section_gap, left[2], left[1] + 2 * statement_h + section_gap),
            (84, 96, 111),
            (248, 249, 251),
            "BDT score regions",
            [
                [("background-like = ", (35, 86, 220), hp.font(hp.TIMES_BOLD, row_size)), ("BDT score left of the threshold", hp.INK, hp.font(hp.TIMES, row_size))],
                [("photon-like = ", (222, 64, 50), hp.font(hp.TIMES_BOLD, row_size)), ("BDT score right of the threshold", hp.INK, hp.font(hp.TIMES, row_size))],
            ],
        ),
    ]
    for box, accent, fill, heading, rows in statement_boxes:
        draw.rounded_rectangle(box, radius=14, fill=(*fill, 255), outline=(*accent, 115), width=2)
        draw.text((box[0] + 34, box[1] + 28), heading, font=hp.font(hp.TIMES_BOLD, heading_size), fill=accent)
        y_rows = box[1] + 108
        for row in rows:
            y_rows += draw_rich_line(draw, row, box[0] + 42, y_rows) + 18

    take = (left[0], left[1] + 2 * statement_h + 2 * section_gap + 8, left[2], left[3])
    draw.line((take[0] + 18, take[1] + 6, take[2] - 18, take[1] + 6), fill=(214, 225, 236, 255), width=2)
    center_x = (take[0] + take[2]) // 2
    lead_font = hp.font(hp.TIMES_BOLD, 50)
    body_font = hp.font(hp.TIMES_BOLD, 50)
    line_gap = 42
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


def write_script() -> None:
    SCRIPT_PATH.write_text(
        """# HP2026 Slide 11 Script - Photon identification via BDT score

The boosted decision tree turns the shower-shape handles into a single photon-ID score, and this slide shows that score. On the left I remind the audience what the candidates have already been through: NCB removes non-collisional background, and preselection gives a clean candidate sample before any BDT scoring.

The plot on the right is the normalized photon-ID BDT score. Signal MC peaks at high score and inclusive MC peaks at low score, and the data sit in between, as expected for a background-rich candidate sample. I split the axis at the threshold: candidates to the right are photon-like, and candidates to the left are background-like. That score split is the identification axis I use, together with isolation, to build the data-driven purity measurement on the next slides.
""",
        encoding="utf-8",
    )


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    img = draw_main_bdt_slide()
    img.convert("RGB").save(MAIN_PNG, "PNG")
    save_header_spec(MAIN_PNG)
    write_script()
    MANIFEST.write_text(
        json.dumps(
            {
                "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
                "google_slides_mutation": False,
                "deck_slide": "HPslides_v1 Slide 11 ('Photon identification via BDT score')",
                "bifurcated_from": "scripts/slides/hp2026/fulltalk/make_hp2026_slide14_bdt_focus_split.py::draw_main_bdt_slide (Codex-owned; not edited)",
                "imports": ["make_hp2026_fulltalk_candidates (base helpers only)"],
                "figure_source": "PPG12 paper Fig. 2 (photon-ID BDT score) -> fig2_bdt_score, annotated with background-like/photon-like regions",
                "header_contract": HP2026_MAIN_HEADER,
                "fixes": [
                    "Title 108 -> 86 and (132,104) -> (132,76), matching slides 9 and 13.",
                    "Title underline y=286 -> y=232.",
                    "Main card uses (132,258,2390,1308), matching the taller neighboring isolation/purity slide band.",
                    "Left-side statement and takeaway typography increased by roughly 10-15%, with taller statement boxes and cleaner spacing.",
                ],
                "labels_note": "Plot still reads 'sPHENIX Internal'; public HP version needs the released/Preliminary-labeled BDT-score figure.",
                "output_png": str(MAIN_PNG.relative_to(ROOT)),
                "companion_script": str(SCRIPT_PATH.relative_to(ROOT)),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(MAIN_PNG)
    print(SCRIPT_PATH)
    print(MANIFEST)


if __name__ == "__main__":
    main()

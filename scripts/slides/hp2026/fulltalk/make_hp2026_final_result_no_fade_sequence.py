#!/usr/bin/env python3
"""Render a no-fade HP2026 final-result overview sequence.

This is a local PNG prototype only. It reuses the PPG12 paper Fig. 8 crop from
the full-talk generator. The plot stays intact in the same style as the current
Slide 18/19 sequence; no blur/fade treatment is used.
"""

from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter

import make_hp2026_fulltalk_candidates as full


ROOT = full.ROOT
OUTDIR = ROOT / "outputs/manual-20260611-slides15-17-header-normalized/final_result"
SCRIPT_DIR = OUTDIR / "speaker_scripts"
W, H = full.W, full.H
HP2026_MAIN_HEADER = {
    "deck": "hp2026_main_talk",
    "title_font_size": 86,
    "subtitle_font_size": None,
    "title_xy": [132, 76],
    "subtitle_xy": None,
    "divider_y": 232,
}

PLOT_PANEL_BOX = (132, 270, 1048, 1290)
RIGHT_CARD_X0 = 1090
RIGHT_CARD_X1 = 2390
FULL_FIG8_VISIBLE_SOURCE_Y0 = 49
TOP_MIDDLE_SOURCE_Y1 = 1378
FULL_FIG8_VISIBLE_SOURCE_HEIGHT = 1827


DATA_BLUE = (67, 132, 217)
PYTHIA_ORANGE = (238, 126, 40)
JETPHOX_MAGENTA = (214, 86, 169)
VOGELSANG_GREEN = (91, 174, 96)
PDF_OLIVE = (117, 143, 89)
PDF_BLUE = (82, 88, 210)
PDF_SLATE = (70, 83, 122)
PDF_STEEL = (58, 122, 169)
PDF_INDIGO = (82, 96, 174)
READOUT_GOLD = (195, 138, 33)
SOFT_BLUE = (235, 244, 252)
SOFT_YELLOW = (255, 249, 232)
SOFT_TEAL = (237, 248, 250)
SOFT_GRAY = (247, 249, 252)
SOFT_SLATE = (246, 248, 252)
SOFT_INDIGO = (244, 246, 255)


def fnt(path: Path, size: int):
    return full.font(path, size)


def text_size(draw: ImageDraw.ImageDraw, text: str, font) -> tuple[int, int]:
    return full.text_box(draw, text, font)


def draw_header(draw: ImageDraw.ImageDraw, title: str, subtitle: str | None = None) -> None:
    draw.rectangle((0, 0, W, 22), fill=full.SPHENIX_BLUE)
    draw.rectangle((0, 22, W, 30), fill=full.PHOTON)
    draw.text(tuple(HP2026_MAIN_HEADER["title_xy"]), title, font=fnt(full.TIMES_BOLD, HP2026_MAIN_HEADER["title_font_size"]), fill=full.INK)
    if subtitle and HP2026_MAIN_HEADER["subtitle_xy"] and HP2026_MAIN_HEADER["subtitle_font_size"]:
        draw.text(
            tuple(HP2026_MAIN_HEADER["subtitle_xy"]),
            subtitle,
            font=fnt(full.TIMES_ITALIC, HP2026_MAIN_HEADER["subtitle_font_size"]),
            fill=full.BLUE,
        )
    y = HP2026_MAIN_HEADER["divider_y"]
    draw.line((132, y, W - 132, y), fill=(221, 226, 232), width=3)


def base_slide(title: str, subtitle: str | None = None) -> Image.Image:
    img = Image.new("RGBA", (W, H), (*full.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw_header(draw, title, subtitle)
    full.add_top_right_sphenix_logo_like_slide2(img)
    return img


def soft_shadow(base: Image.Image, box: tuple[int, int, int, int], radius: int = 10) -> None:
    layer = Image.new("RGBA", base.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(layer, "RGBA")
    d.rounded_rectangle((box[0] + 7, box[1] + 9, box[2] + 7, box[3] + 9), radius=radius, fill=(30, 42, 58, 26))
    layer = layer.filter(ImageFilter.GaussianBlur(14))
    base.alpha_composite(layer)


def panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], fill=(255, 255, 255), outline=full.PANEL_EDGE) -> None:
    draw.rounded_rectangle(box, radius=10, fill=(*fill, 255), outline=(*outline, 255), width=2)


def paste_image_panel(
    base: Image.Image,
    img: Image.Image,
    box: tuple[int, int, int, int],
    *,
    pad: int = 22,
    fill=(255, 255, 255),
) -> tuple[int, int, int, int]:
    draw = ImageDraw.Draw(base, "RGBA")
    soft_shadow(base, box)
    panel(draw, box, fill=fill)
    fitted = full.fit(img.convert("RGBA"), box[2] - box[0] - 2 * pad, box[3] - box[1] - 2 * pad)
    x = box[0] + ((box[2] - box[0]) - fitted.width) // 2
    y = box[1] + ((box[3] - box[1]) - fitted.height) // 2
    base.alpha_composite(fitted, (x, y))
    return (x, y, x + fitted.width, y + fitted.height)


def wrapped_lines(draw: ImageDraw.ImageDraw, text: str, font, max_width: int) -> list[str]:
    words = text.split()
    lines: list[str] = []
    cur = ""
    for word in words:
        trial = word if not cur else f"{cur} {word}"
        if text_size(draw, trial, font)[0] <= max_width:
            cur = trial
        else:
            if cur:
                lines.append(cur)
            cur = word
    if cur:
        lines.append(cur)
    return lines


def draw_text_block(
    draw: ImageDraw.ImageDraw,
    text: str,
    box: tuple[int, int, int, int],
    font,
    *,
    fill=full.MUTED,
    line_gap: int = 8,
    valign: str = "center",
) -> int:
    lines = wrapped_lines(draw, text, font, box[2] - box[0])
    heights = [text_size(draw, line, font)[1] for line in lines]
    total = sum(heights) + line_gap * max(0, len(lines) - 1)
    if valign == "top":
        y = box[1]
    else:
        y = box[1] + ((box[3] - box[1]) - total) // 2
    for line, height in zip(lines, heights):
        draw.text((box[0], y), line, font=font, fill=fill)
        y += height + line_gap
    return y


def draw_label_value(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    label: str,
    body: str,
    accent: tuple[int, int, int],
    max_width: int,
    *,
    body_size: int = 31,
) -> int:
    draw.rounded_rectangle((x, y + 10, x + 20, y + 42), radius=6, fill=(*accent, 255))
    label_font = fnt(full.TIMES_BOLD, body_size)
    body_font = fnt(full.TIMES, body_size)
    lw, _ = text_size(draw, label, label_font)
    draw.text((x + 36, y), label, font=label_font, fill=accent)
    return full.draw_wrapped(
        draw,
        body,
        (x + 48 + lw, y + 1),
        max_width - lw - 48,
        body_font,
        fill=full.MUTED,
        line_gap=7,
    )


def callout_card(
    base: Image.Image,
    box: tuple[int, int, int, int],
    title: str,
    body: str,
    accent: tuple[int, int, int],
    *,
    fill=(255, 255, 255),
    body_size: int = 34,
) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    soft_shadow(base, box)
    panel(draw, box, fill=fill)
    draw.rounded_rectangle((box[0], box[1], box[0] + 14, box[3]), radius=7, fill=(*accent, 255))
    draw.text((box[0] + 42, box[1] + 26), title, font=fnt(full.TIMES_BOLD, 42), fill=full.INK)
    draw_text_block(
        draw,
        body,
        (box[0] + 42, box[1] + 86, box[2] - 36, box[3] - 24),
        fnt(full.TIMES, body_size),
        fill=full.MUTED,
        line_gap=9,
        valign="center",
    )


def wrap_text_for_width(draw: ImageDraw.ImageDraw, text: str, font, width: int) -> list[str]:
    return wrapped_lines(draw, text, font, width)


def callout_row_card(
    base: Image.Image,
    box: tuple[int, int, int, int],
    title: str,
    rows: list[tuple[str, str]],
    accent: tuple[int, int, int],
    *,
    fill=(255, 255, 255),
    title_size: int = 45,
    label_size: int = 36,
    body_size: int = 36,
) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    soft_shadow(base, box)
    panel(draw, box, fill=fill)
    draw.rounded_rectangle((box[0], box[1], box[0] + 14, box[3]), radius=7, fill=(*accent, 255))
    draw.text((box[0] + 42, box[1] + 24), title, font=fnt(full.TIMES_BOLD, title_size), fill=full.INK)

    body_top = box[1] + 92
    body_bottom = box[3] - 24
    row_h = (body_bottom - body_top) / max(1, len(rows))
    label_font = fnt(full.TIMES_BOLD, label_size)
    body_font = fnt(full.TIMES, body_size)
    label_x = box[0] + 44
    body_x = box[0] + 252
    body_w = box[2] - body_x - 38
    for i, (label, body) in enumerate(rows):
        cell_top = int(body_top + i * row_h)
        cell_bottom = int(body_top + (i + 1) * row_h)
        lines = wrap_text_for_width(draw, body, body_font, body_w)
        line_heights = [text_size(draw, line, body_font)[1] for line in lines]
        body_total = sum(line_heights) + 6 * max(0, len(lines) - 1)
        label_h = text_size(draw, label, label_font)[1]
        row_total = max(label_h, body_total)
        y = cell_top + ((cell_bottom - cell_top) - row_total) // 2
        draw.text((label_x, y), label, font=label_font, fill=accent)
        yy = y
        for line, height in zip(lines, line_heights):
            draw.text((body_x, yy), line, font=body_font, fill=full.INK)
            yy += height + 6


def physics_message_card(
    base: Image.Image,
    box: tuple[int, int, int, int],
    bullets: list[str],
    accent: tuple[int, int, int],
    *,
    title_size: int = 45,
    body_size: int = 38,
) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    soft_shadow(base, box)
    panel(draw, box, fill=(255, 255, 255))
    draw.rounded_rectangle((box[0], box[1], box[0] + 14, box[3]), radius=7, fill=(*accent, 255))
    draw.text((box[0] + 42, box[1] + 24), "Physics message", font=fnt(full.TIMES_BOLD, title_size), fill=full.INK)

    bullet_font = fnt(full.TIMES_BOLD, body_size + 8)
    body_font = fnt(full.TIMES, body_size)
    body_x = box[0] + 86
    body_w = box[2] - body_x - 42
    row_gap = 28

    wrapped = []
    row_heights = []
    for bullet in bullets:
        lines = wrap_text_for_width(draw, bullet, body_font, body_w)
        heights = [text_size(draw, line, body_font)[1] for line in lines]
        wrapped.append((lines, heights))
        row_heights.append(sum(heights) + 6 * max(0, len(lines) - 1))
    body_top = box[1] + 98
    body_bottom = box[3] - 34
    total_h = sum(row_heights) + row_gap * max(0, len(row_heights) - 1)
    y = body_top + max(0, (body_bottom - body_top - total_h) // 2)

    for (lines, heights), row_h in zip(wrapped, row_heights):
        draw.text((box[0] + 48, y - 1), "\u2022", font=bullet_font, fill=accent)
        yy = y
        for line, lh in zip(lines, heights):
            draw.text((body_x, yy), line, font=body_font, fill=full.INK)
            yy += lh + 6
        y = yy + row_gap


def physics_message_row_card(
    base: Image.Image,
    box: tuple[int, int, int, int],
    rows: list[tuple[str, str]],
    accent: tuple[int, int, int],
) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    soft_shadow(base, box)
    panel(draw, box, fill=(255, 255, 255))
    draw.rounded_rectangle((box[0], box[1], box[0] + 14, box[3]), radius=7, fill=(*accent, 255))
    draw.text((box[0] + 42, box[1] + 30), "Physics message", font=fnt(full.TIMES_BOLD, 54), fill=full.INK)

    top = box[1] + 118
    bottom = box[3] - 42
    row_h = (bottom - top) / max(1, len(rows))
    label_font = fnt(full.TIMES_BOLD, 45)
    body_font = fnt(full.TIMES, 45)
    label_x = box[0] + 64
    body_x = box[0] + 255
    body_w = box[2] - body_x - 52

    for i, (label, body) in enumerate(rows):
        y0 = int(top + i * row_h)
        y1 = int(top + (i + 1) * row_h)
        if i:
            draw.line((box[0] + 42, y0 - 10, box[2] - 42, y0 - 10), fill=(224, 230, 236, 255), width=2)

        lines = wrap_text_for_width(draw, body, body_font, body_w)
        line_heights = [text_size(draw, line, body_font)[1] for line in lines]
        body_total = sum(line_heights) + 8 * max(0, len(lines) - 1)
        label_h = text_size(draw, label, label_font)[1]
        row_total = max(label_h, body_total)
        y = y0 + ((y1 - y0) - row_total) // 2

        draw.text((label_x, y), label, font=label_font, fill=accent)
        yy = y
        for line, lh in zip(lines, line_heights):
            draw.text((body_x, yy), line, font=body_font, fill=full.INK)
            yy += lh + 8


def result_crops() -> dict[str, Image.Image]:
    src = Image.open(full.figure_path("fig8_cross_section")).convert("RGBA")
    # Crop coordinates are in the original PPG12 Fig. 8 crop. They preserve the
    # paper's axes, labels, legend, and plotted objects; slide code only frames.
    crops = {
        "top": src.crop((0, 0, src.width, 930)),
        "middle": src.crop((0, 950, src.width, 1378)),
        "bottom": src.crop((0, 1370, src.width, src.height)),
        "full": full.crop_visible(src, white_threshold=252, pad=18),
    }
    return {
        key: Image.alpha_composite(
            Image.new("RGBA", tightened.size, (255, 255, 255, 255)),
            tightened,
        )
        for key, im in crops.items()
        for tightened in [full.crop_visible(im, white_threshold=252, pad=12)]
    }


def full_result_plot() -> Image.Image:
    src = Image.open(full.figure_path("fig8_cross_section")).convert("RGBA")
    plot = full.crop_visible(src, white_threshold=252, pad=18)
    return Image.alpha_composite(Image.new("RGBA", plot.size, (255, 255, 255, 255)), plot)


def top_middle_result_plot() -> Image.Image:
    src = Image.open(full.figure_path("fig8_cross_section")).convert("RGBA")
    # Stop at the paper figure's boundary after the theory/data panel. This is
    # a true crop, not a white overlay, so no lower-panel labels can leak into
    # the first build state.
    plot = src.crop((0, 0, src.width, 1378))
    plot = full.crop_visible(plot, white_threshold=252, pad=12)
    return Image.alpha_composite(Image.new("RGBA", plot.size, (255, 255, 255, 255)), plot)


def place_full_result_plot(
    base: Image.Image,
    *,
    mask_bottom: bool = False,
) -> tuple[tuple[int, int, int, int], tuple[int, int, int, int]]:
    plot = full_result_plot()
    panel_box = PLOT_PANEL_BOX
    placed = paste_image_panel(base, plot, panel_box, pad=22)
    if mask_bottom:
        draw = ImageDraw.Draw(base, "RGBA")
        x0, y0, x1, y1 = placed
        # This is a hard reveal mask, not a blur/fade. It keeps the full-figure
        # geometry fixed while withholding the PDF-sensitivity panel until the
        # next click.
        boundary_frac = (TOP_MIDDLE_SOURCE_Y1 - FULL_FIG8_VISIBLE_SOURCE_Y0) / FULL_FIG8_VISIBLE_SOURCE_HEIGHT
        mask_y0 = y0 + int(round((y1 - y0) * boundary_frac))
        draw.rectangle((x0, mask_y0, x1, y1), fill=(255, 255, 255, 255))
    return placed, panel_box


def place_top_middle_result_plot(base: Image.Image) -> tuple[tuple[int, int, int, int], tuple[int, int, int, int]]:
    plot = top_middle_result_plot()
    placed = paste_image_panel(base, plot, PLOT_PANEL_BOX, pad=16)
    return placed, PLOT_PANEL_BOX


def draw_plot_caption(base: Image.Image, box: tuple[int, int, int, int], text: str, accent: tuple[int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    draw.rounded_rectangle(box, radius=8, fill=(247, 250, 253, 255), outline=(218, 226, 235, 255), width=2)
    draw.rounded_rectangle((box[0], box[1], box[0] + 12, box[3]), radius=6, fill=(*accent, 255))
    draw_text_block(draw, text, (box[0] + 34, box[1] + 8, box[2] - 28, box[3] - 8), fnt(full.TIMES_BOLD, 35), fill=full.BLUE)


def draw_theory_data_marker_callouts(
    base: Image.Image,
    placed: tuple[int, int, int, int],
    panel_box: tuple[int, int, int, int],
) -> None:
    """Annotate Fig. 8 theory/data markers in the lower white reveal space."""
    draw = ImageDraw.Draw(base, "RGBA")
    px0, py0, px1, py1 = placed
    pw = px1 - px0
    ph = py1 - py0

    callouts = [
        {
            "box": (panel_box[0] + 40, panel_box[3] - 250, panel_box[0] + 400, panel_box[3] - 130),
            "title": "blue band",
            "body": ["experimental", "systematic"],
            "color": DATA_BLUE,
            "target": (px0 + int(0.152 * pw), py0 + int(0.684 * ph)),
            "anchor": "center",
        },
        {
            "box": (panel_box[0] + 335, panel_box[3] - 98, panel_box[0] + 660, panel_box[3] - 8),
            "title": "pink boxes",
            "body": ["JETPHOX scale variation"],
            "color": JETPHOX_MAGENTA,
            "target": (px0 + int(0.456 * pw), py0 + int(0.648 * ph)),
        },
        {
            "box": (panel_box[0] + 660, panel_box[3] - 202, panel_box[2] - 22, panel_box[3] - 82),
            "title": "green boxes",
            "body": ["Vogelsang scale", "variation"],
            "color": VOGELSANG_GREEN,
            "target": (px0 + int(0.742 * pw), py0 + int(0.512 * ph)),
        },
    ]

    for item in callouts:
        box = item["box"]
        color = item["color"]
        cx = (box[0] + box[2]) // 2
        cy = (box[1] + box[3]) // 2
        anchor = (cx, cy) if item.get("anchor") == "center" else (cx, box[1] - 8)
        full.draw_arrow(draw, anchor, item["target"], fill=color, width=3)

    for item in callouts:
        box = item["box"]
        color = item["color"]
        draw.rounded_rectangle(box, radius=9, fill=(255, 255, 255, 238), outline=(*color, 210), width=2)
        draw.rounded_rectangle((box[0], box[1], box[0] + 10, box[3]), radius=5, fill=(*color, 230))
        title_font = fnt(full.TIMES_BOLD, 34)
        body_font = fnt(full.TIMES, 28)
        draw.text((box[0] + 24, box[1] + 9), item["title"], font=title_font, fill=color)
        body_y = box[1] + (48 if len(item["body"]) == 1 else 54)
        draw.multiline_text(
            (box[0] + 24, body_y),
            "\n".join(item["body"]),
            font=body_font,
            fill=full.MUTED,
            spacing=2,
        )


def slide18_top_middle_overview() -> tuple[Path, Path]:
    img = base_slide(
        "Main result: isolated prompt-photon cross section",
    )
    placed, panel_box = place_full_result_plot(img, mask_bottom=True)
    draw_theory_data_marker_callouts(img, placed, panel_box)
    callout_row_card(
        img,
        (RIGHT_CARD_X0, 276, RIGHT_CARD_X1, 570),
        "Corrected spectrum",
        [
            ("Data", "corrected isolated-photon cross section"),
            ("Scale", "steep falloff with photon energy"),
            ("Corrections", "purity, efficiency, unfolding, energy scale"),
        ],
        DATA_BLUE,
        fill=(255, 255, 255),
        body_size=38,
    )
    callout_row_card(
        img,
        (RIGHT_CARD_X0, 612, RIGHT_CARD_X1, 928),
        "Theory comparison",
        [
            ("Unity", "dashed line marks theory / data = 1"),
            ("Band", "blue band is experimental systematic uncertainty"),
            ("Result", "NLO pQCD agrees within uncertainties"),
        ],
        JETPHOX_MAGENTA,
        fill=(255, 255, 255),
        body_size=38,
    )
    physics_message_card(
        img,
        (RIGHT_CARD_X0, 970, RIGHT_CARD_X1, 1248),
        [
            "Measured cross section agrees with NLO pQCD within uncertainties",
            "Provides the p+p baseline for future heavy-ion photon and \u03b3+jet studies",
        ],
        full.TEAL,
    )
    full.draw_hp2026_identity_footer(img)
    png = OUTDIR / "slide18_final_result_top_middle_overview_no_fade.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        "slide18_final_result_top_middle_overview_no_fade_script.md",
        "I would start the result by reading the top two panels. The top panel is the corrected isolated prompt-photon cross section in 200 GeV p+p collisions. The spectrum falls steeply, which is why the correction chain we just walked through matters. Then the middle panel turns the spectrum into the theory comparison. The dashed line is perfect theory/data agreement, and the blue band is the experimental systematic uncertainty. The main result statement is that the NLO pQCD calculations are compatible with the corrected sPHENIX measurement within uncertainties.",
    )
    return png, script


def slide19_full_pdf_overview() -> tuple[Path, Path]:
    img = base_slide(
        "Main result: PDF sensitivity",
    )
    place_full_result_plot(img, mask_bottom=False)
    callout_row_card(
        img,
        (RIGHT_CARD_X0, 276, RIGHT_CARD_X1, 552),
        "PDF sensitivity",
        [
            ("Method", "repeat JETPHOX with different proton PDFs"),
            ("Readout", "PDF choices produce visible theory/data shape differences"),
        ],
        DATA_BLUE,
        fill=(255, 255, 255),
        body_size=42,
    )
    physics_message_row_card(
        img,
        (RIGHT_CARD_X0, 594, RIGHT_CARD_X1, 1248),
        [
            ("Trend", "proton-PDF choices change the theory/data shape, not just the normalization"),
            ("Spread", "visible, but smaller than the NLO scale variation over most of the measured range"),
        ],
        READOUT_GOLD,
    )
    full.draw_hp2026_identity_footer(img)
    png = OUTDIR / "slide19_final_result_pdf_reveal_no_fade.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        "slide19_final_result_pdf_reveal_no_fade_script.md",
        "On the next click, I would reveal the lower panel and read the figure as a complete result. The bottom panel repeats the JETPHOX calculation with different proton PDF choices. The visible spread means the isolated-photon measurement is not only a rate measurement; it carries interpretable dependence on the pQCD inputs. Together with the top and middle panels, this says that the first sPHENIX isolated prompt-photon cross section establishes the p+p reference for the photon and gamma-plus-jet program.",
    )
    return png, script


def save_script(name: str, body: str) -> Path:
    SCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    path = SCRIPT_DIR / name
    path.write_text(body.strip() + "\n", encoding="utf-8")
    return path


def save_header_spec(png: Path) -> None:
    png.with_suffix(".header.json").write_text(
        json.dumps({"hp2026_main_header": HP2026_MAIN_HEADER}, indent=2) + "\n",
        encoding="utf-8",
    )


def make_contact_sheet(paths: list[Path]) -> Path:
    thumbs = []
    for path in paths:
        im = Image.open(path).convert("RGB")
        im.thumbnail((720, 405), Image.Resampling.LANCZOS)
        thumbs.append((path.name, im.copy()))
    sheet = Image.new("RGB", (760, len(paths) * 465 + 30), (245, 247, 250))
    d = ImageDraw.Draw(sheet)
    y = 20
    for name, im in thumbs:
        sheet.paste(im, (20, y))
        d.text((20, y + 412), name, font=fnt(full.TIMES, 22), fill=(40, 48, 58))
        y += 465
    out = OUTDIR / "slide18_19_final_result_no_fade_contact_sheet.png"
    sheet.save(out, "PNG")
    return out


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    SCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    full.prepare_figures()
    outputs = []
    scripts = []
    for func in (slide18_top_middle_overview, slide19_full_pdf_overview):
        png, script = func()
        save_header_spec(png)
        outputs.append(png)
        scripts.append(script)
    contact = make_contact_sheet(outputs)
    manifest = {
        "created": full.datetime.now().isoformat(timespec="seconds"),
        "source_generator": str(Path(__file__).resolve()),
        "source_figure": str(full.figure_path("fig8_cross_section")),
        "google_slides_mutated": False,
        "design_intent": "No fading; Slides 18 and 19 use identical full Fig. 8 placement. Slide 18 hides only the lower PDF panel with a hard white cover, so the transition reads as adding the bottom row.",
        "outputs": [str(p) for p in outputs],
        "speaker_scripts": [str(p) for p in scripts],
        "contact_sheet": str(contact),
        "hp2026_main_header": HP2026_MAIN_HEADER,
    }
    (OUTDIR / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

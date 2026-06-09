#!/usr/bin/env python3
"""Local-only redesign candidates for HP2026 Slide 16 yield-correction flow."""

from __future__ import annotations

import json
import math
import shutil
import sys
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure


ROOT = next(
    p
    for p in Path(__file__).resolve().parents
    if (p / "README.md").exists() and (p / "scripts").exists() and (p / "src").exists()
)
SCRIPT_DIR = ROOT / "scripts/slides/hp2026/fulltalk"
sys.path.insert(0, str(SCRIPT_DIR))

import make_hp2026_closing_three_candidates as c3  # noqa: E402
import make_hp2026_fulltalk_candidates as ft  # noqa: E402


ARIAL_UNICODE = Path("/System/Library/Fonts/Supplemental/Arial Unicode.ttf")
OUTDIR = ROOT / "outputs/manual-20260608-yield_flow_redesign_refined"
CANDIDATE_DIR = OUTDIR / "internal_candidates"
FINAL_PNG = OUTDIR / "slide16_yield_flow_redesign_refined_purity_card.png"
MANIFEST = OUTDIR / "slide16_yield_flow_redesign_refined_purity_card_manifest.json"


def trim_white_margins(img: Image.Image, tolerance: int = 12, pad: int = 8) -> Image.Image:
    rgb = img.convert("RGB")
    white = Image.new("RGB", rgb.size, "white")
    diff = ImageChops.difference(rgb, white).convert("L")
    bbox = diff.point(lambda value: 255 if value > tolerance else 0).getbbox()
    if bbox is None:
        return img
    x0, y0, x1, y1 = bbox
    return img.crop((max(0, x0 - pad), max(0, y0 - pad), min(img.width, x1 + pad), min(img.height, y1 + pad)))


def base_slide() -> Image.Image:
    img = ft.base_slide(
        "From selected candidates to a corrected photon yield",
        "Sidebands determine purity; purity, efficiency, and unfolding convert Region A candidates to particle level.",
    )
    ft.add_top_right_sphenix_logo_like_slide2(img)
    return img


def draw_card(base: Image.Image, box: tuple[int, int, int, int], accent: tuple[int, int, int], title: str, subtitle: str) -> ImageDraw.ImageDraw:
    draw = ImageDraw.Draw(base, "RGBA")
    ft.shadow(base, box, radius=12)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((box[0] + 18, box[1] + 24, box[0] + 30, box[3] - 24), radius=6, fill=(*accent, 255))
    draw.text((box[0] + 54, box[1] + 24), title, font=ft.font(ft.TIMES_BOLD, 33), fill=ft.INK)
    ft.draw_wrapped(
        draw,
        subtitle,
        (box[0] + 54, box[1] + 66),
        box[2] - box[0] - 90,
        ft.font(ft.TIMES_ITALIC, 22),
        fill=ft.MUTED,
        line_gap=3,
    )
    draw.line((box[0] + 54, box[1] + 110, box[2] - 34, box[1] + 110), fill=(221, 228, 236), width=2)
    return draw


def draw_plain_accent_card(base: Image.Image, box: tuple[int, int, int, int], accent: tuple[int, int, int]) -> ImageDraw.ImageDraw:
    draw = ImageDraw.Draw(base, "RGBA")
    ft.shadow(base, box, radius=12)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((box[0] + 18, box[1] + 24, box[0] + 30, box[3] - 24), radius=6, fill=(*accent, 255))
    return draw


def draw_leakage_label(base: Image.Image, center: tuple[int, int], text: str, color: tuple[int, int, int], angle: float = 0) -> None:
    font = ft.font(ft.TIMES_BOLD, 22)
    tmp = Image.new("RGBA", (180, 46), (255, 255, 255, 0))
    d = ImageDraw.Draw(tmp, "RGBA")
    tw, th = ft.text_box(d, text, font)
    d.rounded_rectangle((6, 6, tw + 24, 38), radius=8, fill=(255, 255, 255, 245), outline=(*color, 185), width=2)
    d.text((16, 10), text, font=font, fill=color)
    tmp = tmp.crop((0, 0, tw + 34, 46)).rotate(-angle, expand=True, resample=Image.Resampling.BICUBIC)
    base.alpha_composite(tmp, (int(center[0] - tmp.width / 2), int(center[1] - tmp.height / 2)))


def draw_leakage_map(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box

    header_font = ft.font(ft.TIMES_BOLD, 35)
    definition_font = ft.font(ft.TIMES, 27)
    header_x, header_y = x0 + 8, y0 + 4
    draw.text((header_x, header_y), "Leakage", font=header_font, fill=ft.BLUE)
    header_w = ft.text_box(draw, "Leakage", header_font)[0]
    draw.text(
        (header_x + header_w + 9, header_y + 7),
        "= truth-matched signal photons in PYTHIA MC found in regions B/C/D",
        font=definition_font,
        fill=ft.INK,
    )

    map_x0, map_y0 = x0 + 8, y0 + 82
    available_w = x1 - x0 - 16
    cell_w = min(350, max(310, int((available_w - 142) / 2)))
    cell_h = 112
    h_gap = available_w - 2 * cell_w
    note_top = y1 - 66
    v_gap = max(54, note_top - map_y0 - 2 * cell_h - 12)
    cells = {
        "C": (map_x0, map_y0),
        "D": (map_x0 + cell_w + h_gap, map_y0),
        "A": (map_x0, map_y0 + cell_h + v_gap),
        "B": (map_x0 + cell_w + h_gap, map_y0 + cell_h + v_gap),
    }
    meta = {
        "A": ("selected sample", "isolated + tight target", ft.PHOTON_DARK, (255, 248, 229)),
        "B": ("non-isolated + tight", "true photon fails isolation", ft.SPHENIX_BLUE, (238, 247, 252)),
        "C": ("isolated + non-tight", "true photon fails tight ID", ft.TEAL, (239, 249, 250)),
        "D": ("non-isolated + non-tight", "true photon fails both axes", ft.MUTED, (245, 247, 250)),
    }
    for region, (cx, cy) in cells.items():
        label, detail, color, fill = meta[region]
        draw.rounded_rectangle((cx, cy, cx + cell_w, cy + cell_h), radius=8, fill=(*fill, 255), outline=(*color, 230), width=3)
        draw.text((cx + 18, cy + 14), region, font=ft.font(ft.TIMES_BOLD, 41), fill=color)
        ft.draw_wrapped(draw, label, (cx + 78, cy + 14), cell_w - 94, ft.font(ft.TIMES_BOLD, 24), fill=ft.INK, line_gap=2)
        ft.draw_wrapped(draw, detail, (cx + 78, cy + 64), cell_w - 94, ft.font(ft.TIMES, 22), fill=ft.MUTED, line_gap=2)

    ax, ay = cells["A"]
    bx, by = cells["B"]
    cx, cy = cells["C"]
    dx, dy = cells["D"]
    c_start = (ax + cell_w // 2, ay - 10)
    c_end = (cx + cell_w // 2, cy + cell_h + 10)
    b_start = (ax + cell_w + 8, ay + cell_h // 2 + 4)
    b_end = (bx - 12, by + cell_h // 2 + 4)
    d_start = (ax + cell_w + 14, ay + 26)
    d_end = (dx + 4, dy + cell_h - 4)
    ft.draw_arrow(draw, c_start, c_end, fill=ft.TEAL, width=4)
    ft.draw_arrow(draw, b_start, b_end, fill=ft.SPHENIX_BLUE, width=4)
    ft.draw_arrow(draw, d_start, d_end, fill=ft.MUTED, width=4)
    b_angle = math.degrees(math.atan2(b_end[1] - b_start[1], b_end[0] - b_start[0]))
    d_angle = math.degrees(math.atan2(d_end[1] - d_start[1], d_end[0] - d_start[0]))
    draw_leakage_label(base, (c_start[0] + 74, (c_start[1] + c_end[1]) // 2), "C leakage", ft.TEAL)
    draw_leakage_label(base, ((b_start[0] + b_end[0]) // 2, b_start[1] - 18), "B leakage", ft.SPHENIX_BLUE, b_angle)
    draw_leakage_label(base, ((d_start[0] + d_end[0]) // 2 + 4, (d_start[1] + d_end[1]) // 2 - 16), "D leakage", ft.MUTED, d_angle)

    note = (x0 + 8, y1 - 76, x1 - 8, y1 - 12)
    icon_x, icon_y = note[0] + 134, note[1] + 15
    draw.rounded_rectangle((icon_x, icon_y, icon_x + 34, icon_y + 34), radius=5, fill=(255, 255, 255, 255), outline=(*ft.SPHENIX_BLUE, 255), width=2)
    draw.line((icon_x + 9, icon_y + 11, icon_x + 25, icon_y + 11), fill=(*ft.SPHENIX_BLUE, 220), width=3)
    draw.line((icon_x + 9, icon_y + 19, icon_x + 24, icon_y + 19), fill=(*ft.SPHENIX_BLUE, 190), width=3)
    draw.line((icon_x + 9, icon_y + 27, icon_x + 20, icon_y + 27), fill=(*ft.SPHENIX_BLUE, 165), width=3)
    note_font = ft.font(ft.TIMES_BOLD, 27)
    body_font = ft.font(ft.TIMES, 25)
    line1 = "Note: ABCD purity requires isolation-tight ID complementarity;"
    line2 = "isolation is not used in BDT training."
    text_left = icon_x + 44
    text_right = note[2] - 18
    line1_w = ft.text_box(draw, line1, note_font)[0]
    line2_w = ft.text_box(draw, line2, body_font)[0]
    draw.text((text_left + max(0, (text_right - text_left - line1_w) // 2), note[1] + 8), line1, font=note_font, fill=ft.BLUE)
    draw.text((text_left + max(0, (text_right - text_left - line2_w) // 2), note[1] + 37), line2, font=body_font, fill=ft.BLUE)


def draw_purity_plot(base: Image.Image, box: tuple[int, int, int, int], title: str = "Leakage-corrected purity") -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    ft.shadow(base, box, radius=12)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((box[0] + 18, box[1] + 24, box[0] + 30, box[3] - 24), radius=6, fill=(*ft.SPHENIX_BLUE, 255))
    header_x = box[0] + 54
    header_y = box[1] + 24
    title_font = ft.font(ft.TIMES_BOLD, 33)
    body_font = ft.font(ft.TIMES, 27)
    draw.text((header_x, header_y), title, font=title_font, fill=ft.INK)
    title_w = ft.text_box(draw, title, title_font)[0]
    draw.text((header_x + title_w + 10, header_y + 8), "= signal fraction carried into the corrected yield.", font=body_font, fill=ft.INK)
    draw.line((box[0] + 54, box[1] + 82, box[2] - 34, box[1] + 82), fill=(221, 228, 236), width=2)

    plot = trim_white_margins(Image.open(ft.figure_path("fig5_purity")).convert("RGBA"), tolerance=10, pad=6)
    content = (box[0] + 54, box[1] + 100, box[2] - 34, box[3] - 28)
    plot_box = (content[0] + 32, content[1] + 6, content[0] + 548, content[3] - 6)
    ft.shadow(base, plot_box, radius=10)
    draw.rounded_rectangle(
        (plot_box[0] - 8, plot_box[1] - 8, plot_box[2] + 8, plot_box[3] + 8),
        radius=8,
        fill=(255, 255, 255, 255),
        outline=(216, 226, 235, 255),
        width=2,
    )
    ft.paste_fit(base, plot, plot_box, anchor="center")
    draw_purity_reading_key(base, content, plot_box)


def draw_purity_reading_key(base: Image.Image, content: tuple[int, int, int, int], plot_box: tuple[int, int, int, int]) -> None:
    key_x0 = plot_box[2] + 34
    key_x1 = content[2] - 8
    gap = 20
    key_h = ((content[3] - 18) - (content[1] + 18) - gap) // 2
    raw_box = (key_x0, content[1] + 18, key_x1, content[1] + 18 + key_h)
    corr_box = (key_x0, raw_box[3] + gap, key_x1, content[3] - 18)

    draw_purity_key_card(
        base,
        raw_box,
        (32, 32, 32),
        "Raw sideband purity",
        "Black points show the ABCD sideband estimate before correcting signal leakage into B/C/D.",
        fill=(246, 247, 249, 255),
    )
    draw_purity_key_card(
        base,
        corr_box,
        ft.SPHENIX_BLUE,
        "Leakage-corrected purity",
        "Blue points remove true-photon leakage using MC. This corrected purity curve is applied to Region A in the yield correction.",
        fill=(236, 247, 255, 255),
    )


def draw_purity_key_card(
    base: Image.Image,
    box: tuple[int, int, int, int],
    color: tuple[int, int, int],
    title: str,
    body: str,
    fill: tuple[int, int, int, int],
) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    draw.rounded_rectangle(box, radius=12, fill=fill, outline=(205, 221, 235, 255), width=2)
    draw.text((box[0] + 24, box[1] + 22), title, font=ft.font(ft.TIMES_BOLD, 29), fill=color)
    ft.draw_wrapped(
        draw,
        body,
        (box[0] + 24, box[1] + 68),
        box[2] - box[0] - 48,
        ft.font(ft.TIMES, 25),
        fill=ft.INK,
        line_gap=4,
    )


def draw_purity_logic_box(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    gap = 14
    raw_h = 188
    raw_box = (box[0], box[1], box[2], box[1] + raw_h)
    lc_box = (box[0], raw_box[3] + gap, box[2], box[3])
    draw_definition_section(
        base,
        raw_box,
        "Raw sideband estimate",
        [
            r"$N_{\mathrm{sig}}^{A,\mathrm{raw}}=N_{\mathrm{raw}}^A-N_{\mathrm{bkg}}^{\mathrm{ABCD}}$",
            r"$N_{\mathrm{bkg}}^{\mathrm{ABCD}}=N_{\mathrm{raw}}^B(N_{\mathrm{raw}}^C/N_{\mathrm{raw}}^D)$",
        ],
        [
            ("N_A raw", "selected A"),
            ("N_bkg", "ABCD estimate"),
            ("B/C/D", "control regions"),
        ],
        fontsize=26,
    )
    draw_definition_section(
        base,
        lc_box,
        "Leakage-corrected estimate",
        [
            r"$b_X=N_{\mathrm{raw}}^X-f_X^{\mathrm{MC}}N_{\mathrm{sig}}^A,\quad X=B,C,D$",
            r"$N_{\mathrm{sig}}^{A,\mathrm{LC}}=N_{\mathrm{raw}}^A-b_B(b_C/b_D)$",
            r"$P_{\mathrm{LC}}(E_T)=N_{\mathrm{sig}}^{A,\mathrm{LC}}/N_{\mathrm{raw}}^A$",
        ],
        [
            ("f_MC", "MC leakage fraction"),
            (r"$b_X$", "leakage-corrected sideband"),
            ("P_LC", "purity carried forward"),
        ],
        fontsize=25,
    )


def draw_definition_section(
    base: Image.Image,
    section_box: tuple[int, int, int, int],
    title: str,
    eqs: list[str],
    terms: list[tuple[str, str]],
    fontsize: int,
) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = section_box
    pad = 12
    draw.rounded_rectangle(
        section_box,
        radius=9,
        fill=(255, 255, 255, 248),
        outline=(216, 226, 235, 255),
        width=2,
    )
    draw.text((x0 + pad + 8, y0 + 9), title, font=ft.font(ft.TIMES_BOLD, 24), fill=ft.INK)
    draw.line((x0 + pad + 8, y0 + 44, x1 - pad - 8, y0 + 44), fill=(228, 234, 240), width=2)

    content_top = y0 + 54
    content_bottom = y1 - 12
    content_h = content_bottom - content_top
    split_x = x0 + int((x1 - x0) * 0.63)
    draw.line((split_x, content_top + 2, split_x, content_bottom - 2), fill=(234, 239, 244), width=2)

    eq_img = math_block_image(eqs, fontsize=fontsize, width=4.8, line_step=0.31)
    paste_math_left(base, eq_img, (x0 + pad + 8, content_top, split_x - 16, content_bottom))

    term_font = ft.font(ft.TIMES_BOLD, 18)
    def_font = ft.font(ft.TIMES, 18)
    pill_x = split_x + 12
    row_h = max(37, content_h // max(1, len(terms)))
    y = content_top + max(0, (content_h - row_h * len(terms)) // 2)
    for label, meaning in terms:
        if label.startswith("$"):
            label_img = math_block_image([label], fontsize=18, width=1.15, line_step=0.30)
            label_w = max(64, min(105, label_img.width))
            pill = (pill_x, y + 2, pill_x + label_w + 16, y + 30)
            draw.rounded_rectangle(pill, radius=6, fill=(246, 249, 252, 255), outline=(215, 225, 234, 255), width=1)
            paste_math_fit(base, label_img, (pill_x + 8, y + 4, pill_x + label_w + 8, y + 27))
        else:
            label_w = ft.text_box(draw, label, term_font)[0]
            pill = (pill_x, y + 2, pill_x + max(label_w + 16, 78), y + 30)
            draw.rounded_rectangle(pill, radius=6, fill=(246, 249, 252, 255), outline=(215, 225, 234, 255), width=1)
            draw.text((pill_x + 8, y + 5), label, font=term_font, fill=ft.INK)
        ft.draw_wrapped(
            draw,
            f"= {meaning}",
            (pill[2] + 10, y + 5),
            max(80, x1 - pill[2] - 24),
            def_font,
            fill=ft.INK,
            line_gap=1,
        )
        y += row_h


def draw_final_purity_carry_forward(base: Image.Image, box: tuple[int, int, int, int], y: int) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    pad = 18
    draw.rounded_rectangle((box[0] + pad, y, box[2] - pad, box[3] - 14), radius=9, fill=(247, 250, 252, 255), outline=(216, 226, 235, 255), width=2)
    draw.text((box[0] + pad + 16, y + 7), "Purity carried forward", font=ft.font(ft.TIMES_BOLD, 20), fill=ft.INK)
    eq_img = math_block_image([r"$P_{\mathrm{LC}}(E_T)=N_{\mathrm{sig}}^{A,\mathrm{LC}}/N_{\mathrm{raw}}^A$"], fontsize=29, width=5.9, line_step=0.32)
    paste_math_left(base, eq_img, (box[0] + pad + 16, y + 29, box[2] - pad - 8, box[3] - 17))


def paste_scaled_formula(
    base: Image.Image,
    renderer,
    xy: tuple[int, int],
    scale: float,
    canvas_size: tuple[int, int],
) -> None:
    layer = Image.new("RGBA", canvas_size, (255, 255, 255, 0))
    renderer(ImageDraw.Draw(layer, "RGBA"), 18, 24)
    bbox = layer.getbbox()
    if not bbox:
        return
    crop = layer.crop(bbox)
    crop = crop.resize((max(1, int(crop.width * scale)), max(1, int(crop.height * scale))), Image.Resampling.LANCZOS)
    base.alpha_composite(crop, xy)


def draw_f_mc_var(draw: ImageDraw.ImageDraw, x: int, y: int, sub: str = "X", size: int = 27, fill=ft.INK) -> int:
    base = ft.font(ft.TIMES_ITALIC, size)
    small = ft.font(ft.TIMES, max(13, int(size * 0.55)))
    draw.text((x, y), "f", font=base, fill=fill)
    bw, _ = ft.text_box(draw, "f", base)
    draw.text((x + bw - 1, y - int(size * 0.36)), "MC", font=small, fill=fill)
    draw.text((x + bw - 1, y + int(size * 0.48)), sub, font=small, fill=fill)
    return x + bw + max(ft.text_box(draw, "MC", small)[0], ft.text_box(draw, sub, small)[0]) + 8


def draw_b_sub(draw: ImageDraw.ImageDraw, x: int, y: int, sub: str = "X", size: int = 28, fill=ft.INK) -> int:
    base = ft.font(ft.TIMES_ITALIC, size)
    small = ft.font(ft.TIMES, max(13, int(size * 0.55)))
    draw.text((x, y), "b", font=base, fill=fill)
    bw, _ = ft.text_box(draw, "b", base)
    draw.text((x + bw - 1, y + int(size * 0.48)), sub, font=small, fill=fill)
    return x + bw + ft.text_box(draw, sub, small)[0] + 8


def draw_raw_sideband_equation(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 28) -> None:
    x = c3.draw_n_var(draw, x, y, "A", "sig", size=size)
    x = c3.math_text(draw, x, y + 2, " = ", size=size)
    x = c3.draw_n_var(draw, x, y, "A", "raw", size=size)
    x = c3.math_text(draw, x, y + 2, " − ", size=size)
    x = c3.draw_n_var(draw, x, y, "B", "raw", size=size)
    x = c3.math_text(draw, x, y + 2, " (", size=size)
    frac_x = x + 4
    num_end = c3.draw_n_var(draw, frac_x + 10, y - 11, "C", "raw", size=size - 4)
    frac_w = max(78, num_end - frac_x + 8)
    draw.line((frac_x, y + 35, frac_x + frac_w, y + 35), fill=ft.INK, width=2)
    c3.draw_n_var(draw, frac_x + 10, y + 38, "D", "raw", size=size - 4)
    draw.text((frac_x + frac_w + 6, y - 3), ")", font=ft.font(ft.TIMES, size + 15), fill=ft.INK)


def draw_tag(draw: ImageDraw.ImageDraw, center_x: int, y: int, text: str, color: tuple[int, int, int], size: int = 19) -> tuple[int, int, int, int]:
    fnt = ft.font(ft.TIMES_BOLD, size)
    tw, th = ft.text_box(draw, text, fnt)
    x0 = int(center_x - tw / 2 - 14)
    box = (x0, y, x0 + tw + 28, y + th + 10)
    draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 248), outline=(*color, 215), width=2)
    draw.text((x0 + 14, y + 4), text, font=fnt, fill=color)
    return box


def draw_raw_sideband_equation_annotated(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 32) -> None:
    x = c3.draw_n_var(draw, x, y, "A", "sig", size=size)
    x = c3.math_text(draw, x, y + 2, " = ", size=size)

    selected_start = x
    x = c3.draw_n_var(draw, x, y, "A", "raw", size=size)
    selected_end = x

    x = c3.math_text(draw, x, y + 2, " − ", size=size)

    sideband_start = x
    x = c3.draw_n_var(draw, x, y, "B", "raw", size=size)
    x = c3.math_text(draw, x, y + 2, " (", size=size)
    x = c3.draw_n_var(draw, x, y, "C", "raw", size=size)
    x = c3.math_text(draw, x, y + 2, " / ", size=size)
    x = c3.draw_n_var(draw, x, y, "D", "raw", size=size)
    x = c3.math_text(draw, x, y + 2, ")", size=size)
    sideband_end = x

    term_top = y - int(size * 0.50)
    term_bottom = y + int(size * 1.18)
    selected_box = (selected_start - 10, term_top, selected_end + 10, term_bottom)
    sideband_box = (sideband_start - 10, term_top, sideband_end + 10, term_bottom)
    draw.rounded_rectangle(selected_box, radius=18, outline=(*ft.PHOTON_DARK, 230), width=4)
    draw.rounded_rectangle(sideband_box, radius=18, outline=(*ft.TEAL, 230), width=4)

    selected_tag = draw_tag(draw, (selected_start + selected_end) // 2, y - 57, "selected Region A", ft.PHOTON_DARK, size=20)
    sideband_tag = draw_tag(draw, (sideband_start + sideband_end) // 2, y - 57, "sideband estimate", ft.TEAL, size=20)
    draw.line(((selected_tag[0] + selected_tag[2]) // 2, selected_tag[3], (selected_box[0] + selected_box[2]) // 2, selected_box[1]), fill=(*ft.PHOTON_DARK, 190), width=2)
    draw.line(((sideband_tag[0] + sideband_tag[2]) // 2, sideband_tag[3], (sideband_box[0] + sideband_box[2]) // 2, sideband_box[1]), fill=(*ft.TEAL, 190), width=2)


def draw_corrected_abcd_equation(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 27) -> None:
    line1_x = x
    line1_x = c3.draw_n_var(draw, line1_x, y, "A", "signal", size=size)
    line1_x = c3.math_text(draw, line1_x, y + 2, " = ", size=size)
    line1_x = c3.draw_n_var(draw, line1_x, y, "A", "raw", size=size)
    line1_x = c3.math_text(draw, line1_x, y + 2, " − [(", size=size)
    line1_x = c3.draw_n_var(draw, line1_x, y, "B", "raw", size=size)
    line1_x = c3.math_text(draw, line1_x, y + 2, " − ", size=size)
    line1_x = draw_f_mc_var(draw, line1_x, y, "B", size=size, fill=ft.TEAL)
    line1_x = c3.draw_n_var(draw, line1_x + 1, y, "A", "signal", size=size)
    c3.math_text(draw, line1_x, y + 2, ")", size=size)

    line2_y = y + 58
    line2_x = x + 142
    line2_x = c3.math_text(draw, line2_x, line2_y + 2, "× ", size=size)
    frac_x = line2_x + 4
    num_x = frac_x + 18
    num_end = c3.draw_n_var(draw, num_x, line2_y - 24, "C", "raw", size=size)
    num_end = c3.math_text(draw, num_end, line2_y - 22, " − ", size=size)
    num_end = draw_f_mc_var(draw, num_end, line2_y - 24, "C", size=size, fill=ft.TEAL)
    num_end = c3.draw_n_var(draw, num_end + 1, line2_y - 24, "A", "signal", size=size)
    frac_w = max(270, num_end - frac_x + 22)
    draw.line((frac_x, line2_y + 20, frac_x + frac_w, line2_y + 20), fill=ft.INK, width=2)
    den_x = frac_x + 18
    den_end = c3.draw_n_var(draw, den_x, line2_y + 25, "D", "raw", size=size)
    den_end = c3.math_text(draw, den_end, line2_y + 27, " − ", size=size)
    den_end = draw_f_mc_var(draw, den_end, line2_y + 25, "D", size=size, fill=ft.TEAL)
    c3.draw_n_var(draw, den_end + 1, line2_y + 25, "A", "signal", size=size)
    draw.text((frac_x + frac_w + 12, line2_y - 21), "]", font=ft.font(ft.TIMES, size + 34), fill=ft.INK)


def math_lines_image(lines: list[str], fontsize: int = 30, dpi: int = 240) -> Image.Image:
    fig = Figure(figsize=(6.8, 2.75), dpi=dpi)
    fig.patch.set_alpha(0)
    canvas = FigureCanvasAgg(fig)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.axis("off")
    if len(lines) == 3:
        placements = [(0.04, 0.78), (0.19, 0.49), (0.21, 0.20)]
    else:
        placements = [(0.02, 0.70), (0.16, 0.28)]
    for line, (tx, ty) in zip(lines, placements):
        fig.text(tx, ty, line, fontsize=fontsize, family="serif", color=tuple(v / 255 for v in ft.INK), va="center")
    canvas.draw()
    img = Image.frombuffer("RGBA", canvas.get_width_height(), canvas.buffer_rgba(), "raw", "RGBA", 0, 1).copy()
    return trim_white_margins(img, tolerance=2, pad=4)


def math_block_image(lines: list[str], fontsize: int = 26, width: float = 4.5, line_step: float = 0.32, dpi: int = 240) -> Image.Image:
    height = 0.70 + line_step * max(1, len(lines))
    fig = Figure(figsize=(width, height), dpi=dpi)
    fig.patch.set_alpha(0)
    canvas = FigureCanvasAgg(fig)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.axis("off")
    top = 0.82
    for i, line in enumerate(lines):
        fig.text(0.03, top - i * line_step, line, fontsize=fontsize, family="serif", color=tuple(v / 255 for v in ft.INK), va="center")
    canvas.draw()
    img = Image.frombuffer("RGBA", canvas.get_width_height(), canvas.buffer_rgba(), "raw", "RGBA", 0, 1).copy()
    return trim_white_margins(img, tolerance=2, pad=4)


def paste_math_fit(base: Image.Image, img: Image.Image, box: tuple[int, int, int, int]) -> None:
    bw = box[2] - box[0]
    bh = box[3] - box[1]
    scale = min(bw / img.width, bh / img.height)
    resized = img.resize((max(1, int(img.width * scale)), max(1, int(img.height * scale))), Image.Resampling.LANCZOS)
    base.alpha_composite(resized, (box[0] + (bw - resized.width) // 2, box[1] + (bh - resized.height) // 2))


def paste_math_left(base: Image.Image, img: Image.Image, box: tuple[int, int, int, int]) -> None:
    bw = box[2] - box[0]
    bh = box[3] - box[1]
    scale = min(bw / img.width, bh / img.height)
    resized = img.resize((max(1, int(img.width * scale)), max(1, int(img.height * scale))), Image.Resampling.LANCZOS)
    base.alpha_composite(resized, (box[0], box[1] + (bh - resized.height) // 2))


def paste_corrected_abcd_equation(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    lines = [
        r"$N^A_{\mathrm{signal}}=N^A_{\mathrm{raw}}-\Big[$",
        r"$(N^B_{\mathrm{raw}}-f^{B,\mathrm{MC}}N^A_{\mathrm{signal}})$",
        r"$\times\frac{N^C_{\mathrm{raw}}-f^{C,\mathrm{MC}}N^A_{\mathrm{signal}}}{N^D_{\mathrm{raw}}-f^{D,\mathrm{MC}}N^A_{\mathrm{signal}}}\Big]$",
    ]
    paste_math_fit(base, math_lines_image(lines, fontsize=34), box)


def draw_leakage_correction_equation(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 28) -> None:
    x = draw_b_sub(draw, x, y, "X", size=size, fill=ft.TEAL)
    x = c3.math_text(draw, x, y + 2, " = ", size=size)
    x = c3.draw_n_var(draw, x, y, "X", "raw", size=size)
    x = c3.math_text(draw, x, y + 2, " − ", size=size)
    x = draw_f_mc_var(draw, x, y, "X", size=size, fill=ft.TEAL)
    x = c3.draw_n_var(draw, x + 2, y, "A", "sig", size=size)
    tail_x = x + 20
    tail_font = ft.font(ft.TIMES_BOLD, size - 6)
    symbol_font = ft.font(ARIAL_UNICODE, size - 6)
    draw.text((tail_x, y + 4), "X ", font=tail_font, fill=ft.MUTED)
    tail_x += ft.text_box(draw, "X ", tail_font)[0]
    draw.text((tail_x, y + 4), "∈", font=symbol_font, fill=ft.MUTED)
    tail_x += ft.text_box(draw, "∈", symbol_font)[0] + 4
    draw.text((tail_x, y + 4), "{B,C,D}", font=tail_font, fill=ft.MUTED)


def draw_purity_equation(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 30) -> None:
    x = c3.math_text(draw, x, y, "P(", size=size)
    x = c3.math_text(draw, x, y, "E", size=size, fill=ft.BLUE, bold=True)
    x = c3.math_text(draw, x - 2, y + int(size * 0.50), "T", size=int(size * 0.55), fill=ft.BLUE, bold=True)
    x = c3.math_text(draw, x + 7, y, ") = ", size=size)
    frac_x = x
    num_end = c3.draw_n_var(draw, frac_x + 18, y - 13, "A", "sig", size=size - 2)
    frac_w = max(106, num_end - frac_x + 16)
    draw.line((frac_x, y + 37, frac_x + frac_w, y + 37), fill=ft.INK, width=2)
    c3.draw_n_var(draw, frac_x + 18, y + 39, "A", "raw", size=size - 2)


def draw_purity_equation_inline(draw: ImageDraw.ImageDraw, x: int, y: int, size: int = 35) -> None:
    x = c3.math_text(draw, x, y, "P(", size=size)
    x = c3.math_text(draw, x, y, "E", size=size, fill=ft.BLUE, bold=True)
    x = c3.math_text(draw, x - 2, y + int(size * 0.50), "T", size=int(size * 0.55), fill=ft.BLUE, bold=True)
    x = c3.math_text(draw, x + 7, y, ") = ", size=size)
    x = c3.draw_n_var(draw, x, y, "A", "sig", size=size)
    x = c3.math_text(draw, x + 2, y + 2, " / ", size=size)
    c3.draw_n_var(draw, x, y, "A", "raw", size=size)


def draw_minimal_math_strip(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    ft.shadow(base, box, radius=12)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.text((box[0] + 30, box[1] + 18), "Leakage Correction", font=ft.font(ft.TIMES_BOLD, 38), fill=ft.INK)

    x0, y0, x1, y1 = box
    gap = 22
    usable_w = x1 - x0 - 60 - 2 * gap
    col_ws = [680, 900, usable_w - 680 - 900]
    cols = [
        ("1", "Raw sideband estimate", ft.PHOTON_DARK),
        ("2", "Leakage correction", ft.TEAL),
        ("3", "Purity carried forward", ft.SPHENIX_BLUE),
    ]
    card_boxes = []
    cx = x0 + 30
    for idx, (num, title, color) in enumerate(cols):
        col_w = col_ws[idx]
        cy = y0 + 72
        draw.rounded_rectangle((cx, cy, cx + col_w, y1 - 24), radius=10, fill=(248, 251, 253, 255), outline=(*color, 175), width=2)
        draw.rounded_rectangle((cx, cy, cx + 12, y1 - 24), radius=4, fill=(*color, 255))
        draw.ellipse((cx + 28, cy + 16, cx + 74, cy + 62), fill=(*color, 255))
        draw.text((cx + 51, cy + 39), num, font=ft.font(ft.TIMES_BOLD, 27), fill=(255, 255, 255), anchor="mm")
        draw.text((cx + 92, cy + 19), title, font=ft.font(ft.TIMES_BOLD, 29), fill=ft.INK)
        card_boxes.append((cx, cy, cx + col_w, y1 - 24))
        cx += col_w + gap

    draw_raw_sideband_equation_annotated(draw, card_boxes[0][0] + 74, card_boxes[0][1] + 104, size=32)
    paste_corrected_abcd_equation(base, (card_boxes[1][0] + 42, card_boxes[1][1] + 78, card_boxes[1][2] - 34, card_boxes[1][3] - 14))
    draw_purity_equation_inline(draw, card_boxes[2][0] + 74, card_boxes[2][1] + 104, size=40)


def draw_bottom_flow(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    ft.shadow(base, box, radius=12)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.text((box[0] + 42, box[1] + 28), "Corrected-yield flow", font=ft.font(ft.TIMES_BOLD, 42), fill=ft.INK)
    ft.draw_wrapped(
        draw,
        "Leakage-corrected purity is applied to Region A; efficiency and unfolding move the yield to particle level.",
        (box[0] + 42, box[1] + 82),
        box[2] - box[0] - 84,
        ft.font(ft.TIMES_ITALIC, 31),
        fill=ft.MUTED,
        line_gap=5,
    )
    nodes = [
        ("Region A\ncandidates", ft.PHOTON_DARK, 340),
        ("apply leakage-corrected\npurity curve", ft.SPHENIX_BLUE, 372),
        ("apply\nefficiency", ft.TEAL, 270),
        ("unfold detector\nresponse", ft.PHOTON, 356),
        ("particle-level yield\n/ cross section", ft.BLUE, 404),
    ]
    arrow_w = 46
    gap = 20
    total_w = sum(w for _, _, w in nodes) + (len(nodes) - 1) * (arrow_w + gap)
    x = box[0] + (box[2] - box[0] - total_w) // 2
    y = box[1] + 174
    h = 112
    for idx, (label, color, w) in enumerate(nodes):
        draw.rounded_rectangle((x, y, x + w, y + h), radius=12, fill=(247, 250, 252, 255), outline=(*color, 230), width=3)
        draw.rounded_rectangle((x, y, x + 15, y + h), radius=5, fill=(*color, 255))
        lines = label.splitlines()
        line_h = 36
        start_y = y + (h - len(lines) * line_h) // 2 - 1
        for j, line in enumerate(lines):
            fnt = ft.font(ft.TIMES_BOLD, 30)
            tw, _ = ft.text_box(draw, line, fnt)
            draw.text((x + 15 + (w - 15 - tw) / 2, start_y + j * line_h), line, font=fnt, fill=ft.INK)
        if idx < len(nodes) - 1:
            ft.draw_arrow(draw, (x + w + 10, y + h // 2), (x + w + arrow_w, y + h // 2), fill=(151, 169, 188), width=6)
        x += w + arrow_w + gap


def render_candidate(kind: str) -> Image.Image:
    img = base_slide()
    if kind == "wide_plot":
        leak_box = (132, 306, 1168, 890)
        purity_box = (1210, 306, 2390, 890)
    elif kind == "balanced":
        leak_box = (132, 306, 1168, 890)
        purity_box = (1210, 306, 2390, 890)
    else:
        leak_box = (132, 306, 1206, 896)
        purity_box = (1248, 306, 2390, 896)

    draw_plain_accent_card(img, leak_box, ft.PHOTON_DARK)
    draw_leakage_map(img, (leak_box[0] + 54, leak_box[1] + 26, leak_box[2] - 34, leak_box[3] - 22))
    draw_purity_plot(img, purity_box)
    flow_y = 914 if kind != "large_diagram" else 920
    draw_bottom_flow(img, (132, flow_y, 2390, 1294))
    ft.draw_hp2026_identity_footer(img)
    return img


def write_contact_sheet(paths: list[Path]) -> Path:
    sheet = Image.new("RGB", (640 * len(paths), 420), (255, 255, 255))
    draw = ImageDraw.Draw(sheet)
    for i, path in enumerate(paths):
        thumb = Image.open(path).convert("RGB").resize((640, 360), Image.Resampling.LANCZOS)
        x = i * 640
        sheet.paste(thumb, (x, 58))
        draw.text((x + 18, 18), path.stem, font=ft.font(ft.TIMES, 22), fill=ft.INK)
    out = OUTDIR / "slide16_yield_flow_redesign_internal_contact_sheet.png"
    sheet.save(out)
    return out


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    CANDIDATE_DIR.mkdir(parents=True, exist_ok=True)

    variants = {
        "candidate01_wide_plot": "wide_plot",
        "candidate02_balanced_selected": "balanced",
        "candidate03_large_diagram": "large_diagram",
    }
    paths: list[Path] = []
    for stem, kind in variants.items():
        img = render_candidate(kind)
        out = CANDIDATE_DIR / f"{stem}.png"
        img.convert("RGB").save(out, "PNG")
        paths.append(out)

    selected = paths[1]
    shutil.copyfile(selected, FINAL_PNG)
    contact = write_contact_sheet(paths)
    MANIFEST.write_text(
        json.dumps(
            {
                "generated_at": datetime.now().isoformat(timespec="seconds"),
                "script": str(Path(__file__).resolve()),
                "source_generator_found": str((SCRIPT_DIR / "make_hp2026_closing_three_candidates.py").resolve()),
                "source_function": "slide10() / draw_purity_equation_panel() / draw_leakage_meaning_card()",
                "source_plot_assets": {
                    "purity": str(ft.figure_path("fig5_purity")),
                },
                "deck_mutation": "none; local PNG prototype only",
                "candidate_paths": [str(p) for p in paths],
                "selected": str(FINAL_PNG),
                "contact_sheet": str(contact),
                "selection_rationale": "two-panel upper row: the sideband/leakage diagram explains the control regions, while the purity panel places the corrected-purity plot next to compact raw and leakage-corrected ABCD definitions; the freed lower region becomes a larger corrected-yield flow ribbon.",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(FINAL_PNG)


if __name__ == "__main__":
    main()

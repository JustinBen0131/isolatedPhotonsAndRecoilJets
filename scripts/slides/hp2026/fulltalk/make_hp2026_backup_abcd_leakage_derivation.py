#!/usr/bin/env python3
"""Local-only HP2026 backup slide explaining the ABCD leakage correction."""

from __future__ import annotations

import json
import math
import sys
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw
from matplotlib import rc_context
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure


ROOT = next(
    p
    for p in Path(__file__).resolve().parents
    if (p / "README.md").exists() and (p / "scripts").exists() and (p / "src").exists()
)
SCRIPT_DIR = ROOT / "scripts/slides/hp2026/fulltalk"
sys.path.insert(0, str(SCRIPT_DIR))
COMMON_DIR = ROOT / "scripts/slides/common"
sys.path.insert(0, str(COMMON_DIR))

import make_hp2026_fulltalk_candidates as ft  # noqa: E402
import make_hp2026_slide16_yield_flow_redesign as slide16  # noqa: E402


OUTDIR = ROOT / "outputs/manual-20260610-backup-abcd-leakage-derivation"
PNG = OUTDIR / "backup_abcd_leakage_correction_derivation.png"
SCRIPT_MD = OUTDIR / "backup_abcd_leakage_correction_derivation_script.md"
MANIFEST = OUTDIR / "backup_abcd_leakage_correction_derivation_manifest.json"
HEADER_JSON = PNG.with_suffix(".header.json")

TITLE = "Backup: ABCD leakage correction"
SUBTITLE = "Subtract true-photon leakage from B/C/D before computing the Region A purity."

HEADER_SPEC = {
    "deck": "hp2026_main_talk",
    "title_font_size": 86,
    "subtitle_font_size": 56,
    "title_xy": [132, 76],
    "subtitle_xy": [136, 190],
    "divider_y": 286,
}


def text_center_y(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], text: str, font, fill) -> None:
    tw, th = ft.text_box(draw, text, font)
    draw.text((box[0] + (box[2] - box[0] - tw) / 2, box[1] + (box[3] - box[1] - th) / 2 - 2), text, font=font, fill=fill)


def draw_slide_shell() -> Image.Image:
    img = Image.new("RGBA", (ft.W, ft.H), (*ft.SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, ft.W, ft.H), fill=(*ft.SOFT_BG, 255))
    draw.rectangle((0, 0, ft.W, 22), fill=(*ft.SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, ft.W, 30), fill=(*ft.PHOTON, 255))
    draw.text(tuple(HEADER_SPEC["title_xy"]), TITLE, font=ft.font(ft.TIMES_BOLD, HEADER_SPEC["title_font_size"]), fill=ft.INK)
    draw.text(tuple(HEADER_SPEC["subtitle_xy"]), SUBTITLE, font=ft.font(ft.TIMES_ITALIC, HEADER_SPEC["subtitle_font_size"]), fill=ft.MUTED)
    draw.line((132, HEADER_SPEC["divider_y"], ft.W - 132, HEADER_SPEC["divider_y"]), fill=(221, 226, 232, 255), width=3)
    ft.add_top_right_sphenix_logo_like_slide2(img)
    return img


def draw_card(base: Image.Image, box: tuple[int, int, int, int], accent: tuple[int, int, int], title: str, subtitle: str = "") -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    ft.shadow(base, box, radius=12)
    draw.rounded_rectangle(box, radius=13, fill=(255, 255, 255, 255), outline=(*ft.PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((box[0], box[1], box[0] + 14, box[3]), radius=7, fill=(*accent, 255))
    draw.text((box[0] + 58, box[1] + 26), title, font=ft.font(ft.TIMES_BOLD, 41), fill=ft.INK)
    divider_y = box[1] + 138
    if subtitle:
        ft.draw_wrapped(draw, subtitle, (box[0] + 58, box[1] + 80), box[2] - box[0] - 104, ft.font(ft.TIMES_ITALIC, 27), fill=ft.MUTED, line_gap=4)
        divider_y = box[1] + 158
    draw.line((box[0] + 58, divider_y, box[2] - 36, divider_y), fill=(221, 228, 236, 255), width=2)


def paste_math(base: Image.Image, lines: list[str], box: tuple[int, int, int, int], fontsize: int = 34, width: float = 7.4) -> None:
    img = slide16.math_block_image(lines, fontsize=fontsize, width=width, line_step=0.34)
    slide16.paste_math_fit(base, img, box)


def math_line_image(line: str, fontsize: int = 34, width: float = 11.0, dpi: int = 240) -> Image.Image:
    with rc_context({"mathtext.fontset": "stix", "font.family": "STIXGeneral"}):
        fig = Figure(figsize=(width, 0.78), dpi=dpi)
        fig.patch.set_alpha(0)
        canvas = FigureCanvasAgg(fig)
        ax = fig.add_axes((0, 0, 1, 1))
        ax.axis("off")
        fig.text(0.015, 0.52, line, fontsize=fontsize, family="STIXGeneral", color=tuple(v / 255 for v in ft.INK), va="center")
        canvas.draw()
        img = Image.frombuffer("RGBA", canvas.get_width_height(), canvas.buffer_rgba(), "raw", "RGBA", 0, 1).copy()
    return slide16.trim_white_margins(img, tolerance=2, pad=4)


def paste_math_lines_centered(
    base: Image.Image,
    lines: list[str],
    box: tuple[int, int, int, int],
    fontsize: int = 34,
    line_gap: int = 8,
    max_scale: float = 1.0,
) -> None:
    rendered = [math_line_image(line, fontsize=fontsize) for line in lines]
    bw = box[2] - box[0]
    bh = box[3] - box[1]
    total_source_h = sum(img.height for img in rendered) + line_gap * (len(rendered) - 1)
    max_source_w = max(img.width for img in rendered)
    scale = min(bw / max_source_w, bh / total_source_h, max_scale)
    pasted: list[Image.Image] = []
    for img in rendered:
        resized = img.resize((max(1, int(img.width * scale)), max(1, int(img.height * scale))), Image.Resampling.LANCZOS)
        pasted.append(resized)
    total_h = sum(img.height for img in pasted) + line_gap * (len(pasted) - 1)
    y = int(box[1] + (bh - total_h) / 2)
    for img in pasted:
        x = int(box[0] + (bw - img.width) / 2)
        base.alpha_composite(img, (x, y))
        y += img.height + line_gap


def paste_math_line_fit(
    base: Image.Image,
    line: str,
    box: tuple[int, int, int, int],
    fontsize: int = 38,
    align: str = "left",
) -> None:
    img = math_line_image(line, fontsize=fontsize)
    bw = box[2] - box[0]
    bh = box[3] - box[1]
    scale = min(bw / img.width, bh / img.height, 1.0)
    resized = img.resize((max(1, int(img.width * scale)), max(1, int(img.height * scale))), Image.Resampling.LANCZOS)
    if align == "center":
        x = int(box[0] + (bw - resized.width) / 2)
    else:
        x = int(box[0])
    y = int(box[1] + (bh - resized.height) / 2)
    base.alpha_composite(resized, (x, y))


def draw_labeled_equation_rows(
    base: Image.Image,
    box: tuple[int, int, int, int],
    rows: list[tuple[str, str]],
    color: tuple[int, int, int],
    fontsize: int = 38,
    label_w: int = 238,
    label_size: int = 25,
) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 240), outline=(*color, 85), width=1)
    row_h = (y1 - y0) / len(rows)
    for idx, (label, eq) in enumerate(rows):
        ry0 = int(y0 + idx * row_h)
        ry1 = int(y0 + (idx + 1) * row_h)
        if idx:
            draw.line((x0 + 16, ry0, x1 - 16, ry0), fill=(226, 234, 240, 255), width=1)
        label_box = (x0 + 18, ry0, x0 + label_w, ry1)
        math_box = (x0 + label_w + 8, ry0 + 4, x1 - 18, ry1 - 4)
        text_center_y(draw, label_box, label, ft.font(ft.TIMES_BOLD, label_size), color)
        paste_math_line_fit(base, eq, math_box, fontsize=fontsize, align="left")


def draw_step_badge(draw: ImageDraw.ImageDraw, center: tuple[int, int], text: str, color: tuple[int, int, int]) -> None:
    r = 28
    draw.ellipse((center[0] - r, center[1] - r, center[0] + r, center[1] + r), fill=(*color, 255))
    draw.text(center, text, font=ft.font(ft.TIMES_BOLD, 31), fill=(255, 255, 255), anchor="mm")


def draw_leakage_label_clean(base: Image.Image, center: tuple[int, int], text: str, color: tuple[int, int, int], angle: float = 0) -> None:
    font = ft.font(ft.TIMES_BOLD, 28)
    tmp = Image.new("RGBA", (240, 58), (255, 255, 255, 0))
    d = ImageDraw.Draw(tmp, "RGBA")
    tw, _ = ft.text_box(d, text, font)
    d.rounded_rectangle((5, 5, tw + 28, 47), radius=9, fill=(255, 255, 255, 248), outline=(*color, 210), width=2)
    d.text((17, 11), text, font=font, fill=color)
    tmp = tmp.crop((0, 0, tw + 36, 58)).rotate(-angle, expand=True, resample=Image.Resampling.BICUBIC)
    base.alpha_composite(tmp, (int(center[0] - tmp.width / 2), int(center[1] - tmp.height / 2)))


def draw_leakage_map_clean(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box

    header_font = ft.font(ft.TIMES_BOLD, 43)
    definition_font = ft.font(ft.TIMES, 32)
    header_x, header_y = x0 + 8, y0 + 4
    draw.text((header_x, header_y), "Leakage =", font=header_font, fill=ft.BLUE)
    header_w = ft.text_box(draw, "Leakage =", header_font)[0]
    draw.text((header_x + header_w + 10, header_y + 8), "truth-matched signal photons in PYTHIA MC", font=definition_font, fill=ft.INK)
    draw.text((header_x + header_w + 10, header_y + 47), "found in B/C/D control regions", font=definition_font, fill=ft.INK)

    map_x0, map_y0 = x0 + 8, y0 + 96
    available_w = x1 - x0 - 16
    cell_w = min(360, max(318, int((available_w - 132) / 2)))
    cell_h = 116
    h_gap = available_w - 2 * cell_w
    note_top = y1 - 76
    v_gap = max(62, note_top - map_y0 - 2 * cell_h - 10)
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
        draw.text((cx + 18, cy + 13), region, font=ft.font(ft.TIMES_BOLD, 46), fill=color)
        ft.draw_wrapped(draw, label, (cx + 82, cy + 15), cell_w - 98, ft.font(ft.TIMES_BOLD, 27), fill=ft.INK, line_gap=2)
        ft.draw_wrapped(draw, detail, (cx + 82, cy + 68), cell_w - 98, ft.font(ft.TIMES, 24), fill=ft.MUTED, line_gap=2)

    ax, ay = cells["A"]
    bx, by = cells["B"]
    cx, cy = cells["C"]
    dx, dy = cells["D"]
    c_start = (ax + cell_w // 2, ay - 10)
    c_end = (cx + cell_w // 2, cy + cell_h + 10)
    b_start = (ax + cell_w + 8, ay + cell_h // 2 + 4)
    b_end = (bx - 12, by + cell_h // 2 + 4)
    d_start = (ax + cell_w + 14, ay + 24)
    d_end = (dx + 4, dy + cell_h - 4)
    ft.draw_arrow(draw, c_start, c_end, fill=ft.TEAL, width=5)
    ft.draw_arrow(draw, b_start, b_end, fill=ft.SPHENIX_BLUE, width=5)
    ft.draw_arrow(draw, d_start, d_end, fill=ft.MUTED, width=5)
    b_angle = math.degrees(math.atan2(b_end[1] - b_start[1], b_end[0] - b_start[0]))
    d_angle = math.degrees(math.atan2(d_end[1] - d_start[1], d_end[0] - d_start[0]))
    d_t = 0.62
    d_label = (int(d_start[0] * (1 - d_t) + d_end[0] * d_t), int(d_start[1] * (1 - d_t) + d_end[1] * d_t) - 10)
    draw_leakage_label_clean(base, (c_start[0] + 28, (c_start[1] + c_end[1]) // 2), "C leakage", ft.TEAL)
    draw_leakage_label_clean(base, ((b_start[0] + b_end[0]) // 2, b_start[1] - 18), "B leakage", ft.SPHENIX_BLUE, b_angle)
    draw_leakage_label_clean(base, d_label, "D leakage", ft.MUTED, d_angle)

    note = (x0 + 8, y1 - 84, x1 - 8, y1 - 12)
    icon_x, icon_y = note[0] + 124, note[1] + 17
    draw.rounded_rectangle((icon_x, icon_y, icon_x + 38, icon_y + 38), radius=5, fill=(255, 255, 255, 255), outline=(*ft.SPHENIX_BLUE, 255), width=2)
    for j, w in enumerate([20, 18, 13]):
        draw.line((icon_x + 9, icon_y + 12 + j * 9, icon_x + 9 + w, icon_y + 12 + j * 9), fill=(*ft.SPHENIX_BLUE, 210 - 25 * j), width=3)
    note_font = ft.font(ft.TIMES_BOLD, 29)
    body_font = ft.font(ft.TIMES, 27)
    line1 = "Note: ABCD purity requires isolation-tight ID complementarity;"
    line2 = "isolation is not used in BDT training."
    text_left = icon_x + 48
    text_right = note[2] - 18
    line1_w = ft.text_box(draw, line1, note_font)[0]
    line2_w = ft.text_box(draw, line2, body_font)[0]
    draw.text((text_left + max(0, (text_right - text_left - line1_w) // 2), note[1] + 7), line1, font=note_font, fill=ft.BLUE)
    draw.text((text_left + max(0, (text_right - text_left - line2_w) // 2), note[1] + 40), line2, font=body_font, fill=ft.BLUE)


def draw_derivation_panel(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    draw_card(
        base,
        box,
        ft.SPHENIX_BLUE,
        "Derivation used for the leakage-corrected purity",
        "First clean the B/C/D sidebands; then use the ABCD estimate.",
    )
    x0, y0, x1, y1 = box
    content = (x0 + 58, y0 + 178, x1 - 42, y1 - 34)
    gap = 18
    row_h = [210, 282, 130]
    colors = [ft.PHOTON_DARK, ft.TEAL, ft.SPHENIX_BLUE]
    titles = ["Raw ABCD estimate", "Leakage-correct the sidebands", "Purity carried forward"]
    subtitles = [
        "If B/C/D were background-only controls:",
        "Replace each raw sideband by a leakage-cleaned count:",
        "",
    ]
    y = content[1]
    rows: list[tuple[int, int, int, int]] = []
    for i, h in enumerate(row_h):
        row = (content[0], y, content[2], y + h)
        rows.append(row)
        tint = (255, 250, 236, 255) if i == 0 else (239, 249, 250, 255) if i == 1 else (239, 247, 252, 255)
        draw.rounded_rectangle(row, radius=12, fill=tint, outline=(*colors[i], 185), width=2)
        draw_step_badge(draw, (row[0] + 42, row[1] + 45), str(i + 1), colors[i])
        draw.text((row[0] + 88, row[1] + 20), titles[i], font=ft.font(ft.TIMES_BOLD, 37), fill=ft.INK)
        if subtitles[i]:
            draw.text((row[0] + 88, row[1] + 67), subtitles[i], font=ft.font(ft.TIMES_ITALIC, 31), fill=ft.MUTED)
        draw.rounded_rectangle((row[0], row[1], row[0] + 11, row[3]), radius=6, fill=(*colors[i], 215))
        y += h + gap

    raw_eq_box = (rows[0][0] + 98, rows[0][1] + 104, rows[0][2] - 28, rows[0][3] - 14)
    lc_eq_box = (rows[1][0] + 98, rows[1][1] + 104, rows[1][2] - 28, rows[1][3] - 20)
    pur_eq_box = (rows[2][0] + 104, rows[2][1] + 58, rows[2][2] - 242, rows[2][3] - 18)

    draw_labeled_equation_rows(
        base,
        raw_eq_box,
        [
            ("A background", r"$N_{\mathrm{bkg}}^{A,\mathrm{ABCD}}=N_{\mathrm{raw}}^B\left(N_{\mathrm{raw}}^C/N_{\mathrm{raw}}^D\right)$"),
            ("raw signal", r"$N_{\mathrm{signal}}^{A,\mathrm{raw}}=N_{\mathrm{raw}}^A-N_{\mathrm{bkg}}^{A,\mathrm{ABCD}}$"),
        ],
        colors[0],
        fontsize=35,
    )
    draw_labeled_equation_rows(
        base,
        lc_eq_box,
        [
            ("leakage-corrected", r"$b_X=N_{\mathrm{raw}}^X-f_X^{\mathrm{MC}}\,N_{\mathrm{signal}}^A,\quad X\in\{B,C,D\}$"),
            ("corrected signal", r"$N_{\mathrm{signal}}^{A,\mathrm{LC}}=N_{\mathrm{raw}}^A-b_B\left(b_C/b_D\right)$"),
        ],
        colors[1],
        fontsize=36,
        label_w=296,
        label_size=26,
    )
    draw_labeled_equation_rows(
        base,
        pur_eq_box,
        [("purity", r"$P_{\mathrm{LC}}(E_T)=N_{\mathrm{signal}}^{A,\mathrm{LC}}/N_{\mathrm{raw}}^A$")],
        colors[2],
        fontsize=42,
    )
    pill = (rows[2][2] - 224, rows[2][1] + 64, rows[2][2] - 28, rows[2][3] - 22)
    draw.rounded_rectangle(pill, radius=12, fill=(255, 255, 255, 235), outline=(*ft.SPHENIX_BLUE, 160), width=2)
    text_center_y(draw, pill, "used on Region A", ft.font(ft.TIMES_BOLD, 26), ft.BLUE)


def draw_meaning_panel(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw_card(
        base,
        box,
        ft.PHOTON_DARK,
        "What the correction is doing",
        "A is selected; B/C/D are controls after true-signal leakage is subtracted.",
    )
    inner = (box[0] + 50, box[1] + 178, box[2] - 34, box[3] - 30)
    draw_leakage_map_clean(base, inner)


def draw_takeaway(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    ft.shadow(base, box, radius=10)
    draw.rounded_rectangle(box, radius=14, fill=(244, 249, 250, 255), outline=(184, 213, 217, 255), width=2)
    lead = "Key point:"
    body = "Region A stays fixed; only the B/C/D controls are cleaned before purity is computed."
    lead_font = ft.font(ft.TIMES_BOLD, 35)
    body_font = ft.font(ft.TIMES, 35)
    lead_w = ft.text_box(draw, lead + " ", lead_font)[0]
    body_w = ft.text_box(draw, body, body_font)[0]
    x = box[0] + (box[2] - box[0] - lead_w - body_w) / 2
    y = box[1] + 24
    draw.text((x, y), lead + " ", font=lead_font, fill=ft.BLUE)
    draw.text((x + lead_w, y), body, font=body_font, fill=ft.INK)


def render() -> Image.Image:
    img = draw_slide_shell()
    draw_meaning_panel(img, (132, 318, 1162, 1198))
    draw_derivation_panel(img, (1192, 318, ft.W - 132, 1198))
    draw_takeaway(img, (190, 1218, ft.W - 190, 1302))
    ft.draw_hp2026_identity_footer(img)
    return img


def write_script() -> None:
    SCRIPT_MD.write_text(
        """# Backup Slide Script - ABCD Leakage Correction

This backup slide is for the leakage-correction detail behind the purity extraction.

The left panel defines leakage physically. Region A is the selected tight-and-isolated sample. Regions B, C, and D are the sideband controls. A true prompt photon can still end up in B, C, or D if it fails isolation, fails tight ID, or fails both axes. Those truth-matched signal photons are what we call signal leakage.

The right panel shows the bookkeeping. In the raw ABCD estimate, B, C, and D are treated as background-only controls, so the background in A is estimated as B times C over D. The correction is to replace each sideband count with a leakage-corrected sideband, b_X, where the MC-predicted true-signal leakage is subtracted. Because that leakage term is proportional to the unknown signal count in A, the equation is solved self-consistently.

After that, the leakage-corrected signal count in A divided by the raw Region A count gives the leakage-corrected purity. That is the purity curve carried into the corrected yield.
""",
        encoding="utf-8",
    )


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    img = render()
    img.convert("RGB").save(PNG, "PNG")
    HEADER_JSON.write_text(json.dumps(HEADER_SPEC, indent=2) + "\n", encoding="utf-8")
    write_script()
    MANIFEST.write_text(
        json.dumps(
            {
                "generated_at": datetime.now().isoformat(timespec="seconds"),
                "script": str(Path(__file__).resolve()),
                "output_png": str(PNG),
                "speaker_script": str(SCRIPT_MD),
                "source_reused": {
                    "leakage_map": "make_hp2026_slide16_yield_flow_redesign.draw_leakage_map",
                    "math_rendering": "make_hp2026_slide16_yield_flow_redesign.math_block_image",
                    "deck_chrome": "make_hp2026_fulltalk_candidates",
                },
                "deck_mutation": "none; local backup-slide PNG candidate only",
                "physics_content": "ABCD raw estimate, MC signal-leakage subtraction in B/C/D, self-consistent leakage-corrected signal, and purity P_LC(E_T).",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(PNG)


if __name__ == "__main__":
    main()

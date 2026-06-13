#!/usr/bin/env python3
"""Render HP2026 full-talk full-slide PNG candidates.

Public-facing PPG12 plot candidates must come from the current PPG12 paper PDF
or another explicitly public/approved source. IAN-derived figures and screenshot
crops are recorded as placeholder-only assets until separately approved.
Generated drawings are only explanatory scaffolding and are recorded as such in
the manifest. Google Slides insertion is a separate approval step.
"""

from __future__ import annotations

import json
import math
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont
from pypdf import PdfReader, PdfWriter


W, H = 2560, 1440

ROOT = next(
    p
    for p in Path(__file__).resolve().parents
    if (p / "README.md").exists() and (p / "scripts").exists() and (p / "src").exists()
)
WORKSPACE = ROOT / "outputs/manual-20260601-hp2026-fulltalk/presentations/hp2026-fulltalk"
OUTPUT = WORKSPACE / "output"
ASSETS = WORKSPACE / "assets"
PAGE_ASSETS = ASSETS / "paper_pages"
FIGURE_ASSETS = ASSETS / "paper_figures"
IAN_PAGE_ASSETS = ASSETS / "ian_pages"
IAN_FIGURE_ASSETS = ASSETS / "ian_figures"
SCRIPT_DIR = OUTPUT / "speaker_scripts"
CONTACT_SHEET = OUTPUT / "hp2026_slides05_14_contact_sheet.png"
TITLE_ASSET_DIR = ROOT / "outputs/manual-20260601-hp2026-title/presentations/hp2026-title-slide/assets"
USER_HEADER_SHOWER_SCREENSHOT = Path("/Users/patsfan753/Desktop/Screenshot 2026-06-02 at 2.39.03 PM.png")
SLIDE08_E11E33_FLOW_SCREENSHOT = Path(
    "/var/folders/l3/f02nw86n5cn0tpf_zstf0ypr0000gn/T/TemporaryItems/NSIRD_screencaptureui_jOuoGO/Screenshot 2026-06-02 at 5.59.58 PM.png"
)

PAPER = ROOT / "usefulDocs/sPHENIX_PPG12_Paper_2026-05-21_current_draft.pdf"
IAN = ROOT / "usefulDocs/PPG12_analysis_note_2026-05-21_v4_current_IAN.pdf"

FONT_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES = FONT_DIR / "Times New Roman.ttf"
TIMES_BOLD = FONT_DIR / "Times New Roman Bold.ttf"
TIMES_ITALIC = FONT_DIR / "Times New Roman Italic.ttf"

INK = (18, 22, 28)
BLUE = (19, 41, 75)
SPHENIX_BLUE = (30, 143, 214)
TEAL = (18, 112, 142)
TEAL_SOFT = (102, 173, 190)
PHOTON = (245, 181, 34)
PHOTON_DARK = (196, 124, 10)
MUTED = (71, 79, 91)
LIGHT_MUTED = (117, 126, 138)
SOFT_BG = (252, 253, 254)
PANEL = (246, 249, 252)
PANEL_EDGE = (218, 226, 235)
CARD_EDGE = (221, 228, 236)
HP2026_FOOTER_RULE = (217, 225, 233)
HP2026_FOOTER_RULE_RGBA = (*HP2026_FOOTER_RULE, 255)
HP2026_MAIN_HEADER = {
    "deck": "hp2026_main_talk",
    "title_font_size": 86,
    "subtitle_font_size": None,
    "title_xy": [132, 76],
    "subtitle_xy": None,
    "divider_y": 232,
}


@dataclass(frozen=True)
class FigureSpec:
    key: str
    label: str
    page: int
    box: tuple[int, int, int, int]
    source: str = "paper plot"


@dataclass(frozen=True)
class IanFigureSpec:
    key: str
    label: str
    page: int
    box: tuple[int, int, int, int]
    source: str = "IAN plot"


FIGURES: dict[str, FigureSpec] = {
    "fig1_shower_shape": FigureSpec(
        "fig1_shower_shape",
        "PPG12 paper Fig. 1: shower-shape distributions",
        6,
        (336, 280, 1984, 1096),
    ),
    "fig2_bdt_score": FigureSpec(
        "fig2_bdt_score",
        "PPG12 paper Fig. 2: photon-ID BDT score",
        7,
        (640, 280, 1680, 1272),
    ),
    "fig3_isolation": FigureSpec(
        "fig3_isolation",
        "PPG12 paper Fig. 3: reconstructed isolation energy",
        9,
        (728, 280, 1588, 1044),
    ),
    "fig4_abcd": FigureSpec(
        "fig4_abcd",
        "PPG12 paper Fig. 4: ABCD purity regions",
        9,
        (728, 1172, 1588, 1800),
    ),
    "fig5_purity": FigureSpec(
        "fig5_purity",
        "PPG12 paper Fig. 5: purity vs photon ET",
        10,
        (728, 280, 1588, 1104),
    ),
    "fig6_efficiencies": FigureSpec(
        "fig6_efficiencies",
        "PPG12 paper Fig. 6: reconstruction, ID, isolation, total efficiencies",
        12,
        (728, 280, 1588, 1104),
    ),
    "fig7_systematics": FigureSpec(
        "fig7_systematics",
        "PPG12 paper Fig. 7: systematic uncertainty breakdown",
        14,
        (416, 280, 1912, 1240),
    ),
    "fig8_cross_section": FigureSpec(
        "fig8_cross_section",
        "PPG12 paper Fig. 8: isolated prompt-photon cross section",
        15,
        (456, 320, 1860, 2196),
    ),
    "fig9_phenix": FigureSpec(
        "fig9_phenix",
        "PPG12 paper Fig. 9: comparison with PHENIX",
        16,
        (456, 560, 1860, 2080),
    ),
}


IAN_FIGURES: dict[str, IanFigureSpec] = {
    "fig15_npb_score_vs_time": IanFigureSpec(
        "fig15_npb_score_vs_time",
        "PPG12 IAN Fig. 15: NPB BDT score vs cluster-MBD time",
        22,
        (705, 1230, 1620, 1988),
    ),
    "fig17_npb_threshold_scan": IanFigureSpec(
        "fig17_npb_threshold_scan",
        "PPG12 IAN Fig. 17 left: NPB purity and signal-retention threshold scan",
        24,
        (420, 318, 1148, 1036),
    ),
    "fig58_preselection_e32e35_lowpt_main": IanFigureSpec(
        "fig58_preselection_e32e35_lowpt_main",
        "PPG12 IAN Fig. 58: E3x2/E3x5 main panel after preselection, 10 < ET < 14 GeV",
        68,
        (875, 440, 1435, 895),
    ),
    "fig20_tight_e32e35_lowpt_main": IanFigureSpec(
        "fig20_tight_e32e35_lowpt_main",
        "PPG12 IAN Fig. 20: E3x2/E3x5 main panel after tight photon-ID BDT selection, 10 < ET < 14 GeV",
        27,
        (875, 440, 1435, 895),
    ),
}


def font(path: Path, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(path), size)


def text_box(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.ImageFont) -> tuple[int, int]:
    box = draw.textbbox((0, 0), text, font=fnt)
    return box[2] - box[0], box[3] - box[1]


def crop_visible(img: Image.Image, white_threshold: int = 248, pad: int = 0) -> Image.Image:
    rgba = img.convert("RGBA")
    pix = rgba.load()
    xs: list[int] = []
    ys: list[int] = []
    for y in range(rgba.height):
        for x in range(rgba.width):
            r, g, b, a = pix[x, y]
            if a > 12 and not (r >= white_threshold and g >= white_threshold and b >= white_threshold):
                xs.append(x)
                ys.append(y)
    if not xs:
        return rgba
    x0 = max(0, min(xs) - pad)
    y0 = max(0, min(ys) - pad)
    x1 = min(rgba.width, max(xs) + 1 + pad)
    y1 = min(rgba.height, max(ys) + 1 + pad)
    return rgba.crop((x0, y0, x1, y1))


def load_sphenix_logo() -> Image.Image | None:
    logo = TITLE_ASSET_DIR / "sphenix-logo-white-bg_0.png"
    if not logo.exists():
        return None
    return crop_visible(Image.open(logo).convert("RGBA"), white_threshold=252)


def white_to_transparent(img: Image.Image, threshold: int = 248) -> Image.Image:
    rgba = img.convert("RGBA")
    out = Image.new("RGBA", rgba.size, (0, 0, 0, 0))
    src = rgba.load()
    dst = out.load()
    for y in range(rgba.height):
        for x in range(rgba.width):
            r, g, b, a = src[x, y]
            if a == 0:
                continue
            whiteness = min(r, g, b)
            if whiteness >= threshold:
                continue
            alpha = a
            if whiteness > threshold - 18:
                alpha = int(a * (threshold - whiteness) / 18)
            dst[x, y] = (r, g, b, alpha)
    return out


def add_top_right_sphenix_logo(base: Image.Image) -> None:
    logo = load_sphenix_logo()
    if logo is not None:
        paste_fit(base, logo, (2128, 58, 2444, 150), anchor="right")


def add_top_right_sphenix_logo_like_slide2(base: Image.Image) -> None:
    logo = load_sphenix_logo()
    if logo is not None:
        paste_fit(base, logo, (2188, 58, 2432, 164), anchor="right")


def draw_hp2026_identity_footer(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    footer_top = 1326
    draw.rectangle((0, footer_top, W, H), fill=(*SOFT_BG, 255))
    draw.line((0, footer_top, W, footer_top), fill=HP2026_FOOTER_RULE_RGBA, width=2)

    illinois_path = TITLE_ASSET_DIR / "illinois_logo_fullcolor_rgb.png"
    hp_path = TITLE_ASSET_DIR / "hp2026_indico_logo.png"
    cy = 1384
    if illinois_path.exists():
        illinois = crop_visible(Image.open(illinois_path).convert("RGBA"), white_threshold=252)
        illinois = fit(illinois, 54, 62)
        base.alpha_composite(illinois, (30, cy - illinois.height // 2))
    draw.text((104, cy - 17), "Justin Bennett", font=font(TIMES, 31), fill=(43, 49, 57))

    center = "Hard Probes 2026 / June 24, 2026"
    center_font = font(TIMES_BOLD, 31)
    cw, ch = text_box(draw, center, center_font)
    hp = None
    if hp_path.exists():
        hp = white_to_transparent(Image.open(hp_path).convert("RGBA"), threshold=246)
        hp = fit(hp, 116, 60)
    group_w = (hp.width + 20 if hp is not None else 0) + cw
    gx = (W - group_w) // 2
    if hp is not None:
        base.alpha_composite(hp, (gx, cy - hp.height // 2))
        gx += hp.width + 20
    draw.text((gx, cy - ch // 2 - 1), center, font=center_font, fill=(43, 49, 57))


def draw_header_shower_icon(base: Image.Image) -> None:
    def draw_gamma_label(draw: ImageDraw.ImageDraw, center: tuple[int, int], size: int = 34) -> None:
        color = (7, 101, 151)
        label_font = font(TIMES_ITALIC, size)
        tw, th = text_box(draw, "γ", label_font)
        draw.text(
            (center[0] - tw / 2, center[1] - th / 2 - 2),
            "γ",
            font=label_font,
            fill=color,
            stroke_width=max(2, round(size * 0.10)),
            stroke_fill=(255, 255, 255, 248),
        )

    def draw_pi0_eta_label(draw: ImageDraw.ImageDraw, center: tuple[int, int], size: int = 28) -> None:
        color = (202, 112, 28)
        main = font(TIMES_ITALIC, 24)
        sup = font(TIMES_ITALIC, 15)
        if size != 28:
            main = font(TIMES_ITALIC, size)
            sup = font(TIMES_ITALIC, max(13, round(size * 0.62)))
        parts = [("π", main, 0), ("0", sup, -round(size * 0.36)), ("/η", main, 0)]
        widths = [text_box(draw, text, fnt)[0] for text, fnt, _ in parts]
        heights = [text_box(draw, text, fnt)[1] for text, fnt, _ in parts]
        cursor = center[0] - sum(widths) / 2
        y = center[1] - max(heights) / 2
        for (text, fnt, dy), width in zip(parts, widths):
            draw.text(
                (cursor, y + dy),
                text,
                font=fnt,
                fill=color,
                stroke_width=max(2, round(size * 0.08)),
                stroke_fill=(255, 255, 255, 248),
            )
            cursor += width

    scale = 4
    icon_w, icon_h = 392, 184
    layer = Image.new("RGBA", (icon_w * scale, icon_h * scale), (0, 0, 0, 0))
    draw = ImageDraw.Draw(layer, "RGBA")

    def spts(points: list[tuple[float, float]]) -> list[tuple[int, int]]:
        return [(round(x * scale), round(y * scale)) for x, y in points]

    def seg(points: list[tuple[float, float]], fill=(237, 241, 242, 238), outline=(126, 134, 140, 225)) -> None:
        draw.polygon(spts(points), fill=fill)
        draw.line(spts(points + [points[0]]), fill=outline, width=round(2.2 * scale), joint="curve")

    def jagged_blob(
        cx: float,
        cy: float,
        rx: float,
        ry: float,
        angle_deg: float,
        fill: tuple[int, int, int, int],
        n: int = 18,
    ) -> None:
        pts = []
        angle = math.radians(angle_deg)
        for i in range(n):
            theta = angle + i * math.tau / n
            jitter = 0.62 + 0.44 * (((i * 5 + 3) % 9) / 8)
            pts.append((cx + math.cos(theta) * rx * jitter, cy + math.sin(theta) * ry * jitter))
        draw.polygon(spts(pts), fill=fill)

    # Calorimeter modules, hand-drawn as two gently curved stacks to echo the
    # standard shower-shape explanatory cartoon without pasting a raster source.
    left_rows = [
        [(24, 54), (174, 42), (174, 72), (34, 84)],
        [(34, 88), (174, 76), (173, 106), (44, 116)],
        [(44, 120), (173, 110), (172, 140), (53, 148)],
        [(53, 152), (172, 144), (171, 174), (62, 181)],
    ]
    right_rows = [
        [(218, 42), (368, 54), (358, 84), (218, 72)],
        [(219, 76), (358, 88), (348, 116), (219, 106)],
        [(220, 110), (348, 120), (338, 148), (220, 140)],
        [(221, 144), (338, 152), (328, 181), (221, 174)],
    ]
    for row in left_rows + right_rows:
        seg(row)

    # Compact photon-like shower.
    blue = [(101, 58), (121, 176), (85, 179)]
    draw.polygon(spts(blue), fill=(26, 152, 211, 238))
    draw.line(spts(blue + [blue[0]]), fill=(7, 101, 151, 240), width=round(1.8 * scale), joint="curve")
    wave = feynman_points((105 * scale, 36 * scale), (104 * scale, 176 * scale), 5.5 * scale, 5.6, 140)
    draw_polyline(draw, wave, (255, 255, 255, 245), round(5.0 * scale))
    draw_polyline(draw, wave, (7, 101, 151, 255), round(2.6 * scale))
    draw_polyline(draw, wave, (83, 198, 235, 255), round(1.2 * scale))
    draw_gamma_label(draw, (round(105 * scale), round(18 * scale)), size=35 * scale)

    # Diffuse split/decay-like shower: a central orange deposit with gray
    # shoulders and small splashes, all clipped visually by the module stack.
    gray_blobs = [
        (286, 118, 56, 42, -6),
        (272, 94, 38, 20, -18),
        (309, 102, 32, 17, 12),
        (256, 124, 38, 17, 22),
        (315, 139, 36, 17, 15),
        (280, 158, 42, 17, 6),
    ]
    for blob in gray_blobs:
        jagged_blob(*blob, fill=(105, 111, 114, 128), n=21)

    orange_core = [
        (278, 57), (293, 74), (294, 112), (306, 150), (296, 179),
        (268, 177), (271, 141), (260, 112), (267, 83),
    ]
    draw.polygon(spts(orange_core), fill=(244, 137, 32, 226))
    jagged_blob(282, 116, 26, 57, -2, fill=(244, 137, 32, 218), n=25)
    jagged_blob(288, 148, 23, 34, 8, fill=(244, 137, 32, 215), n=19)
    hadron_path = spts([(286, 36), (296, 68), (276, 93), (305, 119), (287, 146), (307, 178)])
    draw.line(hadron_path, fill=(255, 250, 238, 240), width=round(6.8 * scale), joint="curve")
    draw.line(hadron_path, fill=(202, 112, 28, 255), width=round(3.2 * scale), joint="curve")
    draw.line(hadron_path, fill=(255, 174, 54, 245), width=round(1.4 * scale), joint="curve")
    branches = [
        [(276, 93), (249, 76)],
        [(276, 93), (244, 109)],
        [(305, 119), (336, 101)],
        [(287, 146), (330, 160)],
    ]
    for branch in branches:
        pts = spts(branch)
        draw.line(pts, fill=(255, 250, 238, 235), width=round(3.8 * scale))
        draw.line(pts, fill=(202, 112, 28, 220), width=round(1.8 * scale))
    draw_pi0_eta_label(draw, (round(286 * scale), round(36 * scale)), size=31 * scale)

    orange_marks = [
        (245, 67, 5, 16, -25), (326, 78, 5, 16, 20), (341, 121, 5, 14, 68),
        (241, 116, 5, 15, -65), (317, 160, 5, 16, -35), (260, 167, 5, 14, 42),
        (332, 151, 4, 12, 24), (252, 97, 4, 11, -42), (294, 67, 3, 10, 8),
    ]
    for cx, cy, rx, ry, angle in orange_marks:
        jagged_blob(cx, cy, rx, ry, angle, fill=(244, 137, 32, 210), n=9)

    layer = layer.resize((392, 184), Image.Resampling.LANCZOS)
    base.alpha_composite(layer, (1808, 62))


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    text: str,
    xy: tuple[int, int],
    max_width: int,
    fnt: ImageFont.ImageFont,
    fill: tuple[int, int, int] = INK,
    line_gap: int = 8,
) -> int:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        trial = word if not current else f"{current} {word}"
        if text_box(draw, trial, fnt)[0] <= max_width:
            current = trial
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)

    x, y = xy
    for line in lines:
        draw.text((x, y), line, font=fnt, fill=fill)
        y += text_box(draw, line, fnt)[1] + line_gap
    return y


def feynman_points(
    start: tuple[float, float],
    end: tuple[float, float],
    amplitude: float,
    cycles: float,
    steps: int,
) -> list[tuple[float, float]]:
    sx, sy = start
    ex, ey = end
    dx, dy = ex - sx, ey - sy
    length = math.hypot(dx, dy)
    if length == 0:
        return [start]
    nx, ny = -dy / length, dx / length
    pts = []
    for i in range(steps + 1):
        t = i / steps
        wave = math.sin(t * cycles * math.tau) * amplitude
        pts.append((sx + dx * t + nx * wave, sy + dy * t + ny * wave))
    return pts


def draw_polyline(draw: ImageDraw.ImageDraw, pts: list[tuple[float, float]], fill, width: int) -> None:
    draw.line([(round(x), round(y)) for x, y in pts], fill=fill, width=width, joint="curve")


def fit(img: Image.Image, max_w: int, max_h: int) -> Image.Image:
    scale = min(max_w / img.width, max_h / img.height)
    size = (max(1, int(img.width * scale)), max(1, int(img.height * scale)))
    return img.resize(size, Image.Resampling.LANCZOS)


def paste_fit(
    base: Image.Image,
    img: Image.Image,
    box: tuple[int, int, int, int],
    anchor: str = "center",
) -> tuple[int, int, int, int]:
    x0, y0, x1, y1 = box
    fitted = fit(img.convert("RGBA"), x1 - x0, y1 - y0)
    if anchor == "left":
        x = x0
    elif anchor == "right":
        x = x1 - fitted.width
    else:
        x = x0 + ((x1 - x0) - fitted.width) // 2
    y = y0 + ((y1 - y0) - fitted.height) // 2
    base.alpha_composite(fitted, (x, y))
    return (x, y, x + fitted.width, y + fitted.height)


def header(draw: ImageDraw.ImageDraw, title: str, subtitle: str | None = None) -> None:
    draw.rectangle((0, 0, W, 22), fill=SPHENIX_BLUE)
    draw.rectangle((0, 22, W, 30), fill=PHOTON)
    title_font = font(TIMES_BOLD, 74)
    draw.text((132, 84), title, font=title_font, fill=INK)
    if subtitle:
        draw.text((136, 184), subtitle, font=font(TIMES_ITALIC, 36), fill=BLUE)
    draw.line((132, 278, W - 132, 278), fill=(221, 226, 232), width=3)


def hp2026_main_header(draw: ImageDraw.ImageDraw, title: str, subtitle: str | None = None) -> None:
    draw.rectangle((0, 0, W, 22), fill=SPHENIX_BLUE)
    draw.rectangle((0, 22, W, 30), fill=PHOTON)
    draw.text(
        tuple(HP2026_MAIN_HEADER["title_xy"]),
        title,
        font=font(TIMES_BOLD, HP2026_MAIN_HEADER["title_font_size"]),
        fill=INK,
    )
    if subtitle and HP2026_MAIN_HEADER.get("subtitle_xy") and HP2026_MAIN_HEADER.get("subtitle_font_size"):
        draw.text(
            tuple(HP2026_MAIN_HEADER["subtitle_xy"]),
            subtitle,
            font=font(TIMES_ITALIC, HP2026_MAIN_HEADER["subtitle_font_size"]),
            fill=MUTED,
        )
    y = HP2026_MAIN_HEADER["divider_y"]
    draw.line((132, y, W - 132, y), fill=(221, 226, 232), width=3)


def base_slide(title: str, subtitle: str | None = None) -> Image.Image:
    img = Image.new("RGBA", (W, H), (*SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    header(draw, title, subtitle)
    return img


def base_slide_hp2026_main(title: str, subtitle: str | None = None) -> Image.Image:
    img = Image.new("RGBA", (W, H), (*SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    hp2026_main_header(draw, title, subtitle)
    return img


def write_hp2026_main_header_spec(png: Path) -> None:
    png.with_suffix(".header.json").write_text(
        json.dumps({"hp2026_main_header": HP2026_MAIN_HEADER}, indent=2) + "\n",
        encoding="utf-8",
    )


def rounded_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], fill=(255, 255, 255), radius=10) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=(*fill, 255), outline=(*PANEL_EDGE, 255), width=2)


def outer_card_sidebar(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], accent: tuple[int, int, int], *, width: int = 12) -> None:
    draw.rounded_rectangle((box[0], box[1], box[0] + width, box[3]), radius=6, fill=(*accent, 235))


def shadow(base: Image.Image, box: tuple[int, int, int, int], radius: int = 10) -> None:
    layer = Image.new("RGBA", base.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(layer, "RGBA")
    d.rounded_rectangle((box[0] + 8, box[1] + 10, box[2] + 8, box[3] + 10), radius=radius, fill=(30, 42, 58, 30))
    layer = layer.filter(ImageFilter.GaussianBlur(12))
    base.alpha_composite(layer)


def figure_path(key: str) -> Path:
    return FIGURE_ASSETS / f"{key}.png"


def ian_figure_path(key: str) -> Path:
    return IAN_FIGURE_ASSETS / f"{key}.png"


def render_pdf_page(page_no: int) -> Path:
    PAGE_ASSETS.mkdir(parents=True, exist_ok=True)
    pdf_out = PAGE_ASSETS / f"paper_page_{page_no:02d}.pdf"
    png_out = PAGE_ASSETS / f"paper_page_{page_no:02d}.png"
    if png_out.exists():
        return png_out

    reader = PdfReader(str(PAPER))
    writer = PdfWriter()
    writer.add_page(reader.pages[page_no - 1])
    with pdf_out.open("wb") as f:
        writer.write(f)
    subprocess.run(
        ["/usr/bin/sips", "-Z", "3000", "-s", "format", "png", str(pdf_out), "--out", str(png_out)],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    return png_out


def render_ian_pdf_page(page_no: int) -> Path:
    IAN_PAGE_ASSETS.mkdir(parents=True, exist_ok=True)
    pdf_out = IAN_PAGE_ASSETS / f"ian_page_{page_no:02d}.pdf"
    png_out = IAN_PAGE_ASSETS / f"ian_page_{page_no:02d}.png"
    if png_out.exists():
        return png_out

    reader = PdfReader(str(IAN))
    writer = PdfWriter()
    writer.add_page(reader.pages[page_no - 1])
    with pdf_out.open("wb") as f:
        writer.write(f)
    subprocess.run(
        ["/usr/bin/sips", "-Z", "3000", "-s", "format", "png", str(pdf_out), "--out", str(png_out)],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    return png_out


def crop_figure(spec: FigureSpec) -> Path:
    FIGURE_ASSETS.mkdir(parents=True, exist_ok=True)
    out = figure_path(spec.key)
    if out.exists():
        return out
    page = Image.open(render_pdf_page(spec.page)).convert("RGBA")
    crop = page.crop(spec.box)
    white = Image.new("RGBA", crop.size, (255, 255, 255, 255))
    white.alpha_composite(crop)
    white.convert("RGB").save(out, "PNG")
    return out


def crop_ian_figure(spec: IanFigureSpec) -> Path:
    IAN_FIGURE_ASSETS.mkdir(parents=True, exist_ok=True)
    out = ian_figure_path(spec.key)
    if out.exists():
        return out
    page = Image.open(render_ian_pdf_page(spec.page)).convert("RGBA")
    crop = page.crop(spec.box)
    white = Image.new("RGBA", crop.size, (255, 255, 255, 255))
    white.alpha_composite(crop)
    white.convert("RGB").save(out, "PNG")
    return out


def prepare_figures() -> dict[str, Path]:
    return {key: crop_figure(spec) for key, spec in FIGURES.items()}


def prepare_ian_figures() -> dict[str, Path]:
    return {key: crop_ian_figure(spec) for key, spec in IAN_FIGURES.items()}


def slide08_e11e33_flow_paths() -> dict[str, Path]:
    asset_dir = ASSETS / "slide08_e11e33_flow"
    asset_dir.mkdir(parents=True, exist_ok=True)
    crops = {
        "no_preselection": {
            "box": (62, 114, 349, 397),
            "label": "Backup Slide 13 crop: Fig. 13 E11/E33 without preselection, 22 < pT < 28 GeV",
        },
        "preselection": {
            "box": (430, 114, 732, 415),
            "label": "Backup Slide 13 crop: Fig. 19 E11/E33 with preselection, 18 < pT < 22 GeV",
        },
        "tight_id": {
            "box": (805, 122, 1126, 443),
            "label": "Backup Slide 13 crop: Fig. 19 E11/E33 with tight selection, 22 < pT < 28 GeV",
        },
    }
    out_paths: dict[str, Path] = {}
    if not SLIDE08_E11E33_FLOW_SCREENSHOT.exists():
        return out_paths
    screenshot = Image.open(SLIDE08_E11E33_FLOW_SCREENSHOT).convert("RGBA")
    for key, spec in crops.items():
        out = asset_dir / f"{key}_e11e33_main_panel.png"
        if not out.exists():
            crop = screenshot.crop(spec["box"])
            crop = crop.resize((crop.width * 2, crop.height * 2), Image.Resampling.LANCZOS)
            crop.save(out, "PNG")
        out_paths[key] = out
    return out_paths


def place_figure(
    base: Image.Image,
    key: str,
    box: tuple[int, int, int, int],
    label: str | None = None,
    anchor: str = "center",
    inset: int = 28,
    accent: tuple[int, int, int] | None = None,
) -> tuple[int, int, int, int]:
    draw = ImageDraw.Draw(base, "RGBA")
    shadow(base, box)
    rounded_panel(draw, box)
    if accent:
        outer_card_sidebar(draw, box, accent)
    img = Image.open(figure_path(key)).convert("RGBA")
    placed = paste_fit(base, img, (box[0] + inset, box[1] + inset, box[2] - inset, box[3] - inset), anchor=anchor)
    if label:
        draw.text((box[0] + 28, box[1] + 18), label, font=font(TIMES_ITALIC, 26), fill=LIGHT_MUTED)
    return placed


def place_figure_tight(
    base: Image.Image,
    key: str,
    box: tuple[int, int, int, int],
    label: str | None = None,
    anchor: str = "center",
    inset: int = 18,
    pad: int = 20,
) -> tuple[int, int, int, int]:
    draw = ImageDraw.Draw(base, "RGBA")
    shadow(base, box)
    rounded_panel(draw, box)
    img = Image.open(figure_path(key)).convert("RGBA")
    img = crop_visible(img, white_threshold=252, pad=pad)
    placed = paste_fit(base, img, (box[0] + inset, box[1] + inset, box[2] - inset, box[3] - inset), anchor=anchor)
    if label:
        draw.text((box[0] + 28, box[1] + 18), label, font=font(TIMES_ITALIC, 26), fill=LIGHT_MUTED)
    return placed


def place_figure_snug_panel(
    base: Image.Image,
    key: str,
    max_box: tuple[int, int, int, int],
    pad: int = 22,
    white_threshold: int = 252,
    crop_pad: int = 18,
    trim_bottom: int = 0,
) -> tuple[tuple[int, int, int, int], tuple[int, int, int, int]]:
    x0, y0, x1, y1 = max_box
    draw = ImageDraw.Draw(base, "RGBA")
    img = Image.open(figure_path(key)).convert("RGBA")
    img = crop_visible(img, white_threshold=white_threshold, pad=crop_pad)
    if trim_bottom > 0 and img.height > trim_bottom + 20:
        img = img.crop((0, 0, img.width, img.height - trim_bottom))
    fitted = fit(img, x1 - x0 - 2 * pad, y1 - y0 - 2 * pad)
    px = x0 + ((x1 - x0) - fitted.width) // 2
    py = y0 + ((y1 - y0) - fitted.height) // 2
    panel = (px - pad, py - pad, px + fitted.width + pad, py + fitted.height + pad)
    shadow(base, panel)
    rounded_panel(draw, panel)
    base.alpha_composite(fitted, (px, py))
    return (px, py, px + fitted.width, py + fitted.height), panel


def place_ian_figure(
    base: Image.Image,
    key: str,
    box: tuple[int, int, int, int],
    label: str | None = None,
    anchor: str = "center",
) -> tuple[int, int, int, int]:
    draw = ImageDraw.Draw(base, "RGBA")
    shadow(base, box)
    rounded_panel(draw, box)
    inset = 28
    img = Image.open(ian_figure_path(key)).convert("RGBA")
    placed = paste_fit(base, img, (box[0] + inset, box[1] + inset, box[2] - inset, box[3] - inset), anchor=anchor)
    if label:
        draw.text((box[0] + 28, box[1] + 18), label, font=font(TIMES_ITALIC, 26), fill=LIGHT_MUTED)
    return placed


def claim_box(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], heading: str, body: str, accent=SPHENIX_BLUE) -> None:
    draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 255), outline=(*CARD_EDGE, 255), width=2)
    draw.rounded_rectangle((box[0], box[1], box[0] + 14, box[3]), radius=6, fill=(*accent, 255))
    draw.text((box[0] + 42, box[1] + 28), heading, font=font(TIMES_BOLD, 37), fill=INK)
    draw_wrapped(draw, body, (box[0] + 42, box[1] + 84), box[2] - box[0] - 78, font(TIMES, 28), fill=MUTED, line_gap=8)


def bridge(draw: ImageDraw.ImageDraw, text: str) -> None:
    draw.rounded_rectangle((132, 1240, W - 132, 1326), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw.text((174, 1264), text, font=font(TIMES_ITALIC, 36), fill=BLUE)


def draw_shower_pair(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    rounded_panel(draw, box, fill=PANEL)
    x0, y0, x1, y1 = box
    draw.text((x0 + 34, y0 + 22), "What the BDT sees", font=font(TIMES_BOLD, 31), fill=BLUE)
    centers = [(x0 + 150, y0 + 150), (x0 + 430, y0 + 150)]
    labels = [("prompt", PHOTON), ("decay", TEAL)]
    for idx, ((cx, cy), (label, color)) in enumerate(zip(centers, labels)):
        draw.text((cx - 64, y0 + 72), label, font=font(TIMES_ITALIC, 28), fill=MUTED)
        for r in range(18, 82, 14):
            alpha = max(20, 130 - r)
            if idx == 0:
                draw.ellipse((cx - r, cy - r, cx + r, cy + r), outline=(*color, alpha), width=4)
            else:
                draw.ellipse((cx - r - 36, cy - r, cx + r - 36, cy + r), outline=(*color, alpha), width=3)
                draw.ellipse((cx - r + 36, cy - r, cx + r + 36, cy + r), outline=(*color, alpha), width=3)
        draw.ellipse((cx - 12, cy - 12, cx + 12, cy + 12), fill=(*color, 210))
    draw.text((x0 + 66, y1 - 62), "single core", font=font(TIMES, 25), fill=MUTED)
    draw.text((x0 + 346, y1 - 62), "merged cores", font=font(TIMES, 25), fill=MUTED)


def draw_variable_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    heading: str,
    formula_parts: list[tuple[str, int, int]],
    body: str,
    accent: tuple[int, int, int],
) -> None:
    draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 255), outline=(*CARD_EDGE, 255), width=2)
    draw.rounded_rectangle((box[0], box[1], box[0] + 12, box[3]), radius=6, fill=(*accent, 255))
    draw.text((box[0] + 34, box[1] + 20), heading, font=font(TIMES_BOLD, 30), fill=INK)
    x = box[0] + 34
    y = box[1] + 66
    for text, size, dy in formula_parts:
        fnt = font(TIMES_ITALIC, size)
        draw.text((x, y + dy), text, font=fnt, fill=BLUE)
        x += text_box(draw, text, fnt)[0]
    draw_wrapped(draw, body, (box[0] + 34, box[1] + 108), box[2] - box[0] - 60, font(TIMES, 24), fill=MUTED, line_gap=6)


def draw_formula_run(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    parts: list[tuple[str, int, int, Path]],
    fill: tuple[int, int, int] = BLUE,
) -> int:
    x, y = xy
    for text, size, dy, font_path in parts:
        fnt = font(font_path, size)
        draw.text((x, y + dy), text, font=fnt, fill=fill)
        x += text_box(draw, text, fnt)[0]
    return x


def callout_label(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    fill: tuple[int, int, int] = BLUE,
    bg: tuple[int, int, int] = (255, 255, 255),
    size: int = 22,
) -> None:
    x, y = xy
    fnt = font(TIMES_BOLD, size)
    tw, th = text_box(draw, text, fnt)
    draw.rounded_rectangle((x - 8, y - 5, x + tw + 8, y + th + 9), radius=5, fill=(*bg, 228))
    draw.text((x, y), text, font=fnt, fill=fill)


def draw_small_tower_grid(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    cell: int,
    values: list[list[float]],
    highlight: tuple[int, int, int, int] | None = None,
    highlight_color: tuple[int, int, int] = PHOTON_DARK,
    secondary: tuple[int, int, int, int] | None = None,
    secondary_color: tuple[int, int, int] = SPHENIX_BLUE,
) -> tuple[int, int, int, int]:
    x0, y0 = xy
    for r, row in enumerate(values):
        for c, v in enumerate(row):
            red = int(253 - 34 * min(v, 1))
            green = int(248 - 136 * min(v, 1))
            blue = int(226 - 214 * min(v, 1))
            draw.rectangle(
                (x0 + c * cell, y0 + r * cell, x0 + (c + 1) * cell, y0 + (r + 1) * cell),
                fill=(red, green, max(28, blue), 255),
                outline=(212, 219, 227, 255),
                width=2,
            )
    if secondary:
        c0, r0, c1, r1 = secondary
        draw.rectangle((x0 + c0 * cell, y0 + r0 * cell, x0 + c1 * cell, y0 + r1 * cell), outline=(*secondary_color, 210), width=4)
    if highlight:
        c0, r0, c1, r1 = highlight
        draw.rectangle((x0 + c0 * cell, y0 + r0 * cell, x0 + c1 * cell, y0 + r1 * cell), outline=(*highlight_color, 240), width=5)
    cx = x0 + (len(values[0]) / 2) * cell
    cy = y0 + (len(values) / 2) * cell
    draw.ellipse((cx - 8, cy - 8, cx + 8, cy + 8), fill=(0, 0, 0, 255))
    return (x0, y0, x0 + len(values[0]) * cell, y0 + len(values) * cell)


def draw_compact_core_panel(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    title = "Is the core compact?"
    tw, _ = text_box(draw, title, font(TIMES_BOLD, 39))
    draw.text((x0 + (x1 - x0 - tw) / 2, y0 + 28), title, font=font(TIMES_BOLD, 39), fill=INK)
    values = [
        [0.03, 0.05, 0.07, 0.05, 0.03],
        [0.05, 0.18, 0.32, 0.18, 0.05],
        [0.07, 0.36, 0.95, 0.42, 0.07],
        [0.05, 0.20, 0.38, 0.19, 0.05],
        [0.03, 0.05, 0.07, 0.05, 0.03],
    ]
    cell = 82
    gx = x0 + (x1 - x0 - 5 * cell) // 2
    gy = y0 + 110
    draw_small_tower_grid(draw, (gx, gy), cell, values, highlight=(1, 1, 3, 3), highlight_color=PHOTON_DARK, secondary=(1, 1, 4, 4), secondary_color=SPHENIX_BLUE)
    cog = (gx + 2.5 * cell, gy + 2.5 * cell)
    draw.line((gx + 2 * cell - 18, gy + cell + 58, gx + cell + 14, gy + cell + 14), fill=(*PHOTON_DARK, 180), width=2)
    draw.line((gx + 3 * cell + 18, gy + 3 * cell + 18, gx + 4 * cell - 10, gy + 4 * cell - 10), fill=(*SPHENIX_BLUE, 170), width=2)
    callout_label(draw, (gx + 2 * cell - 70, gy + cell + 18), "2x2 core", PHOTON_DARK, size=21)
    callout_label(draw, (gx + 3 * cell + 6, gy + 3 * cell + 16), "3x3 local core", SPHENIX_BLUE, size=21)
    callout_label(draw, (round(cog[0] + 20), round(cog[1] - 18)), "CoG", INK, size=21)
    draw.text((gx + 2 * cell + 22, gy + 2 * cell + 50), "seed", font=font(TIMES_ITALIC, 22), fill=INK)
    draw_formula_run(
        draw,
        (x0 + 54, y0 + 548),
        [("et1", 31, 0, TIMES_BOLD), (" = energy fraction in the 2x2 core", 31, 0, TIMES)],
        fill=INK,
    )
    draw_formula_run(
        draw,
        (x0 + 54, y0 + 616),
        [("e11/e33", 31, 0, TIMES_BOLD), (" = center tower over local 3x3 core", 31, 0, TIMES)],
        fill=INK,
    )
    draw_wrapped(draw, "Photon-like: energy stays concentrated in the core.", (x0 + 54, y0 + 718), x1 - x0 - 100, font(TIMES_ITALIC, 31), fill=BLUE, line_gap=8)


def draw_narrow_shoulders_panel(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    title = "Are the shoulders narrow?"
    tw, _ = text_box(draw, title, font(TIMES_BOLD, 39))
    draw.text((x0 + (x1 - x0 - tw) / 2, y0 + 28), title, font=font(TIMES_BOLD, 39), fill=INK)
    compact = [
        [0.01, 0.02, 0.04, 0.02, 0.01],
        [0.02, 0.08, 0.20, 0.08, 0.02],
        [0.04, 0.18, 0.92, 0.18, 0.04],
        [0.02, 0.08, 0.20, 0.08, 0.02],
        [0.01, 0.02, 0.04, 0.02, 0.01],
    ]
    broad = [
        [0.10, 0.18, 0.25, 0.18, 0.10],
        [0.18, 0.34, 0.46, 0.34, 0.18],
        [0.25, 0.48, 0.80, 0.48, 0.25],
        [0.18, 0.34, 0.46, 0.34, 0.18],
        [0.10, 0.18, 0.25, 0.18, 0.10],
    ]
    cell = 60
    left = draw_small_tower_grid(draw, (x0 + 48, y0 + 124), cell, compact)
    right = draw_small_tower_grid(draw, (x0 + 392, y0 + 124), cell, broad)
    for idx, grid in enumerate((left, right)):
        gx0, gy0, gx1, gy1 = grid
        if idx == 0:
            start, end = gx0 + 112, gx1 - 112
        else:
            start, end = gx0 + 24, gx1 - 24
        draw.line((start, gy0 - 26, end, gy0 - 26), fill=(230, 70, 45, 230), width=4)
        draw.line((start, gy0 - 36, start, gy0 - 16), fill=(230, 70, 45, 230), width=4)
        draw.line((end, gy0 - 36, end, gy0 - 16), fill=(230, 70, 45, 230), width=4)
        center_cell = (gx0 + 2 * cell, gy0 + 2 * cell, gx0 + 3 * cell, gy0 + 3 * cell)
        draw.rectangle(center_cell, fill=(255, 255, 255, 210), outline=(160, 168, 176, 210), width=2)
        draw.line((center_cell[0] + 8, center_cell[1] + 8, center_cell[2] - 8, center_cell[3] - 8), fill=(160, 168, 176, 190), width=2)
        draw.line((center_cell[0] + 8, center_cell[3] - 8, center_cell[2] - 8, center_cell[1] + 8), fill=(160, 168, 176, 190), width=2)
        draw.ellipse((gx0 + 2.5 * cell - 9, gy0 + 2.5 * cell - 9, gx0 + 2.5 * cell + 9, gy0 + 2.5 * cell + 9), fill=(0, 0, 0, 255))
    for label, grid in (("photon-like", left), ("broad / multi-tower", right)):
        gx0, _, gx1, _ = grid
        label_font = font(TIMES_ITALIC, 28)
        tw, _ = text_box(draw, label, label_font)
        draw.text((gx0 + (gx1 - gx0 - tw) / 2, y0 + 440), label, font=label_font, fill=MUTED)
    callout_label(draw, (x0 + 286, y0 + 374), "seed removed", LIGHT_MUTED, size=20)
    draw_formula_run(
        draw,
        (x0 + 54, y0 + 548),
        [("w", 31, 0, TIMES_BOLD), ("η", 22, 13, TIMES_BOLD), ("cogX", 18, -12, TIMES_BOLD), (" and w", 31, 0, TIMES_BOLD), ("φ", 22, 13, TIMES_BOLD), ("cogX", 18, -12, TIMES_BOLD), (" = seed-excluded spread", 31, 0, TIMES)],
        fill=INK,
    )
    draw_wrapped(draw, "Tests the surrounding shower after the hottest tower is removed.", (x0 + 54, y0 + 620), x1 - x0 - 100, font(TIMES, 28), fill=MUTED, line_gap=8)
    draw_wrapped(draw, "Photon-like: narrow shoulders around the core.", (x0 + 54, y0 + 718), x1 - x0 - 100, font(TIMES_ITALIC, 31), fill=BLUE, line_gap=8)


def draw_elongation_panel(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    title = "Is the shower stretched or split?"
    tw, _ = text_box(draw, title, font(TIMES_BOLD, 39))
    draw.text((x0 + (x1 - x0 - tw) / 2, y0 + 28), title, font=font(TIMES_BOLD, 39), fill=INK)

    compact = [
        [0.00, 0.02, 0.03, 0.02, 0.00],
        [0.01, 0.07, 0.16, 0.08, 0.01],
        [0.02, 0.34, 0.98, 0.38, 0.02],
        [0.01, 0.22, 0.58, 0.24, 0.01],
        [0.00, 0.03, 0.06, 0.03, 0.00],
    ]
    split = [
        [0.03, 0.16, 0.50, 0.34, 0.05],
        [0.06, 0.34, 0.90, 0.58, 0.10],
        [0.08, 0.20, 0.30, 0.24, 0.14],
        [0.10, 0.50, 0.36, 0.22, 0.08],
        [0.06, 0.82, 0.52, 0.18, 0.04],
    ]
    cell = 58
    left = draw_small_tower_grid(draw, (x0 + 82, y0 + 120), cell, compact)
    right = draw_small_tower_grid(draw, (x0 + 398, y0 + 120), cell, split)
    for idx, grid in enumerate((left, right)):
        gx0, gy0, _, _ = grid
        broad_region = (gx0 + cell, gy0, gx0 + 4 * cell, gy0 + 5 * cell)
        narrow_strip = (gx0 + cell, gy0 + 2 * cell, gx0 + 4 * cell, gy0 + 4 * cell)
        overlay = Image.new("RGBA", base.size, (0, 0, 0, 0))
        overlay_draw = ImageDraw.Draw(overlay, "RGBA")
        overlay_draw.rectangle(broad_region, fill=(*SPHENIX_BLUE, 10), outline=(*SPHENIX_BLUE, 140), width=3)
        overlay_draw.rectangle(narrow_strip, fill=(*PHOTON_DARK, 18), outline=(*PHOTON_DARK, 230), width=4)
        base.alpha_composite(overlay)
        draw.rectangle(broad_region, outline=(*SPHENIX_BLUE, 150), width=3)
        draw.rectangle(narrow_strip, outline=(*PHOTON_DARK, 238), width=4)
        for bx in (gx0 + cell, gx0 + 2 * cell, gx0 + 3 * cell, gx0 + 4 * cell):
            draw.line((bx, gy0, bx, gy0 + 5 * cell), fill=(164, 172, 181, 210), width=2)
        for by in (gy0, gy0 + cell, gy0 + 2 * cell, gy0 + 3 * cell, gy0 + 4 * cell, gy0 + 5 * cell):
            draw.line((gx0 + cell, by, gx0 + 4 * cell, by), fill=(164, 172, 181, 210), width=2)
        if idx == 0:
            callout_label(draw, (gx0 - 58, gy0 + 132), "3x2", PHOTON_DARK, size=20)
            draw.line((gx0 - 2, gy0 + 150, gx0 + cell, gy0 + 3 * cell), fill=(*PHOTON_DARK, 220), width=3)
            callout_label(draw, (gx0 - 58, gy0 + 252), "3x5", SPHENIX_BLUE, size=20)
            draw.line((gx0 - 2, gy0 + 270, gx0 + cell, gy0 + 5 * cell - 4), fill=(*SPHENIX_BLUE, 185), width=3)
    for label, grid in (("compact local", left), ("stretched / split", right)):
        gx0, _, gx1, _ = grid
        label_font = font(TIMES_ITALIC, 29)
        tw, _ = text_box(draw, label, label_font)
        draw.text((gx0 + (gx1 - gx0 - tw) / 2, y0 + 432), label, font=label_font, fill=MUTED)
    draw_formula_run(
        draw,
        (x0 + 54, y0 + 566),
        [("e32/e35", 31, 0, TIMES_BOLD), (" = energy in 3x2 strip over local 3x5 region", 31, 0, TIMES)],
        fill=INK,
    )
    draw_wrapped(draw, "Compares energy confined to the narrow strip with energy spread across the broader local region.", (x0 + 54, y0 + 640), x1 - x0 - 100, font(TIMES, 28), fill=MUTED, line_gap=8)
    draw_wrapped(draw, "Photon-like: compact energy stays confined to the narrow strip.", (x0 + 54, y0 + 736), x1 - x0 - 100, font(TIMES_ITALIC, 31), fill=BLUE, line_gap=8)


def draw_bdt_flow(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    rounded_panel(draw, box, fill=PANEL)
    x0, y0, x1, y1 = box
    draw.text((x0 + 34, y0 + 24), "Shower-shape information", font=font(TIMES_BOLD, 30), fill=BLUE)
    steps = [("EMCal widths", SPHENIX_BLUE), ("energy ratios", TEAL), ("cluster kinematics", PHOTON_DARK), ("BDT score", PHOTON)]
    x = x0 + 46
    y = y0 + 96
    for idx, (label, color) in enumerate(steps):
        w = 164 if idx < 3 else 142
        draw.rounded_rectangle((x, y, x + w, y + 78), radius=8, fill=(255, 255, 255, 255), outline=(*CARD_EDGE, 255), width=2)
        draw.text((x + 18, y + 24), label, font=font(TIMES_BOLD if idx == 3 else TIMES, 22), fill=INK)
        draw.rectangle((x, y, x + 8, y + 78), fill=(*color, 255))
        if idx < len(steps) - 1:
            draw.line((x + w + 12, y + 39, x + w + 48, y + 39), fill=(150, 170, 188), width=3)
            draw.polygon([(x + w + 48, y + 39), (x + w + 36, y + 31), (x + w + 36, y + 47)], fill=(150, 170, 188))
        x += w + 54
    draw.text((x0 + 48, y1 - 54), "High score = photon-like cluster", font=font(TIMES_ITALIC, 27), fill=MUTED)


def draw_arrow(
    draw: ImageDraw.ImageDraw,
    start: tuple[int, int],
    end: tuple[int, int],
    fill: tuple[int, int, int] = (151, 169, 188),
    width: int = 5,
) -> None:
    sx, sy = start
    ex, ey = end
    draw.line((sx, sy, ex, ey), fill=(*fill, 255), width=width)
    angle = math.atan2(ey - sy, ex - sx)
    head = 18
    left = (ex - head * math.cos(angle - math.pi / 7), ey - head * math.sin(angle - math.pi / 7))
    right = (ex - head * math.cos(angle + math.pi / 7), ey - head * math.sin(angle + math.pi / 7))
    draw.polygon([(ex, ey), left, right], fill=(*fill, 255))


def draw_manual_cuts_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    draw.text((x0 + 34, y0 + 26), "Manual rectangular cuts", font=font(TIMES_BOLD, 34), fill=INK)
    draw_wrapped(draw, "Understandable, but each gate is mostly tuned one variable at a time.", (x0 + 34, y0 + 78), x1 - x0 - 68, font(TIMES_ITALIC, 25), fill=MUTED, line_gap=5)
    gates = [
        [("w", 28, 0, TIMES_BOLD), ("η", 19, 11, TIMES_BOLD), ("cogX", 15, -10, TIMES_BOLD), (" < a", 28, 0, TIMES_BOLD)],
        [("E", 28, 0, TIMES_BOLD), ("11", 18, 11, TIMES_BOLD), ("/E", 28, 0, TIMES_BOLD), ("33", 18, 11, TIMES_BOLD), (" > b", 28, 0, TIMES_BOLD)],
        [("E", 28, 0, TIMES_BOLD), ("3x2", 18, 11, TIMES_BOLD), ("/E", 28, 0, TIMES_BOLD), ("3x5", 18, 11, TIMES_BOLD), (" > c", 28, 0, TIMES_BOLD)],
        [("E", 28, 0, TIMES_BOLD), ("T", 18, 11, TIMES_BOLD), ("(1)", 18, -6, TIMES_BOLD), (" in range", 28, 0, TIMES_BOLD)],
    ]
    gy = y0 + 166
    for idx, parts in enumerate(gates):
        y = gy + idx * 92
        draw.rounded_rectangle((x0 + 42, y, x1 - 42, y + 66), radius=8, fill=(247, 249, 251, 255), outline=(218, 226, 235, 255), width=2)
        gate_color = PHOTON if idx % 2 == 0 else SPHENIX_BLUE
        draw.rounded_rectangle((x0 + 42, y, x0 + 54, y + 66), radius=5, fill=(*gate_color, 255))
        draw.ellipse((x0 + 75, y + 18, x0 + 106, y + 49), outline=(*MUTED, 255), width=3)
        draw.line((x0 + 82, y + 34, x0 + 91, y + 44), fill=(*TEAL, 255), width=4)
        draw.line((x0 + 91, y + 44, x0 + 101, y + 23), fill=(*TEAL, 255), width=4)
        draw_formula_run(draw, (x0 + 126, y + 18), parts, fill=INK)
    draw.rounded_rectangle((x0 + 54, y1 - 116, x1 - 54, y1 - 46), radius=8, fill=(255, 248, 229, 255), outline=(238, 220, 172, 255), width=2)
    draw_wrapped(draw, "Fixed gates: clear, but rigid.", (x0 + 78, y1 - 98), x1 - x0 - 156, font(TIMES_BOLD, 25), fill=BLUE, line_gap=2)


def tree_node(draw: ImageDraw.ImageDraw, center: tuple[int, int], text_parts: list[tuple[str, int, int, Path]], width: int = 190) -> tuple[int, int, int, int]:
    cx, cy = center
    box = (cx - width // 2, cy - 34, cx + width // 2, cy + 34)
    draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 255), outline=(*SPHENIX_BLUE, 210), width=3)
    draw_formula_run(draw, (box[0] + 18, box[1] + 18), text_parts, fill=INK)
    return box


def leaf_box(draw: ImageDraw.ImageDraw, center: tuple[int, int], label: str, accent: tuple[int, int, int], width: int = 154) -> None:
    cx, cy = center
    box = (cx - width // 2, cy - 34, cx + width // 2, cy + 42)
    draw.rounded_rectangle(box, radius=8, fill=(255, 255, 255, 255), outline=(*accent, 230), width=3)
    draw_wrapped(draw, label, (box[0] + 14, box[1] + 13), width - 28, font(TIMES_BOLD, 20), fill=INK, line_gap=1)


def draw_single_tree_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    draw.text((x0 + 34, y0 + 26), "One learned decision tree", font=font(TIMES_BOLD, 34), fill=INK)
    draw_wrapped(draw, "A tree learns which shower-shape question to ask next.", (x0 + 34, y0 + 78), x1 - x0 - 68, font(TIMES_ITALIC, 25), fill=MUTED, line_gap=5)
    cx = (x0 + x1) // 2
    root = (cx, y0 + 172)
    left = (cx - 170, y0 + 316)
    right = (cx + 170, y0 + 316)
    ll = (cx - 260, y0 + 470)
    lr = (cx - 80, y0 + 470)
    rl = (cx + 80, y0 + 470)
    rr = (cx + 260, y0 + 470)
    for start, end, label in [(root, left, "yes"), (root, right, "no"), (left, ll, "yes"), (left, lr, "no"), (right, rl, "yes"), (right, rr, "no")]:
        draw.line((start[0], start[1] + 36, end[0], end[1] - 36), fill=(168, 184, 199, 255), width=4)
        mx = (start[0] + end[0]) // 2
        my = (start[1] + end[1]) // 2
        callout_label(draw, (mx - 15, my - 17), label, TEAL if label == "yes" else LIGHT_MUTED, bg=(255, 255, 255), size=18)
    tree_node(draw, root, [("w", 27, 0, TIMES_BOLD), ("η", 18, 10, TIMES_BOLD), ("cogX", 14, -9, TIMES_BOLD), (" < a?", 27, 0, TIMES_BOLD)], width=210)
    tree_node(draw, left, [("E", 26, 0, TIMES_BOLD), ("11", 17, 10, TIMES_BOLD), ("/E", 26, 0, TIMES_BOLD), ("33", 17, 10, TIMES_BOLD), (" > b?", 26, 0, TIMES_BOLD)], width=214)
    tree_node(draw, right, [("E", 26, 0, TIMES_BOLD), ("3x2", 17, 10, TIMES_BOLD), ("/E", 26, 0, TIMES_BOLD), ("3x5", 17, 10, TIMES_BOLD), (" > c?", 26, 0, TIMES_BOLD)], width=218)
    leaf_box(draw, ll, "photon-like leaf", SPHENIX_BLUE)
    leaf_box(draw, lr, "mixed leaf", PHOTON_DARK)
    leaf_box(draw, rl, "mixed leaf", PHOTON_DARK)
    leaf_box(draw, rr, "background-like leaf", (205, 75, 55), width=176)
    draw.rounded_rectangle((x0 + 58, y1 - 116, x1 - 58, y1 - 46), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw_wrapped(draw, "Each final leaf assigns a small signal-like score.", (x0 + 82, y1 - 98), x1 - x0 - 164, font(TIMES_BOLD, 25), fill=BLUE, line_gap=2)


def draw_mini_tree(draw: ImageDraw.ImageDraw, origin: tuple[int, int], scale: float = 1.0, alpha: int = 170) -> None:
    x, y = origin
    pts = [
        (x, y),
        (x - 42 * scale, y + 54 * scale),
        (x + 42 * scale, y + 54 * scale),
        (x - 66 * scale, y + 106 * scale),
        (x - 18 * scale, y + 106 * scale),
        (x + 18 * scale, y + 106 * scale),
        (x + 66 * scale, y + 106 * scale),
    ]
    for a, b in [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)]:
        draw.line((pts[a][0], pts[a][1], pts[b][0], pts[b][1]), fill=(119, 139, 157, alpha), width=max(2, round(3 * scale)))
    for idx, p in enumerate(pts[:3]):
        r = 12 * scale
        draw.ellipse((p[0] - r, p[1] - r, p[0] + r, p[1] + r), fill=(*SPHENIX_BLUE, alpha), outline=(255, 255, 255, alpha), width=2)
    for p in pts[3:]:
        r = 10 * scale
        draw.rounded_rectangle((p[0] - r, p[1] - r, p[0] + r, p[1] + r), radius=3, fill=(*PHOTON, alpha), outline=(255, 255, 255, alpha), width=1)


def draw_boosted_score_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    draw.text((x0 + 34, y0 + 26), "Boosted decision tree", font=font(TIMES_BOLD, 34), fill=INK)
    draw_wrapped(draw, "Many shallow trees vote together; their outputs are combined into one ranking score.", (x0 + 34, y0 + 78), x1 - x0 - 68, font(TIMES_ITALIC, 25), fill=MUTED, line_gap=5)

    flow_y = y0 + 246
    origins = [(x0 + 100, y0 + 214), (x0 + 200, y0 + 196), (x0 + 300, y0 + 214), (x0 + 400, y0 + 198)]
    alphas = [220, 200, 180, 160]
    for idx, (origin, alpha) in enumerate(zip(origins, alphas)):
        draw_mini_tree(draw, origin, scale=0.76 - idx * 0.015, alpha=alpha)
        draw.text((origin[0] - 24, y0 + 318), f"T{idx + 1}", font=font(TIMES_BOLD, 19), fill=LIGHT_MUTED)
        if idx < len(origins) - 1:
            draw.text((origin[0] + 46, y0 + 314), "+", font=font(TIMES_BOLD, 29), fill=LIGHT_MUTED)

    combiner = (x0 + 508, flow_y)
    score_box = (x1 - 174, flow_y - 40, x1 - 48, flow_y + 40)
    draw_arrow(draw, (x0 + 468, flow_y), (combiner[0] - 36, flow_y), fill=(151, 169, 188), width=4)
    draw.ellipse((combiner[0] - 34, combiner[1] - 34, combiner[0] + 34, combiner[1] + 34), fill=(239, 246, 250, 255), outline=(*SPHENIX_BLUE, 230), width=3)
    draw.text((combiner[0] - 18, combiner[1] - 21), "Σ", font=font(TIMES_BOLD, 38), fill=BLUE)
    draw.text((combiner[0] - 58, combiner[1] + 48), "weighted sum", font=font(TIMES_ITALIC, 20), fill=MUTED)
    draw_arrow(draw, (combiner[0] + 42, flow_y), (score_box[0] - 18, flow_y), fill=(151, 169, 188), width=4)
    draw.rounded_rectangle(score_box, radius=13, fill=(*PHOTON, 230), outline=(*PHOTON_DARK, 230), width=3)
    draw.text((score_box[0] + 26, score_box[1] + 14), "BDT", font=font(TIMES_BOLD, 28), fill=INK)
    draw.text((score_box[0] + 22, score_box[1] + 45), "score", font=font(TIMES_BOLD, 23), fill=INK)

    draw.text((x0 + 72, y0 + 358), "many weak learners", font=font(TIMES_ITALIC, 22), fill=MUTED)
    score_caption_font = font(TIMES_ITALIC, 20)
    score_caption = "one ranking score"
    score_caption_w, _ = text_box(draw, score_caption, score_caption_font)
    draw.text((x1 - 72 - score_caption_w, y0 + 358), score_caption, font=score_caption_font, fill=MUTED)

    bar_x0 = x0 + 66
    bar_x1 = x1 - 66
    bar_y = y0 + 440
    draw.rounded_rectangle((bar_x0, bar_y, bar_x1, bar_y + 34), radius=17, fill=(232, 236, 240, 255))
    for i in range(bar_x1 - bar_x0):
        t = i / max(1, bar_x1 - bar_x0 - 1)
        r = int(213 * (1 - t) + SPHENIX_BLUE[0] * t)
        g = int(88 * (1 - t) + SPHENIX_BLUE[1] * t)
        b = int(67 * (1 - t) + SPHENIX_BLUE[2] * t)
        draw.line((bar_x0 + i, bar_y, bar_x0 + i, bar_y + 34), fill=(r, g, b, 210), width=1)
    draw.text((bar_x0, bar_y + 52), "background-like", font=font(TIMES_BOLD, 23), fill=(180, 70, 56))
    draw.text((bar_x1 - 124, bar_y + 52), "photon-like", font=font(TIMES_BOLD, 23), fill=SPHENIX_BLUE)
    draw.text((bar_x0 - 6, bar_y - 36), "0", font=font(TIMES_BOLD, 24), fill=MUTED)
    draw.text((bar_x1 - 10, bar_y - 36), "1", font=font(TIMES_BOLD, 24), fill=MUTED)
    draw.rounded_rectangle((x0 + 58, y1 - 116, x1 - 58, y1 - 46), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw_wrapped(draw, "Score ranks clusters; it is not an absolute photon probability.", (x0 + 82, y1 - 98), x1 - x0 - 164, font(TIMES_BOLD, 21), fill=BLUE, line_gap=1)


def draw_bdt_input_ribbon(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(247, 250, 252, 255), outline=(220, 228, 236, 255), width=2)
    draw.text((x0 + 30, y0 + 27), "Inputs to the learned classifier", font=font(TIMES_BOLD, 32), fill=BLUE)
    draw_wrapped(
        draw,
        "shower-shape cues and event kinematics become one ranking score",
        (x0 + 30, y0 + 70),
        482,
        font(TIMES_ITALIC, 24),
        fill=MUTED,
        line_gap=3,
    )

    def gear(center: tuple[int, int], r: int, accent: tuple[int, int, int]) -> None:
        cx, cy = center
        for idx in range(12):
            angle = idx * math.tau / 12
            tx = cx + math.cos(angle) * (r + 6)
            ty = cy + math.sin(angle) * (r + 6)
            draw.line((cx + math.cos(angle) * (r - 2), cy + math.sin(angle) * (r - 2), tx, ty), fill=(*accent, 185), width=3)
        draw.ellipse((cx - r, cy - r, cx + r, cy + r), fill=(247, 250, 252, 255), outline=(*accent, 230), width=3)
        draw.ellipse((cx - r // 3, cy - r // 3, cx + r // 3, cy + r // 3), fill=(255, 255, 255, 255), outline=(*accent, 190), width=2)

    def paste_tilted_fuel_can(dest: Image.Image, origin: tuple[int, int]) -> tuple[int, int, int, int]:
        scale = 4
        iw, ih = 360, 178
        icon = Image.new("RGBA", (iw * scale, ih * scale), (0, 0, 0, 0))
        idraw = ImageDraw.Draw(icon, "RGBA")

        def pt(p: tuple[int, int]) -> tuple[int, int]:
            return (p[0] * scale, p[1] * scale)

        def pts(points: list[tuple[int, int]]) -> list[tuple[int, int]]:
            return [pt(p) for p in points]

        def box(b: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
            return tuple(v * scale for v in b)

        panel_fill = (247, 250, 252, 255)
        outline = (28, 34, 40, 245)
        can_orange = (248, 102, 22, 255)
        can_light = (255, 139, 43, 255)
        can_dark = (206, 78, 18, 255)
        spout_yellow = (255, 211, 42, 255)
        cap_dark = (55, 65, 62, 255)

        def bezier_points(
            p0: tuple[int, int],
            p1: tuple[int, int],
            p2: tuple[int, int],
            p3: tuple[int, int],
            steps: int = 28,
        ) -> list[tuple[int, int]]:
            curve = []
            for idx in range(steps + 1):
                t = idx / steps
                u = 1.0 - t
                x = u**3 * p0[0] + 3 * u * u * t * p1[0] + 3 * u * t * t * p2[0] + t**3 * p3[0]
                y = u**3 * p0[1] + 3 * u * u * t * p1[1] + 3 * u * t * t * p2[1] + t**3 * p3[1]
                curve.append((round(x), round(y)))
            return curve

        def teardrop(cx: int, top_y: int, w: int, h: int) -> list[tuple[int, int]]:
            left = bezier_points(
                (cx, top_y),
                (round(cx - 0.55 * w), round(top_y + 0.32 * h)),
                (round(cx - 0.56 * w), round(top_y + 0.78 * h)),
                (cx, top_y + h),
                steps=18,
            )
            right = bezier_points(
                (cx, top_y + h),
                (round(cx + 0.56 * w), round(top_y + 0.78 * h)),
                (round(cx + 0.55 * w), round(top_y + 0.32 * h)),
                (cx, top_y),
                steps=18,
            )
            return left + right

        def paste_reference_oil_drop(cx: int, top_y: int, w: int, h: int) -> None:
            pad = 10
            local_w = w + 2 * pad
            local_h = h + 2 * pad
            drop_img = Image.new("RGBA", (local_w * scale, local_h * scale), (0, 0, 0, 0))
            mask = Image.new("L", drop_img.size, 0)
            mask_draw = ImageDraw.Draw(mask)

            def lpt(p: tuple[int, int]) -> tuple[int, int]:
                return (p[0] * scale, p[1] * scale)

            def lpts(points: list[tuple[int, int]]) -> list[tuple[int, int]]:
                return [lpt(p) for p in points]

            dcx = pad + w // 2
            dtop = pad
            dbottom = pad + h
            left_top = bezier_points(
                (dcx, dtop),
                (round(dcx - 0.18 * w), round(dtop + 0.14 * h)),
                (round(dcx - 0.53 * w), round(dtop + 0.34 * h)),
                (round(dcx - 0.52 * w), round(dtop + 0.58 * h)),
                steps=18,
            )
            left_bottom = bezier_points(
                left_top[-1],
                (round(dcx - 0.52 * w), round(dtop + 0.86 * h)),
                (round(dcx - 0.31 * w), round(dbottom + 0.02 * h)),
                (dcx, dbottom),
                steps=18,
            )
            right_bottom = bezier_points(
                (dcx, dbottom),
                (round(dcx + 0.31 * w), round(dbottom + 0.02 * h)),
                (round(dcx + 0.52 * w), round(dtop + 0.86 * h)),
                (round(dcx + 0.52 * w), round(dtop + 0.58 * h)),
                steps=18,
            )
            right_top = bezier_points(
                right_bottom[-1],
                (round(dcx + 0.53 * w), round(dtop + 0.34 * h)),
                (round(dcx + 0.18 * w), round(dtop + 0.14 * h)),
                (dcx, dtop),
                steps=18,
            )
            shape = left_top + left_bottom + right_bottom + right_top
            mask_draw.polygon(lpts(shape), fill=255)
            mask = mask.filter(ImageFilter.GaussianBlur(0.38 * scale))

            gradient = Image.new("RGBA", drop_img.size, (0, 0, 0, 0))
            grad_draw = ImageDraw.Draw(gradient, "RGBA")
            for yy in range(gradient.height):
                t = yy / max(1, gradient.height - 1)
                r = round(255 * (1 - t) + 245 * t)
                g = round(205 * (1 - t) + 166 * t)
                b = round(54 * (1 - t) + 31 * t)
                grad_draw.line((0, yy, gradient.width, yy), fill=(r, g, b, 255))
            gradient.putalpha(mask)

            edge = Image.new("RGBA", drop_img.size, (0, 0, 0, 0))
            edge_draw = ImageDraw.Draw(edge, "RGBA")
            edge_draw.line(lpts(shape + [shape[0]]), fill=(232, 157, 28, 105), width=2 * scale, joint="curve")
            edge = edge.filter(ImageFilter.GaussianBlur(0.18 * scale))
            gradient.alpha_composite(edge)

            highlight_draw = ImageDraw.Draw(gradient, "RGBA")
            highlight = bezier_points(
                (pad + round(0.34 * w), pad + round(0.50 * h)),
                (pad + round(0.17 * w), pad + round(0.62 * h)),
                (pad + round(0.20 * w), pad + round(0.78 * h)),
                (pad + round(0.38 * w), pad + round(0.76 * h)),
                steps=16,
            )
            highlight_draw.line(lpts(highlight), fill=(255, 255, 255, 190), width=4 * scale, joint="curve")
            highlight_draw.line(lpts(highlight), fill=(255, 255, 255, 235), width=2 * scale, joint="curve")

            paste_x = round((cx - w / 2 - pad) * scale)
            paste_y = (top_y - pad) * scale
            icon.alpha_composite(gradient, (paste_x, paste_y))

        # Curved right-facing pouring spout, drawn first so the collar can lock it to the can.
        spout_path = pts(bezier_points((235, 55), (258, 21), (306, 16), (343, 42)))
        idraw.line(spout_path, fill=outline, width=18 * scale, joint="curve")
        idraw.line(spout_path, fill=spout_yellow, width=12 * scale, joint="curve")
        tip_angle = math.atan2(42 - 16, 343 - 306)
        tip_tangent = (math.cos(tip_angle), math.sin(tip_angle))
        tip_normal = (-tip_tangent[1], tip_tangent[0])
        tip_center = (341.5, 42.0)

        def tip_point(normal_scale: float, tangent_scale: float) -> tuple[int, int]:
            return (
                round(tip_center[0] + tip_normal[0] * normal_scale + tip_tangent[0] * tangent_scale),
                round(tip_center[1] + tip_normal[1] * normal_scale + tip_tangent[1] * tangent_scale),
            )

        def tip_oval(normal_radius: float, tangent_radius: float, tangent_offset: float = 0.0, n: int = 40) -> list[tuple[int, int]]:
            oval = []
            for idx in range(n):
                theta = math.tau * idx / n
                normal_scale = math.cos(theta) * normal_radius
                tangent_scale = math.sin(theta) * tangent_radius + tangent_offset
                oval.append(tip_point(normal_scale, tangent_scale))
            return oval

        cap_outer = tip_oval(8.5, 4.0, tangent_offset=0.6)
        cap_inner = tip_oval(5.2, 2.1, tangent_offset=0.9)
        idraw.polygon(pts(cap_outer), fill=outline)
        idraw.polygon(pts(cap_inner), fill=spout_yellow)
        idraw.line(pts(cap_outer + [cap_outer[0]]), fill=outline, width=2 * scale, joint="curve")
        idraw.line(pts(cap_inner + [cap_inner[0]]), fill=(255, 241, 92, 255), width=1 * scale, joint="curve")
        paste_reference_oil_drop(338, 62, 27, 41)

        body_poly = [(56, 56), (205, 54), (232, 70), (235, 138), (82, 156), (58, 144)]
        idraw.polygon(pts(body_poly), fill=can_orange)
        idraw.line(pts(body_poly + [body_poly[0]]), fill=outline, width=4 * scale, joint="curve")
        idraw.polygon(pts([(56, 56), (82, 56), (84, 156), (58, 144)]), fill=can_light)
        idraw.line(pts([(82, 61), (82, 146)]), fill=can_dark, width=3 * scale)

        handle_ridge = [(106, 58), (130, 36), (218, 39), (229, 58), (211, 68), (122, 66)]
        handle_cutout = [(145, 46), (199, 47), (192, 59), (153, 58)]
        idraw.polygon(pts(handle_ridge), fill=can_orange)
        idraw.line(pts(handle_ridge + [handle_ridge[0]]), fill=outline, width=3 * scale, joint="curve")
        idraw.polygon(pts(handle_cutout), fill=panel_fill)
        idraw.line(pts(handle_cutout + [handle_cutout[0]]), fill=outline, width=2 * scale, joint="curve")

        collar_axis = [(216, 60), (243, 76)]
        idraw.line(pts(collar_axis), fill=outline, width=24 * scale)
        idraw.line(pts(collar_axis), fill=cap_dark, width=18 * scale)
        for offset in (-4, 5):
            idraw.line(
                pts([(219 + offset, 57 + offset // 3), (238 + offset, 69 + offset // 3)]),
                fill=(119, 130, 126, 210),
                width=2 * scale,
            )
        resample = Image.Resampling.LANCZOS if hasattr(Image, "Resampling") else Image.LANCZOS
        icon = icon.resize((iw, ih), resample)
        icon = icon.rotate(5, expand=True, resample=Image.Resampling.BICUBIC if hasattr(Image, "Resampling") else Image.BICUBIC)
        icon = icon.resize((int(icon.width * 0.80), int(icon.height * 0.80)), resample)
        final_draw = ImageDraw.Draw(icon, "RGBA")
        label_fill = (255, 246, 219, 255)
        label_stroke = (100, 30, 7, 245)
        label_lines = [
            ("feature", font(TIMES_BOLD, 24), 71),
            ("inputs", font(TIMES_BOLD, 24), 98),
        ]
        label_center_x = 128
        for label, label_font, label_y in label_lines:
            tw, _ = text_box(final_draw, label, label_font)
            final_draw.text(
                (label_center_x - tw / 2, label_y),
                label,
                font=label_font,
                fill=label_fill,
                stroke_width=1,
                stroke_fill=label_stroke,
            )
        dest.paste(icon, origin, icon)
        right = origin[0] + icon.width
        return (right - 44, origin[1] + 58, right - 4, origin[1] + 88)

    nozzle = paste_tilted_fuel_can(base, (x0 + 520, y0 - 4))

    engine = (x1 - 420, y0 + 25, x1 - 112, y0 + 132)

    stream_x0 = nozzle[2] + 2
    stream_x1 = engine[0] - 8
    stream_y0 = y0 + 51
    stream_y1 = y0 + 109
    draw.rounded_rectangle((stream_x0, stream_y0, stream_x1, stream_y1), radius=29, fill=(*PHOTON, 45), outline=(*PHOTON_DARK, 140), width=2)
    for yy in (y0 + 66, y0 + 82, y0 + 98):
        draw.line((stream_x0 + 14, yy, stream_x1 - 34, yy - 6), fill=(*PHOTON_DARK, 86), width=2)
    draw.polygon([(stream_x1, y0 + 80), (stream_x1 - 22, y0 + 66), (stream_x1 - 22, y0 + 94)], fill=(*PHOTON_DARK, 180))

    token_specs = [
        (
            [("et1", 23, 0, TIMES_BOLD)],
            PHOTON_DARK,
            76,
            0.13,
        ),
        (
            [("E", 22, 0, TIMES_BOLD), ("11", 14, 9, TIMES_BOLD), ("/E", 22, 0, TIMES_BOLD), ("33", 14, 9, TIMES_BOLD)],
            PHOTON_DARK,
            128,
            0.31,
        ),
        (
            [("w", 23, 0, TIMES_BOLD), ("η|φ", 14, 9, TIMES_BOLD), ("cogx", 12, -8, TIMES_BOLD)],
            SPHENIX_BLUE,
            134,
            0.49,
        ),
        (
            [("E", 22, 0, TIMES_BOLD), ("3x2", 14, 9, TIMES_BOLD), ("/E", 22, 0, TIMES_BOLD), ("3x5", 14, 9, TIMES_BOLD)],
            TEAL,
            144,
            0.67,
        ),
        (
            [("η", 25, 0, TIMES_BOLD)],
            LIGHT_MUTED,
            64,
            0.82,
        ),
        (
            [("z", 23, 0, TIMES_BOLD), ("vtx", 14, 9, TIMES_BOLD)],
            LIGHT_MUTED,
            84,
            0.94,
        ),
    ]
    token_y_lanes = [y0 + 55, y0 + 68, y0 + 52, y0 + 68, y0 + 54, y0 + 67]
    available_w = stream_x1 - stream_x0 - 104
    token_x0 = stream_x0 + 52
    token_step = available_w / max(1, len(token_specs) - 1)
    tokens = []
    for idx, ((parts, color, bw, phase), by) in enumerate(zip(token_specs, token_y_lanes)):
        bx = round(token_x0 + idx * token_step - bw / 2 + (phase - 0.5) * 14)
        bx = max(stream_x0 + 20, min(bx, stream_x1 - bw - 26))
        tokens.append((bx, by, bw, parts, color))
    for bx, by, bw, parts, color in tokens:
        draw.rounded_rectangle((bx, by, bx + bw, by + 42), radius=21, fill=(255, 255, 255, 245), outline=(*color, 230), width=3)
        formula_w = sum(text_box(draw, text, font(font_path, size))[0] for text, size, _, font_path in parts)
        draw_formula_run(draw, (round(bx + (bw - formula_w) / 2), by + 9), parts, fill=INK)

    shadow(base, engine, radius=10)
    draw.rounded_rectangle(engine, radius=14, fill=(255, 255, 255, 255), outline=(*BLUE, 230), width=3)
    draw.rectangle((engine[0], engine[1] + 8, engine[0] + 14, engine[3] - 8), fill=(*BLUE, 235))
    draw.text((engine[0] + 32, engine[1] + 18), "BDT engine", font=font(TIMES_BOLD, 27), fill=INK)
    draw.text((engine[0] + 32, engine[1] + 58), "learned classifier", font=font(TIMES_ITALIC, 19), fill=MUTED)
    gear((engine[0] + 210, y0 + 78), 25, SPHENIX_BLUE)
    gear((engine[0] + 254, y0 + 78), 20, PHOTON_DARK)
    draw.rounded_rectangle((engine[0] + 188, y0 + 104, engine[0] + 276, y0 + 114), radius=5, fill=(*PHOTON, 120))


def draw_pipeline_stage(
    base: Image.Image,
    box: tuple[int, int, int, int],
    label: str,
    title: str,
    subtitle: str,
    accent: tuple[int, int, int],
) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    shadow(base, box)
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    draw.rounded_rectangle((x0, y0, x0 + 14, y1), radius=7, fill=(*accent, 255))
    draw.text((x0 + 34, y0 + 24), label, font=font(TIMES_BOLD, 20), fill=accent)
    draw.text((x0 + 34, y0 + 54), title, font=font(TIMES_BOLD, 29), fill=INK)
    draw_wrapped(draw, subtitle, (x0 + 34, y0 + 94), x1 - x0 - 68, font(TIMES_ITALIC, 21), fill=MUTED, line_gap=3)


def draw_emcal_cluster_stage(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    draw_pipeline_stage(
        base,
        box,
        "A",
        "What the BDT sees",
        "measured EMCal cluster features",
        SPHENIX_BLUE,
    )

    gx0, gy0 = x0 + 54, y0 + 162
    cell = 42
    values = [
        [0.10, 0.18, 0.26, 0.16, 0.08],
        [0.16, 0.36, 0.62, 0.32, 0.14],
        [0.22, 0.58, 1.00, 0.50, 0.18],
        [0.10, 0.32, 0.48, 0.27, 0.10],
        [0.04, 0.10, 0.16, 0.09, 0.04],
    ]
    for row in range(5):
        for col in range(5):
            x = gx0 + col * cell
            y = gy0 + row * cell
            t = values[row][col]
            r = int(240 * (1 - t) + PHOTON[0] * t)
            g = int(247 * (1 - t) + PHOTON[1] * t)
            b = int(251 * (1 - t) + PHOTON[2] * t)
            draw.rounded_rectangle((x, y, x + cell - 4, y + cell - 4), radius=5, fill=(r, g, b, 255), outline=(201, 214, 226, 255), width=2)
    draw.rectangle((gx0 + cell, gy0 + cell, gx0 + 4 * cell - 4, gy0 + 4 * cell - 4), outline=(*SPHENIX_BLUE, 200), width=3)
    draw.rectangle((gx0 + 2 * cell, gy0 + 2 * cell, gx0 + 3 * cell - 4, gy0 + 3 * cell - 4), outline=(*PHOTON_DARK, 240), width=4)
    callout_label(draw, (gx0 + 8, gy0 - 34), "EMCal tower cluster", BLUE, size=20)

    callouts = [
        ((gx0 + 3 * cell, gy0 + 2 * cell), "core compactness", PHOTON_DARK),
        ((gx0 + 4 * cell, gy0 + 1 * cell), "shoulders", SPHENIX_BLUE),
        ((gx0 + 4 * cell, gy0 + 4 * cell), "elongation / split", TEAL),
    ]
    for idx, (start, label, color) in enumerate(callouts):
        lx = x0 + 316
        ly = gy0 + 24 + idx * 74
        draw.line((start[0], start[1], lx - 14, ly + 18), fill=(*color, 180), width=3)
        draw.rounded_rectangle((lx, ly, x1 - 36, ly + 42), radius=8, fill=(247, 250, 252, 255), outline=(*color, 190), width=2)
        draw.text((lx + 16, ly + 11), label, font=font(TIMES_BOLD, 19), fill=INK)

    fv = (x0 + 34, y0 + 452, x1 - 34, y0 + 612)
    draw.rounded_rectangle(fv, radius=10, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw.text((fv[0] + 22, fv[1] + 18), "Feature vector", font=font(TIMES_BOLD, 25), fill=BLUE)
    draw_formula_run(
        draw,
        (fv[0] + 22, fv[1] + 62),
        [
            ("x = [ core, shoulders,", 24, 0, TIMES_BOLD),
        ],
        fill=INK,
    )
    draw_formula_run(
        draw,
        (fv[0] + 22, fv[1] + 100),
        [
            ("elongation, E", 24, 0, TIMES_BOLD),
            ("T", 15, 9, TIMES_BOLD),
            (", η, z", 24, 0, TIMES_BOLD),
            ("vtx", 15, 9, TIMES_BOLD),
            (" ]", 24, 0, TIMES_BOLD),
        ],
        fill=INK,
    )

    inset = (x0 + 34, y1 - 178, x1 - 34, y1 - 34)
    draw.rounded_rectangle(inset, radius=10, fill=(255, 248, 229, 255), outline=(238, 220, 172, 255), width=2)
    draw.text((inset[0] + 18, inset[1] + 16), "Manual cuts contrast", font=font(TIMES_BOLD, 22), fill=BLUE)
    draw_wrapped(draw, "fixed gates are clear, but rigid", (inset[0] + 18, inset[1] + 50), inset[2] - inset[0] - 36, font(TIMES_ITALIC, 20), fill=MUTED, line_gap=2)
    mini_y = inset[1] + 92
    for idx, label in enumerate(["width < a", "ratio > b", "ratio > c"]):
        x = inset[0] + 22 + idx * 142
        draw.rounded_rectangle((x, mini_y, x + 118, mini_y + 34), radius=7, fill=(255, 255, 255, 210), outline=(*PHOTON_DARK, 160), width=2)
        tw, _ = text_box(draw, label, font(TIMES_BOLD, 16))
        draw.text((x + (118 - tw) / 2, mini_y + 8), label, font=font(TIMES_BOLD, 16), fill=INK)


def draw_learned_tree_stage(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    draw_pipeline_stage(
        base,
        box,
        "B",
        "Learned shower-shape questions",
        "one tree follows conditional cuts",
        PHOTON_DARK,
    )
    cx = (x0 + x1) // 2
    root = (cx, y0 + 188)
    left = (cx - 130, y0 + 330)
    right = (cx + 130, y0 + 330)
    ll = (cx - 168, y0 + 500)
    lr = (cx, y0 + 500)
    rr = (cx + 168, y0 + 500)
    leaves = [(cx - 186, y0 + 654), (cx, y0 + 654), (cx + 186, y0 + 654)]

    for start, end, label in [
        (root, left, "yes"),
        (root, right, "no"),
        (left, ll, "yes"),
        (left, lr, "no"),
        (right, rr, "no"),
        (ll, leaves[0], "+"),
        (lr, leaves[1], "0"),
        (rr, leaves[2], "-"),
    ]:
        draw.line((start[0], start[1] + 36, end[0], end[1] - 36), fill=(168, 184, 199, 255), width=4)
        mx = (start[0] + end[0]) // 2
        my = (start[1] + end[1]) // 2
        if label in {"yes", "no"}:
            callout_label(draw, (mx - 16, my - 18), label, TEAL if label == "yes" else LIGHT_MUTED, bg=(255, 255, 255), size=17)

    tree_node(draw, root, [("w", 25, 0, TIMES_BOLD), ("η", 16, 10, TIMES_BOLD), ("cogX", 13, -8, TIMES_BOLD), (" < a?", 25, 0, TIMES_BOLD)], width=202)
    tree_node(draw, left, [("E", 24, 0, TIMES_BOLD), ("11", 15, 9, TIMES_BOLD), ("/E", 24, 0, TIMES_BOLD), ("33", 15, 9, TIMES_BOLD), (" > b?", 24, 0, TIMES_BOLD)], width=205)
    tree_node(draw, right, [("E", 24, 0, TIMES_BOLD), ("3x2", 15, 9, TIMES_BOLD), ("/E", 24, 0, TIMES_BOLD), ("3x5", 15, 9, TIMES_BOLD), (" > c?", 24, 0, TIMES_BOLD)], width=208)
    tree_node(draw, ll, [("core", 24, 0, TIMES_BOLD), (" + narrow?", 24, 0, TIMES_BOLD)], width=204)
    tree_node(draw, lr, [("mixed", 24, 0, TIMES_BOLD)], width=154)
    tree_node(draw, rr, [("broad?", 24, 0, TIMES_BOLD)], width=154)

    leaf_box(draw, leaves[0], "photon-like leaf", SPHENIX_BLUE, width=158)
    leaf_box(draw, leaves[1], "mixed leaf", PHOTON_DARK, width=132)
    leaf_box(draw, leaves[2], "background-like leaf", (205, 75, 55), width=182)
    draw.rounded_rectangle((x0 + 38, y1 - 84, x1 - 38, y1 - 30), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw_wrapped(draw, "Each leaf gives a small signal-like or background-like score.", (x0 + 58, y1 - 70), x1 - x0 - 116, font(TIMES_BOLD, 20), fill=BLUE, line_gap=1)


def draw_boosted_ensemble_stage(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    draw_pipeline_stage(
        base,
        box,
        "C",
        "Many shallow trees combine",
        "boosting adds weak learners",
        SPHENIX_BLUE,
    )
    origins = [(x0 + 70 + idx * 92, y0 + 246) for idx in range(5)]
    for idx, origin in enumerate(origins):
        draw_mini_tree(draw, origin, scale=0.55, alpha=220 - idx * 14)
        draw.text((origin[0] - 18, origin[1] + 92), f"T{idx + 1}", font=font(TIMES_BOLD, 18), fill=LIGHT_MUTED)
        if idx < 4:
            draw.text((origin[0] + 52, origin[1] + 90), "+", font=font(TIMES_BOLD, 24), fill=LIGHT_MUTED)

    bus_y = y0 + 500
    bus_x0 = x0 + 72
    bus_x1 = x1 - 72
    for origin in origins:
        draw.line((origin[0], origin[1] + 78, origin[0], bus_y), fill=(168, 184, 199, 175), width=3)
        draw.polygon([(origin[0], bus_y), (origin[0] - 7, bus_y - 12), (origin[0] + 7, bus_y - 12)], fill=(168, 184, 199, 175))
    draw.line((bus_x0, bus_y, bus_x1, bus_y), fill=(168, 184, 199, 205), width=4)
    draw.text((x0 + 124, bus_y + 18), "tree outputs add with learned weights", font=font(TIMES_ITALIC, 21), fill=MUTED)

    combiner = ((x0 + x1) // 2, y0 + 600)
    draw_arrow(draw, (combiner[0], bus_y + 40), (combiner[0], combiner[1] - 52), fill=(151, 169, 188), width=4)
    draw.ellipse((combiner[0] - 48, combiner[1] - 48, combiner[0] + 48, combiner[1] + 48), fill=(239, 246, 250, 255), outline=(*SPHENIX_BLUE, 235), width=4)
    draw.text((combiner[0] - 26, combiner[1] - 28), "Σ", font=font(TIMES_BOLD, 54), fill=BLUE)
    draw.text((combiner[0] - 64, combiner[1] + 62), "weighted sum", font=font(TIMES_ITALIC, 22), fill=MUTED)
    draw.rounded_rectangle((x0 + 46, y1 - 96, x1 - 46, y1 - 34), radius=8, fill=(255, 248, 229, 255), outline=(238, 220, 172, 255), width=2)
    draw_wrapped(draw, "many trees vote / add before one output is formed", (x0 + 68, y1 - 78), x1 - x0 - 136, font(TIMES_BOLD, 20), fill=BLUE, line_gap=2)


def draw_score_output_stage(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    draw_pipeline_stage(
        base,
        box,
        "D",
        "Photon-ID ranking score",
        "one output used for the tight-ID region",
        PHOTON,
    )
    draw.rounded_rectangle((x0 + 58, y0 + 196, x1 - 58, y0 + 292), radius=14, fill=(255, 248, 229, 255), outline=(*PHOTON_DARK, 230), width=3)
    draw.text((x0 + 92, y0 + 220), "BDT score", font=font(TIMES_BOLD, 36), fill=INK)
    draw.text((x0 + 94, y0 + 262), "photon-like ranking", font=font(TIMES_ITALIC, 23), fill=MUTED)

    bar_x0 = x0 + 72
    bar_x1 = x1 - 72
    bar_y = y0 + 414
    draw.text((bar_x0 - 6, bar_y - 56), "0", font=font(TIMES_BOLD, 24), fill=MUTED)
    draw.text((bar_x1 - 10, bar_y - 56), "1", font=font(TIMES_BOLD, 24), fill=MUTED)
    draw.rounded_rectangle((bar_x0, bar_y, bar_x1, bar_y + 42), radius=21, fill=(232, 236, 240, 255))
    for i in range(bar_x1 - bar_x0):
        t = i / max(1, bar_x1 - bar_x0 - 1)
        r = int(213 * (1 - t) + SPHENIX_BLUE[0] * t)
        g = int(88 * (1 - t) + SPHENIX_BLUE[1] * t)
        b = int(67 * (1 - t) + SPHENIX_BLUE[2] * t)
        draw.line((bar_x0 + i, bar_y, bar_x0 + i, bar_y + 42), fill=(r, g, b, 220), width=1)
    draw.text((bar_x0, bar_y + 66), "background-like", font=font(TIMES_BOLD, 22), fill=(180, 70, 56))
    draw.text((bar_x1 - 126, bar_y + 66), "photon-like", font=font(TIMES_BOLD, 22), fill=SPHENIX_BLUE)
    marker_x = bar_x0 + int(0.82 * (bar_x1 - bar_x0))
    draw.line((marker_x, bar_y - 24, marker_x, bar_y + 58), fill=(*INK, 255), width=4)
    draw.polygon([(marker_x, bar_y - 24), (marker_x - 13, bar_y - 46), (marker_x + 13, bar_y - 46)], fill=(*INK, 255))
    callout_label(draw, (marker_x - 60, bar_y - 94), "tight ID", BLUE, size=23)
    draw.rounded_rectangle((bar_x0 + 126, bar_y - 72, bar_x0 + 278, bar_y - 36), radius=7, fill=(255, 255, 255, 230), outline=(238, 220, 172, 255), width=2)
    draw.text((bar_x0 + 144, bar_y - 65), "non-tight", font=font(TIMES_BOLD, 17), fill=MUTED)

    note = (x0 + 50, y1 - 148, x1 - 50, y1 - 34)
    draw.rounded_rectangle(note, radius=9, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw_wrapped(draw, "The score ranks clusters; it is not a calibrated photon probability.", (note[0] + 24, note[1] + 26), note[2] - note[0] - 48, font(TIMES_BOLD, 22), fill=BLUE, line_gap=3)


def draw_isolation_cone(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    rounded_panel(draw, box, fill=PANEL)
    x0, y0, x1, y1 = box
    cx, cy = x0 + 260, y0 + 160
    draw.text((x0 + 34, y0 + 24), "Isolation requirement", font=font(TIMES_BOLD, 31), fill=BLUE)
    draw.pieslice((cx - 122, cy - 122, cx + 122, cy + 122), start=330, end=30, fill=(*PHOTON, 34), outline=(*PHOTON_DARK, 180), width=4)
    wave = feynman_points((x0 + 84, cy + 36), (cx + 132, cy - 14), 12, 5.4, 160)
    draw_polyline(draw, wave, (*PHOTON, 255), 8)
    draw_polyline(draw, wave, (*PHOTON_DARK, 210), 3)
    draw.ellipse((cx + 122, cy - 28, cx + 156, cy + 6), fill=(*PHOTON, 230), outline=(*PHOTON_DARK, 220), width=2)
    for angle in (70, 112, 246, 292):
        ex = cx + math.cos(math.radians(angle)) * 116
        ey = cy + math.sin(math.radians(angle)) * 116
        draw.ellipse((ex - 12, ey - 12, ex + 12, ey + 12), fill=(*TEAL_SOFT, 145))
    draw.text((x0 + 94, y1 - 72), "quiet cone keeps the prompt-photon tag interpretable", font=font(TIMES_ITALIC, 25), fill=MUTED)


def draw_isolation_population_panel(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    shadow(base, box)
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    outer_card_sidebar(draw, box, SPHENIX_BLUE)
    draw.text((x0 + 34, y0 + 26), "How to read the isolation distribution", font=font(TIMES_BOLD, 34), fill=INK)
    draw_wrapped(
        draw,
        "Isolation is the second selection axis after tight photon ID: quiet candidates sit near low reconstructed isolation energy.",
        (x0 + 34, y0 + 68),
        x1 - x0 - 68,
        font(TIMES_ITALIC, 24),
        fill=BLUE,
        line_gap=4,
    )

    rows = [
        {
            "label": "black points",
            "title": "tight-ID data",
            "body": "candidate distribution after preselection and tight photon ID",
            "accent": (24, 24, 24),
            "fill": (248, 249, 251),
            "marker": "points",
        },
        {
            "label": "red shade",
            "title": "non-tight-ID data",
            "body": "background-enriched sideband shape, normalized in the high-isolation tail",
            "accent": (207, 76, 62),
            "fill": (255, 239, 237),
            "marker": "shade",
        },
        {
            "label": "blue shade",
            "title": "tight-ID signal MC",
            "body": "prompt-photon template concentrated in the isolated, low-activity region",
            "accent": (84, 104, 224),
            "fill": (238, 241, 255),
            "marker": "shade",
        },
    ]
    y = y0 + 146
    for row in rows:
        ry0, ry1 = y, y + 112
        accent = row["accent"]
        draw.rounded_rectangle((x0 + 34, ry0, x1 - 34, ry1), radius=10, fill=(*row["fill"], 255), outline=(218, 226, 235, 255), width=2)
        if row["marker"] == "points":
            for px, py in ((x0 + 68, ry0 + 38), (x0 + 90, ry0 + 56), (x0 + 68, ry0 + 74)):
                draw.ellipse((px - 6, py - 6, px + 6, py + 6), fill=(*accent, 255))
                draw.line((px, py - 16, px, py + 16), fill=(*accent, 160), width=2)
        else:
            draw.rounded_rectangle((x0 + 54, ry0 + 31, x0 + 108, ry0 + 85), radius=5, fill=(*accent, 70), outline=(*accent, 220), width=3)
        draw.text((x0 + 134, ry0 + 20), row["label"], font=font(TIMES_BOLD, 24), fill=accent)
        draw.text((x0 + 318, ry0 + 18), row["title"], font=font(TIMES_BOLD, 28), fill=INK)
        draw_wrapped(draw, row["body"], (x0 + 318, ry0 + 56), x1 - x0 - 366, font(TIMES, 23), fill=MUTED, line_gap=3)
        y += 128


def draw_et_iso_formula(draw: ImageDraw.ImageDraw, xy: tuple[int, int]) -> None:
    x, y = xy
    base_font = font(TIMES_BOLD, 27)
    small_font = font(TIMES_BOLD, 16)

    def put(text: str, fnt: ImageFont.ImageFont, dx: int = 0, dy: int = 0, fill=INK) -> int:
        nonlocal x
        draw.text((x + dx, y + dy), text, font=fnt, fill=fill)
        w, _ = text_box(draw, text, fnt)
        x += w + dx
        return w

    put("E", base_font)
    draw.text((x - 1, y + 18), "T", font=small_font, fill=INK)
    draw.text((x + 10, y - 7), "iso,reco", font=small_font, fill=INK)
    x += 76
    put(" < 0.49 + 0.037 E", base_font)
    draw.text((x - 1, y + 18), "T", font=small_font, fill=INK)


def draw_isolation_cut_logic_panel(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    shadow(base, box)
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    outer_card_sidebar(draw, box, PHOTON)
    draw.text((x0 + 34, y0 + 24), "Isolation cut defines the second sideband axis", font=font(TIMES_BOLD, 31), fill=INK)

    left = (x0 + 34, y0 + 84, x0 + 452, y1 - 34)
    right = (x0 + 486, y0 + 84, x1 - 34, y1 - 34)

    draw.rounded_rectangle(left, radius=10, fill=(255, 248, 229, 255), outline=(*PHOTON, 220), width=2)
    draw.text((left[0] + 24, left[1] + 20), "isolated if", font=font(TIMES_BOLD, 27), fill=PHOTON_DARK)
    draw_et_iso_formula(draw, (left[0] + 24, left[1] + 66))
    draw_wrapped(
        draw,
        "threshold chosen for ~80% isolation efficiency",
        (left[0] + 24, left[1] + 124),
        left[2] - left[0] - 48,
        font(TIMES, 24),
        fill=MUTED,
        line_gap=4,
    )
    draw_wrapped(
        draw,
        "non-isolated sideband: cut + 0.8 GeV",
        (left[0] + 24, left[1] + 190),
        left[2] - left[0] - 48,
        font(TIMES_ITALIC, 23),
        fill=BLUE,
        line_gap=3,
    )

    draw.rounded_rectangle(right, radius=10, fill=(247, 250, 252, 255), outline=(*PANEL_EDGE, 255), width=2)
    draw.text((right[0] + 24, right[1] + 18), "ABCD purity regions", font=font(TIMES_BOLD, 28), fill=BLUE)
    cell_w, cell_h = 124, 78
    grid_w, grid_h = 2 * cell_w, 2 * cell_h
    gx0 = right[2] - 16 - grid_w
    gy0 = right[1] + 94

    column_font = font(TIMES_BOLD, 20)
    row_font = font(TIMES_BOLD, 20)
    for label, center_x, color in (
        ("isolated", gx0 + cell_w / 2, PHOTON_DARK),
        ("non-isolated", gx0 + cell_w + cell_w / 2, TEAL),
    ):
        lw, lh = text_box(draw, label, column_font)
        draw.text((center_x - lw / 2, gy0 - 36), label, font=column_font, fill=color)
    for label, center_y, color in (
        ("non-tight", gy0 + cell_h / 2, TEAL),
        ("tight", gy0 + cell_h + cell_h / 2, PHOTON_DARK),
    ):
        lw, lh = text_box(draw, label, row_font)
        draw.text((gx0 - 14 - lw, center_y - lh / 2 - 1), label, font=row_font, fill=color)
    cells = [
        ("C", (gx0, gy0), (239, 249, 250), TEAL),
        ("D", (gx0 + cell_w, gy0), (245, 247, 250), MUTED),
        ("A", (gx0, gy0 + cell_h), (255, 246, 221), PHOTON_DARK),
        ("B", (gx0 + cell_w, gy0 + cell_h), (238, 247, 252), SPHENIX_BLUE),
    ]
    for letter, (cx, cy), fill, accent in cells:
        draw.rounded_rectangle((cx, cy, cx + cell_w, cy + cell_h), radius=9, fill=(*fill, 255), outline=(*accent, 230), width=3)
        letter_font = font(TIMES_BOLD, 46)
        lw, lh = text_box(draw, letter, letter_font)
        draw.text((cx + (cell_w - lw) / 2, cy + (cell_h - lh) / 2 - 2), letter, font=letter_font, fill=accent)


def draw_analysis_flow(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    rounded_panel(draw, box, fill=PANEL)
    steps = [
        ("tight + iso", PHOTON),
        ("purity", SPHENIX_BLUE),
        ("efficiency", TEAL),
        ("unfolding", PHOTON_DARK),
        ("cross section", BLUE),
    ]
    x0, y0, x1, y1 = box
    x = x0 + 34
    y = y0 + 42
    for idx, (label, color) in enumerate(steps):
        w = 220 if idx < 4 else 260
        draw.rounded_rectangle((x, y, x + w, y + 70), radius=8, fill=(255, 255, 255, 255), outline=(*CARD_EDGE, 255), width=2)
        draw.rectangle((x, y, x + 10, y + 70), fill=(*color, 255))
        tw, _ = text_box(draw, label, font(TIMES_BOLD, 26))
        draw.text((x + (w - tw) / 2 + 4, y + 22), label, font=font(TIMES_BOLD, 26), fill=INK)
        if idx < len(steps) - 1:
            draw.line((x + w + 18, y + 35, x + w + 70, y + 35), fill=(158, 178, 194, 255), width=3)
            draw.polygon([(x + w + 70, y + 35), (x + w + 56, y + 26), (x + w + 56, y + 44)], fill=(158, 178, 194, 255))
        x += w + 78


def draw_closing_visual(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    rounded_panel(draw, box, fill=PANEL)
    x0, y0, x1, y1 = box
    draw.text((x0 + 36, y0 + 26), "Baseline to future hard probes", font=font(TIMES_BOLD, 31), fill=BLUE)
    labels = [
        ("p+p photons", PHOTON, TIMES_BOLD, INK),
        ("gamma-jet", SPHENIX_BLUE, TIMES, MUTED),
        ("heavy-ion", TEAL, TIMES, MUTED),
    ]
    cy = y0 + 116
    centers = [x0 + 142, x0 + 472, x0 + 802]
    for idx, (label, color, font_path, text_fill) in enumerate(labels):
        cx = centers[idx]
        draw.ellipse((cx - 22, cy - 22, cx + 22, cy + 22), fill=(*color, 230))
        tw, _ = text_box(draw, label, font(font_path, 28))
        draw.text((cx - tw / 2, cy + 38), label, font=font(font_path, 28), fill=text_fill)
        if idx < len(labels) - 1:
            x_start = cx + 38
            x_end = centers[idx + 1] - 42
            draw.line((x_start, cy, x_end, cy), fill=(165, 184, 198), width=3)
            draw.polygon([(x_end, cy), (x_end - 15, cy - 9), (x_end - 15, cy + 9)], fill=(165, 184, 198))


def place_ian_plot_card(
    base: Image.Image,
    key: str,
    box: tuple[int, int, int, int],
    title: str,
    subtitle: str,
) -> tuple[int, int, int, int]:
    draw = ImageDraw.Draw(base, "RGBA")
    shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    draw.text((box[0] + 28, box[1] + 20), title, font=font(TIMES_BOLD, 30), fill=INK)
    draw_wrapped(draw, subtitle, (box[0] + 28, box[1] + 62), box[2] - box[0] - 56, font(TIMES_ITALIC, 22), fill=MUTED, line_gap=3)
    img = Image.open(ian_figure_path(key)).convert("RGBA")
    return paste_fit(base, img, (box[0] + 24, box[1] + 104, box[2] - 24, box[3] - 24), anchor="center")


def draw_npb_training_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    x0, y0, x1, y1 = box
    draw.text((x0 + 30, y0 + 24), "What the NPB score is trained to reject", font=font(TIMES_BOLD, 30), fill=INK)
    draw.line((x0 + 30, y0 + 68, x1 - 30, y0 + 68), fill=(222, 229, 236), width=2)

    rows = [
        ("signal label", "PYTHIA8 physics clusters", PHOTON),
        ("NPB label", "data clusters with anomalous timing and no recoil jet", (205, 75, 55)),
        ("inputs", "shower shape + cluster kinematics; timing validates the tag", SPHENIX_BLUE),
    ]
    y = y0 + 98
    for heading, body, accent in rows:
        draw.rounded_rectangle((x0 + 30, y, x1 - 30, y + 92), radius=10, fill=(247, 250, 252, 255), outline=(219, 228, 236, 255), width=2)
        draw.rounded_rectangle((x0 + 30, y, x0 + 42, y + 92), radius=6, fill=(*accent, 255))
        draw.text((x0 + 62, y + 16), heading, font=font(TIMES_BOLD, 24), fill=BLUE)
        draw_wrapped(draw, body, (x0 + 62, y + 47), x1 - x0 - 112, font(TIMES, 22), fill=MUTED, line_gap=2)
        y += 112

    draw.rounded_rectangle((x0 + 30, y1 - 98, x1 - 30, y1 - 24), radius=8, fill=(255, 248, 229, 255), outline=(238, 220, 172, 255), width=2)
    draw_wrapped(
        draw,
        "Nominal cut: NPB > 0.5, chosen to be mild for physics clusters.",
        (x0 + 56, y1 - 78),
        x1 - x0 - 112,
        font(TIMES_BOLD, 23),
        fill=BLUE,
        line_gap=2,
    )


def draw_badge(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], heading: str, body: str, accent: tuple[int, int, int]) -> None:
    draw.rounded_rectangle(box, radius=9, fill=(255, 255, 255, 245), outline=(*accent, 190), width=2)
    draw.rounded_rectangle((box[0], box[1], box[0] + 10, box[3]), radius=5, fill=(*accent, 255))
    draw.text((box[0] + 26, box[1] + 10), heading, font=font(TIMES_BOLD, 23), fill=BLUE)
    draw_wrapped(draw, body, (box[0] + 26, box[1] + 42), box[2] - box[0] - 46, font(TIMES, 21), fill=MUTED, line_gap=2)


def draw_npb_time_evidence_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    draw.text((x0 + 30, y0 + 22), "Non-collision signature", font=font(TIMES_BOLD, 32), fill=INK)
    draw_wrapped(
        draw,
        "Timing validates what the NPB score is seeing.",
        (x0 + 30, y0 + 64),
        x1 - x0 - 60,
        font(TIMES_ITALIC, 23),
        fill=MUTED,
        line_gap=3,
    )
    plot = Image.open(ian_figure_path("fig15_npb_score_vs_time")).convert("RGBA")
    plot = plot.crop((0, 0, plot.width, int(plot.height * 0.90)))
    plot = crop_visible(plot, white_threshold=252, pad=12)
    plot_box = (x0 + 34, y0 + 104, x0 + 752, y1 - 28)
    placed = paste_fit(base, plot, plot_box, anchor="center")
    px0, py0, px1, py1 = placed

    plot_w = px1 - px0
    plot_h = py1 - py0
    # Axis-frame coordinates measured from the cropped IAN heatmap asset:
    # left y-axis, right heatmap frame before the color bar, top frame, bottom frame.
    axis_left = px0 + round(plot_w * 110 / 849)
    axis_right = px0 + round(plot_w * 722 / 849)
    axis_top = py0 + round(plot_h * 15 / 682)
    axis_bottom = py0 + round(plot_h * 590 / 682)
    gate_y = round(axis_top + 0.5 * (axis_bottom - axis_top))
    gate_x0 = axis_left + 8
    gate_x1 = axis_right
    draw.line((gate_x0, gate_y, gate_x1, gate_y), fill=(*PHOTON_DARK, 210), width=3)
    gate_font = font(TIMES_BOLD, 22)
    draw.text(
        (gate_x0 + 18, gate_y - 30),
        "pre: NPB > 0.5",
        font=gate_font,
        fill=PHOTON_DARK,
        stroke_width=2,
        stroke_fill=(255, 255, 255, 245),
    )

    callout_x0 = x0 + 778
    callout_x1 = x1 - 30
    draw_badge(
        draw,
        (callout_x0, y0 + 142, callout_x1, y0 + 248),
        "high NPB score",
        "in-time, physics-like clusters",
        SPHENIX_BLUE,
    )
    draw_arrow(draw, (callout_x0 - 12, y0 + 196), (px0 + int(0.70 * (px1 - px0)), py0 + int(0.30 * (py1 - py0))), fill=SPHENIX_BLUE, width=3)
    draw_badge(
        draw,
        (callout_x0, y0 + 282, callout_x1, y0 + 388),
        "low NPB score",
        "out-of-time / non-collision component",
        (175, 70, 56),
    )
    draw_arrow(draw, (callout_x0 - 12, y0 + 336), (px0 + int(0.36 * (px1 - px0)), py0 + int(0.74 * (py1 - py0))), fill=(175, 70, 56), width=3)
    draw_badge(
        draw,
        (callout_x0, y0 + 422, callout_x1, y0 + 528),
        "nominal gate",
        "NPB > 0.5 before photon ID",
        PHOTON_DARK,
    )


def draw_npb_threshold_evidence_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    draw.text((x0 + 30, y0 + 22), "Mild gate, high retention", font=font(TIMES_BOLD, 32), fill=INK)
    draw_wrapped(
        draw,
        "The threshold is chosen to clean the pool without sculpting it.",
        (x0 + 30, y0 + 64),
        x1 - x0 - 60,
        font(TIMES_ITALIC, 23),
        fill=MUTED,
        line_gap=3,
    )
    plot = Image.open(ian_figure_path("fig17_npb_threshold_scan")).convert("RGBA")
    plot = crop_visible(plot, white_threshold=252, pad=12)
    placed = paste_fit(base, plot, (x0 + 42, y0 + 116, x0 + 724, y1 - 42), anchor="center")
    px0, py0, px1, py1 = placed

    badge_x0 = x0 + 766
    badge_x1 = x1 - 42
    badge_h = 124
    draw_badge(
        draw,
        (badge_x0, y0 + 206, badge_x1, y0 + 206 + badge_h),
        "about 99% retained",
        "physics-cluster signal retention stays high",
        SPHENIX_BLUE,
    )
    draw_badge(
        draw,
        (badge_x0, y0 + 360, badge_x1, y0 + 360 + badge_h),
        "about 99% clean",
        "cluster-pool hygiene, not photon purity",
        PHOTON_DARK,
    )


def draw_npb_sideband_flow_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    draw.text((x0 + 30, y0 + 22), "Why this comes before photon ID", font=font(TIMES_BOLD, 32), fill=INK)
    draw_wrapped(
        draw,
        "Preselection defines the candidate pool before tight and non-tight photon-ID regions are interpreted as related sidebands.",
        (x0 + 30, y0 + 64),
        x1 - x0 - 60,
        font(TIMES_ITALIC, 23),
        fill=MUTED,
        line_gap=3,
    )

    steps = [
        ("raw EMCal\nclusters", "contains collision-like clusters plus artifacts", (205, 75, 55)),
        ("NPB +\npreselection", "remove out-of-time / pathological clusters", SPHENIX_BLUE),
        ("physics-candidate\nparent sample", "cleaner population for photon-ID sidebands", PHOTON_DARK),
        ("tight ID\nnon-tight ID", "same parent sample, different purity", TEAL),
    ]
    step_w = 390
    step_h = 118
    gap = 48
    total_w = 4 * step_w + 3 * gap
    sx = x0 + (x1 - x0 - total_w) // 2
    sy = y0 + 126
    centers = []
    for idx, (heading, body, accent) in enumerate(steps):
        bx0 = sx + idx * (step_w + gap)
        bx1 = bx0 + step_w
        box_i = (bx0, sy, bx1, sy + step_h)
        draw.rounded_rectangle(box_i, radius=13, fill=(247, 250, 252, 255), outline=(*accent, 220), width=3)
        draw.rounded_rectangle((bx0, sy, bx0 + 12, sy + step_h), radius=6, fill=(*accent, 255))
        lines = heading.split("\n")
        draw.text((bx0 + 34, sy + 18), lines[0], font=font(TIMES_BOLD, 27), fill=INK)
        if len(lines) > 1:
            draw.text((bx0 + 34, sy + 50), lines[1], font=font(TIMES_BOLD, 27), fill=INK)
        draw_wrapped(draw, body, (bx0 + 34, sy + 82), step_w - 62, font(TIMES, 19), fill=MUTED, line_gap=1)
        centers.append((bx1, sy + step_h // 2, bx0, sy + step_h // 2))
        if idx < len(steps) - 1:
            draw_arrow(draw, (bx1 + 12, sy + step_h // 2), (bx1 + gap - 12, sy + step_h // 2), fill=(158, 178, 194), width=4)

    note = (x0 + 36, y1 - 78, x1 - 36, y1 - 24)
    draw.rounded_rectangle(note, radius=8, fill=(255, 248, 229, 255), outline=(238, 220, 172, 255), width=2)
    note_x = note[0] + 24
    note_y = note[1] + 14
    lead = "Sideband method logic:"
    lead_font = font(TIMES_BOLD, 25)
    body_font = font(TIMES, 25)
    draw.text((note_x, note_y), lead, font=lead_font, fill=BLUE)
    lead_w, _ = text_box(draw, lead, lead_font)
    draw_wrapped(
        draw,
        "non-tight should be a background-enriched version of the same photon-like parent population, not a disconnected junk population.",
        (note_x + lead_w + 8, note_y),
        note[2] - (note_x + lead_w + 8) - 24,
        body_font,
        fill=BLUE,
        line_gap=2,
    )


def draw_preselection_species_diagram(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    draw.text((x0 + 34, y0 + 24), "Why this cut belongs before tight photon ID", font=font(TIMES_BOLD, 32), fill=INK)
    draw.text((x0 + 34, y0 + 68), "Preselection is species control, not purity maximization.", font=font(TIMES_ITALIC, 27), fill=BLUE)

    lane_y = y0 + 130
    rejected = (x0 + 52, lane_y, x0 + 352, lane_y + 106)
    gate = (x0 + 430, lane_y - 12, x0 + 820, lane_y + 118)
    parent = (x0 + 930, lane_y - 16, x0 + 1378, lane_y + 122)
    tight = (x0 + 1538, lane_y - 24, x0 + 1856, lane_y + 60)
    nontight = (x0 + 1538, lane_y + 86, x0 + 1856, lane_y + 170)
    sideband = (x0 + 1930, lane_y - 2, x1 - 54, lane_y + 148)

    draw.rounded_rectangle(rejected, radius=12, fill=(252, 238, 235, 255), outline=(232, 190, 184, 255), width=2)
    draw.line((rejected[0] + 26, rejected[1] + 28, rejected[0] + 72, rejected[1] + 74), fill=(205, 75, 55, 255), width=5)
    draw.line((rejected[0] + 72, rejected[1] + 28, rejected[0] + 26, rejected[1] + 74), fill=(205, 75, 55, 255), width=5)
    draw_wrapped(draw, "non-collision or pathological clusters", (rejected[0] + 96, rejected[1] + 24), rejected[2] - rejected[0] - 122, font(TIMES_BOLD, 23), fill=(118, 49, 44), line_gap=3)

    draw_arrow(draw, (rejected[2] + 20, rejected[1] + 52), (gate[0] - 18, gate[1] + 68), fill=(168, 184, 199), width=4)
    draw.rounded_rectangle(gate, radius=14, fill=(239, 246, 250, 255), outline=(*SPHENIX_BLUE, 230), width=3)
    draw.text((gate[0] + 28, gate[1] + 24), "preselection gate", font=font(TIMES_BOLD, 27), fill=BLUE)
    draw_wrapped(draw, "NPB > 0.5 plus compact-shape sanity cuts", (gate[0] + 28, gate[1] + 62), gate[2] - gate[0] - 56, font(TIMES, 22), fill=MUTED, line_gap=2)

    draw_arrow(draw, (gate[2] + 18, gate[1] + 68), (parent[0] - 18, parent[1] + 70), fill=(168, 184, 199), width=4)
    draw.rounded_rectangle(parent, radius=14, fill=(255, 249, 232, 255), outline=(*PHOTON_DARK, 230), width=3)
    draw.text((parent[0] + 34, parent[1] + 28), "same photon-like", font=font(TIMES_BOLD, 31), fill=INK)
    draw.text((parent[0] + 34, parent[1] + 69), "parent population", font=font(TIMES_BOLD, 31), fill=INK)

    split_start = (parent[2] + 22, parent[1] + 70)
    draw.line((split_start[0], split_start[1], tight[0] - 36, tight[1] + 42), fill=(168, 184, 199, 255), width=4)
    draw.polygon([(tight[0] - 36, tight[1] + 42), (tight[0] - 52, tight[1] + 33), (tight[0] - 52, tight[1] + 51)], fill=(168, 184, 199, 255))
    draw.line((split_start[0], split_start[1], nontight[0] - 36, nontight[1] + 42), fill=(168, 184, 199, 255), width=4)
    draw.polygon([(nontight[0] - 36, nontight[1] + 42), (nontight[0] - 52, nontight[1] + 33), (nontight[0] - 52, nontight[1] + 51)], fill=(168, 184, 199, 255))

    draw.rounded_rectangle(tight, radius=12, fill=(245, 250, 255, 255), outline=(*SPHENIX_BLUE, 230), width=3)
    draw.text((tight[0] + 26, tight[1] + 18), "tight ID", font=font(TIMES_BOLD, 27), fill=BLUE)
    draw.text((tight[0] + 26, tight[1] + 50), "higher purity", font=font(TIMES_ITALIC, 21), fill=MUTED)

    draw.rounded_rectangle(nontight, radius=12, fill=(248, 251, 250, 255), outline=(*TEAL, 220), width=3)
    draw.text((nontight[0] + 26, nontight[1] + 18), "non-tight ID", font=font(TIMES_BOLD, 27), fill=TEAL)
    draw.text((nontight[0] + 26, nontight[1] + 50), "background-enriched", font=font(TIMES_ITALIC, 21), fill=MUTED)

    draw.rounded_rectangle(sideband, radius=12, fill=(247, 250, 252, 255), outline=(219, 228, 236, 255), width=2)
    draw_wrapped(
        draw,
        "Sideband logic is strongest when tight and non-tight are related populations with different purity.",
        (sideband[0] + 28, sideband[1] + 28),
        sideband[2] - sideband[0] - 56,
        font(TIMES_BOLD, 24),
        fill=INK,
        line_gap=4,
    )


def draw_public_ncb_photonid_flow(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    draw.text((x0, y0), "Mechanism: clean, prepare, then classify", font=font(TIMES_BOLD, 37), fill=BLUE)
    draw_wrapped(
        draw,
        "NCB controls what enters the candidate pool; photon-ID BDT defines the axis used for tight and non-tight sidebands.",
        (x0, y0 + 42),
        x1 - x0,
        font(TIMES_ITALIC, 25),
        fill=MUTED,
        line_gap=2,
    )

    stages = [
        (
            "1",
            "Clean the pool",
            "NCB cleaning BDT",
            "removes non-collisional and pathological clusters",
            (232, 244, 252),
            SPHENIX_BLUE,
            None,
        ),
        (
            "2",
            "Prepare feature space",
            "preselection shape gates",
            "keep the shower-shape variables interpretable and reduce ID-isolation correlation",
            (255, 247, 224),
            PHOTON_DARK,
            None,
        ),
        (
            "3",
            "Define the ID axis",
            "photon-ID BDT score",
            "ranks prompt-like showers, then splits into tight and non-tight regions",
            (236, 248, 245),
            TEAL,
            "split",
        ),
    ]
    top = y0 + 90
    card_h = y1 - top
    gap = 58
    step_w = (x1 - x0 - gap * 2) // 3
    sx = x0
    for i, (num, eyebrow, title, body, fill, accent, extra) in enumerate(stages):
        bx = sx + i * (step_w + gap)
        card = (bx, top, bx + step_w, top + card_h)
        draw.rounded_rectangle(card, radius=12, fill=(*fill, 255), outline=(*CARD_EDGE, 255), width=2)
        draw.rectangle((card[0], card[1], card[0] + 14, card[3]), fill=(*accent, 255))
        badge = (card[0] + 34, card[1] + 24, card[0] + 90, card[1] + 80)
        draw.ellipse(badge, fill=(255, 255, 255, 255), outline=(*accent, 255), width=3)
        num_font = font(TIMES_BOLD, 31)
        nb = draw.textbbox((0, 0), num, font=num_font)
        nx = (badge[0] + badge[2] - nb[0] - nb[2]) / 2
        ny = (badge[1] + badge[3] - nb[1] - nb[3]) / 2
        draw.text((nx, ny), num, font=num_font, fill=accent)
        draw.text((card[0] + 110, card[1] + 25), eyebrow, font=font(TIMES_BOLD, 29), fill=accent)
        draw.text((card[0] + 34, card[1] + 84), title, font=font(TIMES_BOLD, 33), fill=INK)
        if extra == "split":
            draw.text((card[0] + 34, card[1] + 126), "score creates:", font=font(TIMES, 24), fill=MUTED)
            pill_y = card[1] + 120
            pill_w = 190
            tight = (card[0] + 214, pill_y, card[0] + 214 + pill_w, pill_y + 44)
            nontight = (tight[2] + 20, pill_y, tight[2] + 20 + pill_w + 50, pill_y + 44)
            draw.rounded_rectangle(tight, radius=8, fill=(255, 255, 255, 255), outline=(*TEAL, 235), width=3)
            draw.rounded_rectangle(nontight, radius=8, fill=(255, 255, 255, 255), outline=(92, 112, 142, 235), width=3)
            draw.text((tight[0] + 22, tight[1] + 8), "tight ID", font=font(TIMES_BOLD, 24), fill=TEAL)
            draw.text((nontight[0] + 22, nontight[1] + 8), "non-tight ID", font=font(TIMES_BOLD, 24), fill=(92, 112, 142))
        else:
            draw_wrapped(draw, body, (card[0] + 34, card[1] + 126), step_w - 68, font(TIMES, 24), fill=MUTED, line_gap=2)
        if i < len(stages) - 1:
            mid_y = top + card_h // 2
            draw_arrow(draw, (card[2] + 14, mid_y), (card[2] + gap - 14, mid_y), fill=(142, 161, 181), width=6)

def draw_fig1_public_story_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    shadow(base, box, radius=12)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    x0, y0, x1, y1 = box
    draw.text((x0 + 30, y0 + 20), "After NCB + preselection: shower-shape space", font=font(TIMES_BOLD, 30), fill=BLUE)
    draw_wrapped(
        draw,
        "The cleaned/preselected sample is where shower-shape variables become useful photon-ID features.",
        (x0 + 30, y0 + 60),
        x1 - x0 - 60,
        font(TIMES_ITALIC, 25),
        fill=MUTED,
        line_gap=1,
    )
    plot = Image.open(figure_path("fig1_shower_shape")).convert("RGBA")
    plot = crop_visible(plot, white_threshold=252, pad=8)
    paste_fit(base, plot, (x0 + 20, y0 + 116, x1 - 20, y1 - 68), anchor="center")

    note_y = y1 - 48
    icon_x = x0 + 36
    icon_y = note_y - 2
    draw.rounded_rectangle((icon_x, icon_y, icon_x + 30, icon_y + 30), radius=4, fill=(255, 255, 255, 0), outline=(*SPHENIX_BLUE, 255), width=2)
    draw.line((icon_x + 8, icon_y + 10, icon_x + 22, icon_y + 10), fill=(*SPHENIX_BLUE, 210), width=2)
    draw.line((icon_x + 8, icon_y + 17, icon_x + 21, icon_y + 17), fill=(*SPHENIX_BLUE, 180), width=2)
    draw.line((icon_x + 8, icon_y + 24, icon_x + 17, icon_y + 24), fill=(*SPHENIX_BLUE, 150), width=2)
    note_font = font(TIMES_BOLD, 22)
    body_font = font(TIMES, 22)
    note_x = icon_x + 42
    draw.text((note_x, note_y), "Note:", font=note_font, fill=BLUE)
    note_w, _ = text_box(draw, "Note:", note_font)
    draw_wrapped(
        draw,
        "data still background-rich, but the shape variables now carry interpretable prompt-vs-background separation.",
        (note_x + note_w + 12, note_y + 1),
        x1 - (note_x + note_w + 12) - 36,
        body_font,
        fill=MUTED,
        line_gap=0,
    )


def draw_fig2_public_story_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    shadow(base, box, radius=12)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    x0, y0, x1, y1 = box
    draw.text((x0 + 30, y0 + 20), "Photon-ID BDT: one score for prompt-like showers", font=font(TIMES_BOLD, 29), fill=BLUE)
    draw_wrapped(
        draw,
        "The score separates prompt-like showers from inclusive-jet background after NCB and preselection.",
        (x0 + 30, y0 + 58),
        x1 - x0 - 60,
        font(TIMES_ITALIC, 25),
        fill=MUTED,
        line_gap=1,
    )
    plot = Image.open(figure_path("fig2_bdt_score")).convert("RGBA")
    plot = crop_visible(plot, white_threshold=252, pad=10)
    paste_fit(base, plot, (x0 + 42, y0 + 118, x1 - 42, y1 - 36), anchor="center")


def slide07_public_ncb_preselection_photonid_bdt() -> tuple[Path, Path]:
    img = base_slide(
        "NCB cleaning prepares the photon-ID BDT",
        "First clean the candidate pool; then use shower shapes to define the tight and non-tight photon-ID axis.",
    )
    add_top_right_sphenix_logo_like_slide2(img)
    draw_fig1_public_story_card(img, (132, 306, 1352, 984))
    draw_fig2_public_story_card(img, (1392, 306, 2390, 984))
    draw_public_ncb_photonid_flow(img, (92, 1015, 2480, 1295))
    draw_hp2026_identity_footer(img)

    png = OUTPUT / "hp2026_slide07_public_ncb_preselection_photonid_bdt.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        7,
        "NCB cleaning prepares the photon-ID BDT",
        "Now I want to make the analysis sequence very explicit, because there are two different BDT ideas that can sound similar if I say them too quickly. The paper calls the first one NCB, for non-collisional background. That is the cleaning step. It removes pathological clusters from beam-pipe interactions and other non-collision-like backgrounds before I ever define the photon-ID sideband regions.\n\nOn the left, Fig. 1 shows the paper-supported view after NCB and preselection. These are shower-shape variables: w eta on the left and E3 by 2 over E3 by 5 on the right. At this stage the data are still background-rich, which is why they follow the inclusive-jet MC more closely than the prompt-photon MC. But the important point is that the variables are now interpretable: prompt photons tend to occupy the narrow, compact-shower regions, while background clusters populate broader or less compact shapes.\n\nThe mechanism strip at the bottom is the clean separation of jobs. First, NCB cleaning controls what enters the candidate pool. Second, preselection keeps the shower-shape feature space usable and reduces the ID-isolation correlation. Third, the photon-ID BDT defines the actual identification axis, splitting the selected candidate population into tight and non-tight ID regions.\n\nOn the right, Fig. 2 shows that second BDT. This is the photon-identification score after NCB and preselection. The prompt-photon signal MC moves toward high score, while the inclusive-jet background peaks toward low score. That score is what defines the tight-ID and non-tight-ID regions. So the clean takeaway is: first I clean the candidate pool, then I use shower shapes to build the ID axis, and only after that do isolation and sidebands enter the purity measurement.",
        stem="hp2026_slide07_public_ncb_preselection_photonid_bdt_script",
    )
    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "google_slides_mutation": False,
        "candidate_role": "public-safe replacement candidate for current Slides 7 and 8",
        "png": str(png.relative_to(ROOT)),
        "speaker_script": str(script.relative_to(ROOT)),
        "size": [W, H],
        "mode": "RGB",
        "source_policy": {
            "plot_sources": [
                {
                    "key": "fig1_shower_shape",
                    "source_pdf": str(PAPER.relative_to(ROOT)),
                    "source_type": "paper plot",
                    "public_status": "current-paper-source_pending_public_label_check",
                },
                {
                    "key": "fig2_bdt_score",
                    "source_pdf": str(PAPER.relative_to(ROOT)),
                    "source_type": "paper plot",
                    "public_status": "current-paper-source_pending_public_label_check",
                },
            ],
            "generated_graphics_role": "explanatory scaffolding only; no generated physics-result plots",
            "excluded_sources": "No IAN plots, backup-slide screenshot crops, or internal ROOT/data-generated physics plots are used.",
            "public_readiness_caveat": "The embedded paper plots still show sPHENIX Internal in the current PDF and must be replaced by Preliminary/approved public-label versions or documented approval before final public use.",
        },
    }
    manifest_path = OUTPUT / "hp2026_slide07_public_ncb_preselection_photonid_bdt_manifest.json"
    with manifest_path.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return png, script


def draw_public_bdt_score_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    draw.text((x0 + 32, y0 + 24), "Photon-ID score after preselection", font=font(TIMES_BOLD, 34), fill=INK)
    draw_wrapped(
        draw,
        "Public paper figure: prompt-photon MC populates the high-score region; inclusive-jet MC peaks toward low score.",
        (x0 + 32, y0 + 68),
        x1 - x0 - 64,
        font(TIMES_ITALIC, 24),
        fill=MUTED,
        line_gap=4,
    )
    plot = Image.open(figure_path("fig2_bdt_score")).convert("RGBA")
    placed = paste_fit(base, plot, (x0 + 36, y0 + 130, x1 - 36, y1 - 112), anchor="center")

    px0, py0, px1, py1 = placed
    # The plotted threshold is ET dependent; this is a qualitative high-score
    # guide, not an extra derived result curve.
    guide_x = px0 + int(0.79 * (px1 - px0))
    draw.line((guide_x, py0 + 210, guide_x, py1 - 108), fill=(*SPHENIX_BLUE, 230), width=4)
    draw.polygon([(guide_x, py0 + 202), (guide_x - 12, py0 + 180), (guide_x + 12, py0 + 180)], fill=(*SPHENIX_BLUE, 230))
    callout_label(draw, (guide_x - 96, py0 + 26), "tight-ID side", BLUE, bg=(255, 255, 255), size=23)

    note = (x0 + 36, y1 - 86, x1 - 36, y1 - 26)
    draw.rounded_rectangle(note, radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw_wrapped(
        draw,
        "The score is a reproducible ID axis, not the physics result itself.",
        (note[0] + 22, note[1] + 14),
        note[2] - note[0] - 44,
        font(TIMES_BOLD, 25),
        fill=BLUE,
        line_gap=2,
    )


def draw_tight_id_before_after_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    draw.text((x0 + 30, y0 + 24), "Selection flow in the compact-core observable", font=font(TIMES_BOLD, 32), fill=INK)
    draw_wrapped(
        draw,
        "The same E11/E33 shape check becomes progressively more prompt-like as selection is applied.",
        (x0 + 30, y0 + 66),
        x1 - x0 - 60,
        font(TIMES_ITALIC, 23),
        fill=MUTED,
        line_gap=4,
    )

    paths = slide08_e11e33_flow_paths()
    gutter = 28
    arrow_w = 36
    top = y0 + 118
    label_h = 42
    plot_top = top + label_h + 12
    plot_bottom = y1 - 78
    col_w = (x1 - x0 - 2 * gutter - 2 * arrow_w - 2 * 18) // 3
    columns = [
        (x0 + gutter, "no preselection", "no_preselection", PHOTON_DARK),
        (x0 + gutter + col_w + arrow_w + 18, "preselection", "preselection", TEAL),
        (x0 + gutter + 2 * (col_w + arrow_w + 18), "tight ID", "tight_id", SPHENIX_BLUE),
    ]
    for idx, (cx0, label, key, accent) in enumerate(columns):
        cx1 = cx0 + col_w
        draw.rounded_rectangle((cx0, top, cx1, top + label_h), radius=8, fill=(247, 250, 252, 255), outline=(*accent, 190), width=2)
        tw, _ = text_box(draw, label, font(TIMES_BOLD, 22))
        draw.text((cx0 + (cx1 - cx0 - tw) / 2, top + 11), label, font=font(TIMES_BOLD, 22), fill=accent)
        draw.rounded_rectangle((cx0, plot_top, cx1, plot_bottom), radius=8, fill=(255, 255, 255, 255), outline=(224, 231, 238, 255), width=2)
        if key in paths:
            plot = Image.open(paths[key]).convert("RGBA")
            paste_fit(base, plot, (cx0 + 10, plot_top + 10, cx1 - 10, plot_bottom - 10), anchor="center")
        else:
            draw_wrapped(draw, "backup plot crop unavailable", (cx0 + 20, plot_top + 90), col_w - 40, font(TIMES_ITALIC, 22), fill=MUTED)
        if idx < 2:
            ax = cx1 + 9
            ay = (plot_top + plot_bottom) // 2
            draw.line((ax, ay, ax + arrow_w - 12, ay), fill=(154, 168, 184, 205), width=4)
            draw.polygon([(ax + arrow_w - 12, ay), (ax + arrow_w - 26, ay - 10), (ax + arrow_w - 26, ay + 10)], fill=(154, 168, 184, 205))

    draw.rounded_rectangle((x0 + 30, y1 - 60, x1 - 30, y1 - 26), radius=7, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=1)
    bottom = "The BDT selection is visible as a shape-distribution transformation, not only as a score cut."
    tw, _ = text_box(draw, bottom, font(TIMES_BOLD, 20))
    draw.text((x0 + (x1 - x0 - tw) / 2, y1 - 54), bottom, font=font(TIMES_BOLD, 20), fill=BLUE)


def draw_tight_id_training_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    draw.text((x0 + 30, y0 + 24), "What the BDT learns", font=font(TIMES_BOLD, 31), fill=INK)
    rows = [
        ("signal", "truth-matched direct or fragmentation photons", PHOTON),
        ("background", "complement of prompt definition in inclusive-jet PYTHIA MC", (205, 75, 55)),
        ("inputs", "shower shapes plus cluster kinematics and vertex", SPHENIX_BLUE),
    ]
    y = y0 + 84
    for heading, body, accent in rows:
        draw.rounded_rectangle((x0 + 30, y, x1 - 30, y + 82), radius=9, fill=(247, 250, 252, 255), outline=(219, 228, 236, 255), width=2)
        draw.rounded_rectangle((x0 + 30, y, x0 + 42, y + 82), radius=5, fill=(*accent, 255))
        draw.text((x0 + 60, y + 14), heading, font=font(TIMES_BOLD, 23), fill=BLUE)
        draw_wrapped(draw, body, (x0 + 60, y + 43), x1 - x0 - 106, font(TIMES, 21), fill=MUTED, line_gap=2)
        y += 98


def draw_tight_id_threshold_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0, x1, y1 = box
    shadow(base, box)
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*PANEL_EDGE, 255), width=2)
    draw.text((x0 + 30, y0 + 24), "Tight-ID working point", font=font(TIMES_BOLD, 31), fill=INK)
    eq = (x0 + 30, y0 + 82, x1 - 30, y0 + 154)
    draw.rounded_rectangle(eq, radius=9, fill=(255, 248, 229, 255), outline=(238, 220, 172, 255), width=2)
    draw_formula_run(
        draw,
        (eq[0] + 24, eq[1] + 19),
        [
            ("ID BDT > 0.8156 - 0.00156 E", 29, 0, TIMES_BOLD),
            ("T", 18, 12, TIMES_BOLD),
        ],
        fill=INK,
    )
    draw_wrapped(
        draw,
        "chosen for about 80% identification efficiency relative to preselection",
        (x0 + 38, y0 + 178),
        x1 - x0 - 76,
        font(TIMES_BOLD, 24),
        fill=BLUE,
        line_gap=5,
    )
    draw_wrapped(
        draw,
        "The tight region becomes the high-purity ID slice used with isolation in the sideband method.",
        (x0 + 38, y0 + 258),
        x1 - x0 - 76,
        font(TIMES, 23),
        fill=MUTED,
        line_gap=6,
    )


def slugify(value: str) -> str:
    cleaned = []
    for ch in value.lower():
        if ch.isalnum():
            cleaned.append(ch)
        elif ch in {" ", "-", "_"}:
            cleaned.append("_")
    slug = "".join(cleaned)
    while "__" in slug:
        slug = slug.replace("__", "_")
    return slug.strip("_")[:72]


def save_script(slide_no: int, title: str, body: str, stem: str | None = None) -> Path:
    SCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    if stem is None:
        stem = f"hp2026_slide{slide_no:02d}_{slugify(title)}_script"
    path = SCRIPT_DIR / f"{stem}.md"
    path.write_text(f"# HP2026 Slide {slide_no} Script - {title}\n\n{body.strip()}\n", encoding="utf-8")
    return path


def write_contact_sheet(outputs: list[tuple[int, Path, Path]]) -> Path:
    thumb_w, thumb_h = 640, 360
    label_h = 52
    rows = math.ceil(len(outputs) / 2)
    sheet = Image.new("RGB", (2 * thumb_w, rows * (thumb_h + label_h)), "white")
    draw = ImageDraw.Draw(sheet)
    label_font = font(TIMES, 24)
    for idx, (_, png, _) in enumerate(outputs):
        col = idx % 2
        row = idx // 2
        x = col * thumb_w
        y = row * (thumb_h + label_h)
        draw.text((x + 16, y + 14), png.stem, font=label_font, fill=INK)
        slide = Image.open(png).convert("RGB").resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        sheet.paste(slide, (x, y + label_h))
    sheet.save(CONTACT_SHEET, "PNG")
    return CONTACT_SHEET


def build_slide05_shower_shape_variables(soften_panel_indices: set[int] | None = None) -> Image.Image:
    soften_panel_indices = soften_panel_indices or set()
    img = base_slide(
        "What the shower-shape variables measure",
        "A photon-like EMCal cluster has a concentrated core, narrow shoulders, and little evidence of elongation or splitting.",
    )
    add_top_right_sphenix_logo_like_slide2(img)
    draw_header_shower_icon(img)
    draw = ImageDraw.Draw(img, "RGBA")
    panel_top = 324
    panel_bottom = 1138
    gap = 38
    panel_w = (W - 2 * 132 - 2 * gap) // 3
    boxes = [
        (132, panel_top, 132 + panel_w, panel_bottom),
        (132 + panel_w + gap, panel_top, 132 + 2 * panel_w + gap, panel_bottom),
        (132 + 2 * (panel_w + gap), panel_top, W - 132, panel_bottom),
    ]
    draw_compact_core_panel(img, boxes[0])
    draw_narrow_shoulders_panel(img, boxes[1])
    draw_elongation_panel(img, boxes[2])
    for idx in sorted(soften_panel_indices):
        soften_focus_region(img, boxes[idx])
    draw.rounded_rectangle((132, 1192, W - 132, 1308), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw_wrapped(
        draw,
        "Together, these variables ask whether the EMCal energy looks like one compact photon, rather than a broad, elongated, or split decay-like shower.",
        (176, 1218),
        W - 352,
        font(TIMES_ITALIC, 34),
        fill=BLUE,
        line_gap=6,
    )
    draw_hp2026_identity_footer(img)
    return img


def slide05_focus_core_compactness() -> tuple[Path, Path]:
    img = build_slide05_shower_shape_variables({1, 2})
    png = OUTPUT / "hp2026_slide05A_shower_shape_focus_core.png"
    pdf = OUTPUT / "hp2026_slide05A_shower_shape_focus_core.pdf"
    img.convert("RGB").save(png, "PNG")
    img.convert("RGB").save(pdf, "PDF", resolution=300.0)
    script = save_script(
        5,
        "What the shower-shape variables measure - frame A",
        "For the first click, I would start with the core compactness panel. At the EMCal level, this is the first physical question: once I have an electromagnetic cluster, is the energy concentrated near the seed? et1 and E11 over E33 are both asking whether the candidate looks like one compact photon-like deposit rather than a diffuse or multi-prong shower.",
        stem="hp2026_slide05A_shower_shape_focus_core_script",
    )
    return png, script


def slide05_focus_shoulders() -> tuple[Path, Path]:
    img = build_slide05_shower_shape_variables({2})
    png = OUTPUT / "hp2026_slide05B_shower_shape_focus_shoulders.png"
    pdf = OUTPUT / "hp2026_slide05B_shower_shape_focus_shoulders.pdf"
    img.convert("RGB").save(png, "PNG")
    img.convert("RGB").save(pdf, "PDF", resolution=300.0)
    script = save_script(
        5,
        "What the shower-shape variables measure - frame B",
        "On the next click, I would add the shoulder-width panel. This is still a local EMCal question, but now the hottest tower is visually removed from the question and we look at the surrounding shower. The width variables in eta and phi tell us whether the shoulders around the seed stay narrow, which is what we expect for a single photon-like cluster.",
        stem="hp2026_slide05B_shower_shape_focus_shoulders_script",
    )
    return png, script


def slide04() -> tuple[Path, Path]:
    img = build_slide05_shower_shape_variables(set())
    png = OUTPUT / "hp2026_slide05_shower_shape_variables_intro.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        5,
        "What the shower-shape variables measure",
        "Now I want to make the photon-identification step concrete before showing any of the BDT or result plots. At the EMCal level, this analysis is asking a very local question: once I have an electromagnetic cluster, where does the energy sit inside the nearby tower pattern?\n\nI would break the shower-shape variables into three physical ideas. First, on the left, is the core compact? et1 asks how much energy is in the two-by-two core around the center of gravity, while e11 over e33 asks how dominant the center tower is inside the local three-by-three core. A prompt-photon-like shower should keep its energy concentrated near the seed.\n\nSecond, in the middle, are the shoulders narrow? The width variables in eta and phi are seed-excluded, so the hottest tower is visually removed from the question. What remains is the surrounding shower. A single photon should have narrow shoulders, while a broader or multi-tower pattern is less photon-like.\n\nThird, on the right, is the shower stretched or split? e32 over e35 compares the energy in a narrow three-by-two strip with a broader local three-by-five region. A compact photon-like cluster should keep most of its energy inside the narrow strip, while an elongated or split decay-like cluster spills more energy into the broader region.\n\nThe takeaway is that these are not arbitrary variables. Together, they ask whether the EMCal energy looks like one compact photon, rather than a broad, elongated, or split decay-like shower. With that picture in mind, the next slide can show how these handles become a quantitative photon-ID selection.",
        stem="hp2026_slide05_shower_shape_variables_intro_script",
    )
    return png, script


def soften_focus_region(base: Image.Image, region: tuple[int, int, int, int]) -> None:
    crop = base.crop(region).convert("RGBA")
    crop = crop.filter(ImageFilter.GaussianBlur(radius=1.8))
    veil = Image.new("RGBA", crop.size, (247, 250, 252, 155))
    crop.alpha_composite(veil)
    edge = ImageDraw.Draw(crop, "RGBA")
    edge.rounded_rectangle(
        (2, 2, crop.width - 3, crop.height - 3),
        radius=12,
        outline=(213, 226, 235, 185),
        width=2,
    )
    base.alpha_composite(crop, (region[0], region[1]))


def soften_plot_pad_region(base: Image.Image, placed: tuple[int, int, int, int], pad_index: int) -> None:
    x0, y0, x1, y1 = placed
    height = y1 - y0
    # The exact paper crop is one tall three-pad figure. These fractions are
    # tied to the detected horizontal pad divider lines in the fitted Fig. 8
    # image: top/middle at y ~= 459 px and middle/bottom at y ~= 703 px for
    # the 647x936 rendered crop used on the slide.
    fractions = {
        1: (0.490, 0.751),
        2: (0.751, 1.000),
    }
    if pad_index not in fractions:
        return
    fy0, fy1 = fractions[pad_index]
    region = (
        x0,
        y0 + round(height * fy0),
        x1,
        y0 + round(height * fy1),
    )
    overlay = Image.new("RGBA", base.size, (0, 0, 0, 0))
    overlay_draw = ImageDraw.Draw(overlay, "RGBA")
    overlay_draw.rectangle(region, fill=(248, 250, 252, 142))
    base.alpha_composite(overlay)


def build_slide06_bdt_conceptual(soften_card_indices: set[int] | None = None) -> Image.Image:
    soften_card_indices = soften_card_indices or set()
    img = base_slide(
        "How a BDT turns shower shapes into a photon-ID score",
        "A learned sequence of simple shower-shape questions replaces hand-tuned one-variable cuts.",
    )
    add_top_right_sphenix_logo_like_slide2(img)
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((132, 272, W - 132, 284), fill=(*SOFT_BG, 255))
    draw.line((132, 264, W - 132, 264), fill=(221, 226, 232, 255), width=3)
    draw_bdt_input_ribbon(img, (132, 284, W - 132, 444))
    card_top = 474
    card_bottom = 1132
    gap = 44
    card_w = (W - 2 * 132 - 2 * gap) // 3
    cards = [
        (132, card_top, 132 + card_w, card_bottom),
        (132 + card_w + gap, card_top, 132 + 2 * card_w + gap, card_bottom),
        (132 + 2 * (card_w + gap), card_top, W - 132, card_bottom),
    ]
    draw_manual_cuts_card(img, cards[0])
    draw_single_tree_card(img, cards[1])
    draw_boosted_score_card(img, cards[2])
    for left_box, right_box in zip(cards[:-1], cards[1:]):
        draw_arrow(draw, (left_box[2] + 8, (card_top + card_bottom) // 2), (right_box[0] - 10, (card_top + card_bottom) // 2), fill=(151, 169, 188), width=5)
    for idx in sorted(soften_card_indices):
        soften_focus_region(img, cards[idx])

    draw.rounded_rectangle((132, 1160, W - 132, 1296), radius=10, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw_wrapped(
        draw,
        "The BDT learns correlated shower-shape patterns, not just isolated thresholds.",
        (174, 1184),
        W - 348,
        font(TIMES_BOLD, 34),
        fill=BLUE,
        line_gap=4,
    )
    draw_wrapped(
        draw,
        "That lets photon ID use core compactness, shoulders, elongation, and kinematics together.",
        (174, 1238),
        W - 348,
        font(TIMES_ITALIC, 30),
        fill=MUTED,
        line_gap=4,
    )
    draw_hp2026_identity_footer(img)
    return img


def slide06_bdt_focus_manual_cuts() -> tuple[Path, Path]:
    img = build_slide06_bdt_conceptual({1, 2})
    png = OUTPUT / "hp2026_slide06A_bdt_focus_manual_cuts.png"
    pdf = OUTPUT / "hp2026_slide06A_bdt_focus_manual_cuts.pdf"
    img.convert("RGB").save(png, "PNG")
    img.convert("RGB").save(pdf, "PDF", resolution=300.0)
    script = save_script(
        6,
        "How a BDT turns shower shapes into a photon-ID score - frame A",
        "On the first click for this build, I would start with the left panel only. The familiar baseline is manual rectangular cuts: one fixed gate on a width, then another fixed gate on an energy ratio, then another fixed gate on the next shower-shape variable. That is understandable and reproducible, but it treats the variables mostly one at a time. The point of this first frame is to establish what the BDT is replacing, before asking the audience to look at the tree or the ensemble.",
        stem="hp2026_slide06A_bdt_focus_manual_cuts_script",
    )
    return png, script


def slide06_bdt_focus_tree_build() -> tuple[Path, Path]:
    img = build_slide06_bdt_conceptual({2})
    png = OUTPUT / "hp2026_slide06B_bdt_focus_tree_build.png"
    pdf = OUTPUT / "hp2026_slide06B_bdt_focus_tree_build.pdf"
    img.convert("RGB").save(png, "PNG")
    img.convert("RGB").save(pdf, "PDF", resolution=300.0)
    script = save_script(
        6,
        "How a BDT turns shower shapes into a photon-ID score - frame B",
        "On the next click, I would bring in the middle panel. This is the conceptual shift from fixed gates to a learned decision tree. The tree is still asking simple shower-shape questions, but it can ask them conditionally. Depending on the first answer, the next question can be different. Each final leaf then assigns a small signal-like or background-like score. So the method is more flexible than hand-tuned rectangular cuts, but it is still built from concrete EMCal shower information.",
        stem="hp2026_slide06B_bdt_focus_tree_build_script",
    )
    return png, script


def slide06_bdt_conceptual() -> tuple[Path, Path]:
    img = build_slide06_bdt_conceptual(set())
    png = OUTPUT / "hp2026_slide06_bdt_conceptual_photon_id_score.png"
    pdf = OUTPUT / "hp2026_slide06_bdt_conceptual_photon_id_score.pdf"
    img.convert("RGB").save(png, "PNG")
    img.convert("RGB").save(pdf, "PDF", resolution=300.0)
    script = save_script(
        6,
        "How a BDT turns shower shapes into a photon-ID score",
        "Now that the shower-shape variables have a physical meaning, the next question is how the analysis combines them into one photon-identification decision. The top band is the compact picture: measured cluster features like core compactness, shoulders, elongation, cluster transverse energy, eta, and the event vertex are passed into a learned classifier. The output is one photon-ID score that ranks the cluster from more background-like to more photon-like.\n\nThe left side shows the most familiar way to do this: manual rectangular cuts. You can cut on a width, then a core-energy ratio, then a strip-energy ratio, and so on. That is very understandable, but it is also rigid because each gate is mostly tuned one variable at a time.\n\nThe middle panel is the idea of a single decision tree. Instead of applying the same cut list to every cluster, the tree learns which question to ask next. If the width is narrow, it may then ask about the core ratio. If the width is broad, it may ask a different shower-shape question. At the bottom, each final leaf assigns a small signal-like or background-like score. The important point is that this is still built from simple, readable shower-shape questions.\n\nThe right side is the boosting step. A boosted decision tree is a committee of many shallow trees. Each tree is weak by itself, but later trees focus on clusters that earlier trees handled poorly. Their outputs are combined into one score that ranks clusters from background-like to photon-like. I would not call that score a photon probability; it is a photon-like score used to define the identification region.\n\nSo the takeaway is that the BDT learns correlated shower-shape patterns, not just isolated thresholds. It lets photon ID use the core, shoulder, elongation, and kinematic information together before we later pair the ID axis with isolation.",
        stem="hp2026_slide06_bdt_conceptual_photon_id_score_script",
    )
    return png, script


def slide07_npb_preselection_species_control() -> tuple[Path, Path]:
    img = base_slide(
        "Preselection makes sidebands comparable",
        "A mild NPB gate removes non-collision-like clusters before tight and non-tight photon-ID regions are defined.",
    )
    add_top_right_sphenix_logo_like_slide2(img)
    draw_npb_time_evidence_card(img, (132, 320, 1192, 930))
    draw_npb_threshold_evidence_card(img, (1234, 320, 2390, 930))
    draw_npb_sideband_flow_card(img, (132, 962, 2390, 1308))
    draw_hp2026_identity_footer(img)

    png = OUTPUT / "hp2026_slide07_npb_preselection_species_control.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        7,
        "Preselection makes sidebands comparable",
        "Now I want to separate two ideas that are easy to blur together. The BDT idea from the previous slide is a way of turning detector information into one score. But before the tight photon-ID BDT defines tight and non-tight photon-ID regions, the analysis first applies a non-physical-background score, or NPB score.\n\nThe point of this NPB score is not to claim photon purity. It is candidate-pool hygiene. It removes clusters that should not be part of the photon-ID sideband problem in the first place: out-of-time, non-collision-like, or otherwise pathological EMCal clusters.\n\nOn the left, the validation plot shows the key signature. High NPB-score clusters are concentrated around the in-time collision region. Low-score clusters carry the out-of-time component. So this is not just an arbitrary extra cut; the score is tied to the timing structure of the unwanted cluster population.\n\nOn the right, the threshold scan explains why the nominal requirement, NPB greater than 0.5, is deliberately mild. At that point, the physics-cluster retention stays about 99 percent, and the candidate pool is also about 99 percent clean in the NPB sense. Tightening the cut further starts to cost good physics clusters, so the goal is not to squeeze out every possible background. The goal is to remove disconnected junk before photon ID begins.\n\nThe bottom flow is the main logic I want the audience to take away. Raw EMCal candidates contain both collision-like clusters and artifacts. The NPB plus preselection gate removes the out-of-time and pathological component. After that, the photon-ID BDT can split a cleaner physics-candidate parent sample into tight and non-tight regions. That matters because the sideband method is strongest when tight and non-tight are related populations with different purity, not disconnected populations. So this slide sets up the next one: now that the candidate pool is clean, the tight photon-ID BDT can define the actual ID axis.",
        stem="hp2026_slide07_npb_preselection_species_control_script",
    )
    return png, script


def slide05() -> tuple[Path, Path]:
    img = base_slide(
        "Tight photon ID selects prompt-like showers",
        "After preselection, the photon-ID BDT defines the tight-ID axis used with isolation and sidebands.",
    )
    add_top_right_sphenix_logo(img)
    draw_public_bdt_score_card(img, (132, 320, 1132, 1228))
    draw_tight_id_before_after_card(img, (1190, 320, 2390, 880))
    draw_tight_id_training_card(img, (1190, 918, 1774, 1294))
    draw_tight_id_threshold_card(img, (1814, 918, 2390, 1294))
    draw_hp2026_identity_footer(img)

    png = OUTPUT / "hp2026_slide08_tight_photon_id_performance.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        8,
        "Tight photon ID selects prompt-like showers",
        "Now that the NPB preselection has made the candidate population photon-like, this slide is where the actual tight photon-ID selection happens. The main plot on the left is the public paper figure for the photon-ID BDT score after preselection. The red prompt-photon simulation moves toward high score, while the blue inclusive-jet background is concentrated toward low score. So the BDT is doing exactly what we want it to do: it turns the correlated EMCal shower-shape information into one reproducible identification axis.\n\nThe part I would emphasize is that this is not just a black-box score. The training target is very concrete. The signal clusters are truth-matched direct or fragmentation photons, and the background class is the complement of the prompt definition in inclusive-jet PYTHIA MC. The inputs are the shower-shape variables from the earlier slide, plus cluster transverse energy, eta, and the event vertex. So the score is learning the EMCal shower pattern of a prompt-like photon candidate.\n\nOn the upper right, I would use the E11 over E33 flow as the physical check. Before preselection, the data are visibly not described by a simple prompt-photon plus inclusive-jet mixture, and the NPB-tagged component sits in a distinct low-core region. After preselection, that pathological component is removed and the data become much more interpretable as a physics-candidate population. After tight ID, the same compact-core observable is pushed toward the prompt-like template. That is the important visual point: the selection is not just a score threshold; it changes the shower-shape distribution in the direction the physics picture says it should.\n\nThe working point is also explicit. The tight requirement is an ET-dependent BDT threshold, chosen to keep about 80 percent identification efficiency relative to preselection. That gives us a high-purity ID slice, but it does not claim the selected sample is perfectly pure. That is why the next step is isolation, and then the ID and isolation axes together define the sidebands for the data-driven purity measurement.",
        stem="hp2026_slide08_tight_photon_id_performance_script",
    )
    return png, script


def slide06() -> tuple[Path, Path]:
    title = "Isolation defines the photon sample"
    subtitle = "Quiet prompt-like candidates separate from nearby jet activity."
    img = base_slide_hp2026_main(title, subtitle)
    draw = ImageDraw.Draw(img, "RGBA")
    add_top_right_sphenix_logo_like_slide2(img)
    place_figure(img, "fig3_isolation", (132, 308, 1446, 1248), inset=18, accent=TEAL)
    draw_isolation_population_panel(img, (1504, 308, 2390, 838))
    draw_isolation_cut_logic_panel(img, (1504, 868, 2390, 1248))
    draw_hp2026_identity_footer(img)
    png = OUTPUT / "hp2026_slide09_isolation_physics_clean.png"
    img.convert("RGB").save(png, "PNG")
    write_hp2026_main_header_spec(png)
    script = save_script(
        9,
        title,
        "Now that the photon-ID BDT has defined a prompt-like cluster, the next question is whether the area around that cluster is quiet. That is what isolation measures. A prompt photon should leave the hard scattering and not carry a lot of nearby hadronic activity with it, while fragmentation photons and neutral-meson backgrounds tend to live inside a busier jet environment.\n\nThe plot on the left is the reconstructed isolation-energy distribution. The black points are the tight-ID data candidate sample after preselection. The red shaded distribution is non-tight-ID data, so it is a background-enriched sideband shape. The blue shaded distribution is the tight-ID prompt-photon signal MC. The useful visual point is that the prompt-photon template is concentrated at low isolation energy, while the background-enriched shape carries a longer high-isolation tail.\n\nThe lower-right panel makes the cut definition explicit. A candidate is called isolated when the reconstructed isolation energy is below 0.49 plus 0.037 times the photon transverse energy, and that threshold is chosen to keep about 80 percent isolation efficiency. The non-isolated sideband is separated from the cut by at least 0.8 GeV, so the sideband is not just the edge of the selected region.\n\nThis is also where the ABCD structure starts to become concrete. Tight versus non-tight photon ID gives one axis, and isolated versus non-isolated gives the second axis. Region A is the selected tight-and-isolated sample, while B, C, and D are the sideband regions that constrain the residual background. So the main message is that isolation cleans the sample, but more importantly it creates the second controlled axis needed for the data-driven purity measurement.",
    )
    return png, script


def slide07() -> tuple[Path, Path]:
    img = base_slide("Purity is measured from data, not assumed", "The residual background is constrained by the ID and isolation sidebands.")
    draw = ImageDraw.Draw(img, "RGBA")
    place_figure(img, "fig4_abcd", (132, 334, 1280, 1124))
    claim_box(draw, (1370, 334, 2390, 548), "Region A is the selected sample", "Tight-ID and isolated candidates are the signal region, but the sample still contains residual background.", accent=PHOTON)
    claim_box(draw, (1370, 600, 2390, 814), "Sidebands constrain the background", "The non-tight and non-isolated regions estimate what background remains in the signal region.", accent=SPHENIX_BLUE)
    claim_box(draw, (1370, 866, 2390, 1124), "This is the clean audience story", "The purity correction is data-driven: select a clean sample, measure contamination, then correct the yield before unfolding.", accent=TEAL)
    bridge(draw, "Next: combine purity, efficiency, and detector-response corrections into the final yield.")
    png = OUTPUT / "hp2026_slide10_data_driven_purity_abcd.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        10,
        "Purity is measured from data, not assumed",
        "This slide is the bridge from selection to measurement. The selected sample is region A: tight photon ID and isolated. But region A is not assumed to be pure. The sideband regions use the same two axes, ID and isolation, to estimate the residual background from data. This is an important point for the Hard Probes audience because it makes the result more than a simulation-based selection. The analysis measures how much of the selected candidate yield is prompt-photon signal, then applies that correction before the detector-response corrections.",
    )
    return png, script


def slide08() -> tuple[Path, Path]:
    img = base_slide("Corrections: purity, efficiency, and unfolding", "Candidate counts become a particle-level cross section through explicit correction steps.")
    draw = ImageDraw.Draw(img, "RGBA")
    place_figure(img, "fig5_purity", (132, 324, 1190, 1002))
    place_figure(img, "fig6_efficiencies", (1282, 324, 2390, 1002))
    draw_analysis_flow(img, (218, 1064, 2304, 1208))
    bridge(draw, "The result slide is now interpretable: the plotted spectrum is corrected to particle level.")
    png = OUTPUT / "hp2026_slide11_corrections_purity_efficiency_unfolding.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        11,
        "Corrections: purity, efficiency, and unfolding",
        "At this point the audience has seen how the candidate is selected and how the remaining background is estimated. This slide shows the correction accounting. The purity tells us what fraction of the selected tight-and-isolated sample is signal. The efficiency curves show the reconstruction, identification, isolation, and total efficiency corrections. Then the purity-corrected yield is unfolded for detector response. The reason to include this slide is not to derive the equation in detail, but to make the final cross section credible: every step between candidate counts and particle-level yield is visible.",
    )
    return png, script


def slide09() -> tuple[Path, Path]:
    img = base_slide("Systematics are controlled across the analysis chain", "The uncertainty budget says which analysis steps matter most in each ET region.")
    draw = ImageDraw.Draw(img, "RGBA")
    place_figure(img, "fig7_systematics", (132, 324, 1628, 1138))
    claim_box(draw, (1700, 324, 2390, 560), "Low ET", "Purity and background-closure effects dominate where the residual background is largest.", accent=PHOTON)
    claim_box(draw, (1700, 624, 2390, 860), "High ET", "Energy scale and resolution become increasingly important as the spectrum is steep and statistics shrink.", accent=SPHENIX_BLUE)
    claim_box(draw, (1700, 924, 2390, 1138), "Why this slide earns its place", "It shows the final uncertainty is a controlled budget, not a single opaque error bar.", accent=TEAL)
    bridge(draw, "Next: show the corrected isolated prompt-photon cross section.")
    png = OUTPUT / "hp2026_slide12_systematics_controlled_chain.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        12,
        "Systematics are controlled across the analysis chain",
        "Before the result, I want one slide that shows the uncertainty budget. This figure breaks down the relative systematic uncertainties across photon ET. The useful way to narrate it is by region. At lower ET, the sample has more residual background, so purity and closure-related terms matter strongly. At higher ET, the spectrum is steep and the photon energy scale and resolution become more important. This slide lets the audience see that the final error bars are not opaque. They come from identifiable analysis steps that we have already walked through.",
    )
    return png, script


def build_slide13_main_cross_section(
    soften_card_indices: set[int] | None = None,
) -> Image.Image:
    soften_card_indices = soften_card_indices or set()
    img = base_slide(
        "Main result: isolated prompt-photon cross section",
        "Read the result top-to-bottom: cross section, theory agreement, then PDF sensitivity.",
    )
    draw = ImageDraw.Draw(img, "RGBA")
    add_top_right_sphenix_logo_like_slide2(img)

    placed, panel = place_figure_snug_panel(img, "fig8_cross_section", (60, 306, 928, 1286), pad=22, crop_pad=18, trim_bottom=42)

    DATA_BLUE = (67, 132, 217)
    PYTHIA_ORANGE = (238, 126, 40)
    JETPHOX_MAGENTA = (214, 86, 169)
    VOGELSANG_GREEN = (91, 174, 96)
    PDF_RED = (242, 126, 119)
    PDF_OLIVE = (117, 143, 89)
    PDF_BLUE = (82, 88, 210)

    def color_note(
        x: int,
        y: int,
        label: str,
        text: str,
        color: tuple[int, int, int],
        max_width: int,
    ) -> int:
        draw.rounded_rectangle((x, y + 7, x + 18, y + 35), radius=6, fill=(*color, 255))
        label_font = font(TIMES_BOLD, 28)
        label_w, _ = text_box(draw, label, label_font)
        draw.text((x + 34, y), label, font=label_font, fill=color)
        return draw_wrapped(
            draw,
            text,
            (x + 54 + label_w, y + 1),
            max_width - label_w - 54,
            font(TIMES, 27),
            fill=MUTED,
            line_gap=6,
        )

    def panel_card(
        box: tuple[int, int, int, int],
        number: str,
        heading: str,
        notes: list[tuple[str, str, tuple[int, int, int]]],
        read_text: str,
        accent: tuple[int, int, int],
        target_y: int,
    ) -> None:
        x0, y0, x1, y1 = box
        draw.line((panel[2] + 12, target_y, x0 - 16, y0 + 64), fill=(*accent, 125), width=3)
        draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 255), outline=(*CARD_EDGE, 255), width=2)
        draw.rounded_rectangle((x0, y0, x0 + 14, y1), radius=7, fill=(*accent, 255))
        draw.ellipse((x0 + 42, y0 + 36, x0 + 94, y0 + 88), fill=(*accent, 255))
        nw, nh = text_box(draw, number, font(TIMES_BOLD, 32))
        draw.text((x0 + 68 - nw / 2, y0 + 45), number, font=font(TIMES_BOLD, 32), fill=(255, 255, 255))
        draw.text((x0 + 124, y0 + 36), heading, font=font(TIMES_BOLD, 38), fill=INK)
        y = y0 + 92
        for label, note, color in notes:
            y = color_note(x0 + 124, y, label, note, color, x1 - x0 - 170) + 6
        draw.line((x0 + 124, y + 6, x1 - 54, y + 6), fill=(226, 232, 238, 255), width=2)
        draw_wrapped(draw, read_text, (x0 + 124, y + 20), x1 - x0 - 180, font(TIMES_BOLD, 28), fill=BLUE, line_gap=5)

    card_x0, card_x1 = 990, 2390
    cards = [
        (card_x0, 318, card_x1, 610),
        (card_x0, 642, card_x1, 950),
        (card_x0, 982, card_x1, 1276),
    ]
    panel_card(
        cards[0],
        "1",
        "Top panel: measured spectrum",
        [
            ("Blue data:", "corrected sPHENIX isolated-photon cross section.", DATA_BLUE),
            ("Orange PYTHIA8:", "Detroit-tune Monte Carlo prediction.", PYTHIA_ORANGE),
            ("Magenta / green NLO:", "JETPHOX and Vogelsang pQCD, with the same truth-level isolation.", JETPHOX_MAGENTA),
        ],
        "64.4 inverse picobarns of p+p at √s = 200 GeV; the spectrum spans about two orders of magnitude.",
        PHOTON,
        placed[1] + 150,
    )
    panel_card(
        cards[1],
        "2",
        "Middle panel: theory / data",
        [
            ("Dashed 1.0:", "perfect theory/data agreement reference.", (60, 60, 60)),
            ("Blue band:", "experimental systematic uncertainty drawn around unity.", DATA_BLUE),
            ("Colored predictions:", "NLO pQCD remains compatible; PYTHIA8 slightly overpredicts.", PYTHIA_ORANGE),
        ],
        "Agreement is read against the 1.0 line; the measurement is consistent with NLO pQCD within quoted uncertainties.",
        SPHENIX_BLUE,
        placed[1] + 540,
    )
    panel_card(
        cards[2],
        "3",
        "Bottom panel: PDF sensitivity",
        [
            ("PDF curves:", "JETPHOX repeated with CT18NLO, NNPDF4.0, CTEQ6.6, and MSHT20NLO.", TEAL),
            ("Shaded bands:", "PDF uncertainties around each PDF-set prediction.", PDF_OLIVE),
            ("Baseline value:", "the p+p result is a reference with interpretable theory dependence.", PDF_BLUE),
        ],
        "Same observable, varied proton PDFs: the lower panel illustrates sensitivity to the proton PDF.",
        TEAL,
        placed[1] + 778,
    )
    for idx in sorted(soften_card_indices):
        soften_plot_pad_region(img, placed, idx)
    for idx in sorted(soften_card_indices):
        soften_focus_region(img, cards[idx])
    draw_hp2026_identity_footer(img)
    return img


def slide13_focus_spectrum() -> tuple[Path, Path]:
    img = build_slide13_main_cross_section({1, 2})
    png = OUTPUT / "hp2026_slide13A_main_cross_section_focus_spectrum.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        13,
        "Main result: isolated prompt-photon cross section - frame A",
        "I would start this result slide with the top panel. This is the corrected isolated prompt-photon cross section in p+p collisions at 200 GeV. The key thing to notice first is the scale of the result: over the measured transverse-energy range, the spectrum falls by about two orders of magnitude. That is why the purity, efficiency, unfolding, and energy-scale pieces we just walked through are not just technical details. They are what let this steep spectrum become a real cross-section measurement.",
        stem="hp2026_slide13A_main_cross_section_focus_spectrum_script",
    )
    return png, script


def slide13_focus_theory_data() -> tuple[Path, Path]:
    img = build_slide13_main_cross_section({2})
    png = OUTPUT / "hp2026_slide13B_main_cross_section_focus_theory_data.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        13,
        "Main result: isolated prompt-photon cross section - frame B",
        "On the next step, I would move from the cross section itself to the theory-over-data panel. The dashed line is the reference for perfect agreement, and the blue band is the experimental systematic uncertainty around unity. The important statement is that the NLO pQCD calculations are compatible with the measurement within uncertainties. PYTHIA sits somewhat high, but I would not overplay that as the main message. The main point is that the corrected sPHENIX result lands in the expected perturbative-QCD range.",
        stem="hp2026_slide13B_main_cross_section_focus_theory_data_script",
    )
    return png, script


def slide10() -> tuple[Path, Path]:
    img = build_slide13_main_cross_section(set())

    png = OUTPUT / "hp2026_slide13_main_cross_section_full_glory.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        13,
        "Main result: isolated prompt-photon cross section",
        "This is the main result of the talk. At this point, we have already seen how the candidate is defined, how the background is constrained, and how the corrections are applied, so now we can read the figure as the physics result rather than as an isolated plot.\n\nStarting from the top panel, the black points are the corrected isolated prompt-photon cross section in p+p collisions at 200 GeV. The spectrum falls quickly across the measured transverse-energy range, over roughly two orders of magnitude, which is exactly why the energy scale, purity, and unfolding pieces mattered earlier.\n\nMoving to the first ratio panel, the key statement is that the NLO pQCD calculations agree with the measurement within uncertainties. PYTHIA is slightly high, but the important public-facing message is not that one curve wins; it is that the corrected sPHENIX result is quantitatively in the expected perturbative-QCD range.\n\nFinally, the lower ratio panel shows that the result is not just a rate measurement. The theory comparison has visible dependence on the proton PDF choice, so isolated photons remain a clean hard-scattering observable with sensitivity to the ingredients of the pQCD calculation.\n\nWith the result itself visible, the final slide places it in RHIC context and closes the loop on the uncertainty budget.",
        stem="hp2026_slide13_main_cross_section_full_glory_script",
    )
    return png, script


def slide11() -> tuple[Path, Path]:
    img = base_slide(
        "RHIC comparison and uncertainty budget",
        "Fig. 9 checks consistency with PHENIX; Fig. 7 shows which systematics control the measurement.",
    )
    draw = ImageDraw.Draw(img, "RGBA")
    add_top_right_sphenix_logo_like_slide2(img)

    fig9_placed, fig9_panel = place_figure_snug_panel(
        img,
        "fig9_phenix",
        (94, 308, 1016, 1262),
        pad=20,
        crop_pad=12,
        trim_bottom=0,
    )

    def info_card(
        box: tuple[int, int, int, int],
        title: str,
        accent: tuple[int, int, int],
        body_lines: list[tuple[str, str, tuple[int, int, int]]],
        title_size: int = 32,
        label_size: int = 24,
        body_size: int = 23,
        line_step: int = 58,
    ) -> None:
        x0, y0, x1, y1 = box
        draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 255), outline=(*CARD_EDGE, 255), width=2)
        draw.rounded_rectangle((x0, y0, x0 + 12, y1), radius=6, fill=(*accent, 255))
        draw.text((x0 + 34, y0 + 24), title, font=font(TIMES_BOLD, title_size), fill=INK)
        y = y0 + 78
        for label, body, color in body_lines:
            draw.ellipse((x0 + 42, y + 7, x0 + 62, y + 27), fill=(*color, 255))
            label_font = font(TIMES_BOLD, label_size)
            label_w, _ = text_box(draw, label, label_font)
            draw.text((x0 + 78, y), label, font=label_font, fill=color)
            draw_wrapped(
                draw,
                body,
                (x0 + 92 + label_w, y + 1),
                x1 - x0 - label_w - 128,
                font(TIMES, body_size),
                fill=MUTED,
                line_gap=2,
            )
            y += line_step

    comparison_box = (1062, 308, 2390, 642)
    info_card(
        comparison_box,
        "PHENIX comparison",
        PHOTON,
        [
            ("Blue:", "sPHENIX isolated prompt photons, |η| < 0.7, 64.4 inverse picobarns at √s = 200 GeV.", (67, 132, 217)),
            ("Pink:", "published PHENIX direct photons, |η| < 0.25, with no isolation requirement.", (202, 78, 177)),
            ("Purple:", "PHENIX after bin-width and acceptance correction; the lower panel is corrected PHENIX divided by sPHENIX.", (142, 92, 220)),
        ],
        title_size=38,
        label_size=28,
        body_size=26,
        line_step=78,
    )

    fig7_placed, fig7_panel = place_figure_snug_panel(
        img,
        "fig7_systematics",
        (1062, 680, 2010, 1288),
        pad=18,
        crop_pad=10,
        trim_bottom=0,
    )
    syst_box = (2038, 680, 2390, 1288)
    draw.rounded_rectangle(syst_box, radius=10, fill=(255, 255, 255, 255), outline=(218, 226, 235, 255), width=2)
    draw.rounded_rectangle((syst_box[0], syst_box[1], syst_box[0] + 9, syst_box[3]), radius=5, fill=(55, 121, 139, 255))
    draw.text((syst_box[0] + 34, syst_box[1] + 28), "Systematics", font=font(TIMES_BOLD, 43), fill=INK)
    syst_points = [
        ("Colored lines", "one source varied"),
        ("Black envelope", "total systematic, incl. luminosity"),
        ("Y-axis", "fractional shift from nominal"),
    ]
    y = syst_box[1] + 120
    for idx, (head, body) in enumerate(syst_points):
        color = [PHOTON, INK, (55, 121, 139)][idx]
        draw.rounded_rectangle((syst_box[0] + 42, y + 8, syst_box[0] + 66, y + 32), radius=6, fill=(*color, 255))
        draw.text((syst_box[0] + 82, y), head, font=font(TIMES_BOLD, 30), fill=INK)
        draw_wrapped(
            draw,
            body,
            (syst_box[0] + 82, y + 40),
            syst_box[2] - syst_box[0] - 116,
            font(TIMES, 27),
            fill=MUTED,
            line_gap=4,
        )
        y += 156

    draw.line((fig9_panel[2] + 18, fig9_panel[1] + 280, comparison_box[0] - 18, comparison_box[1] + 118), fill=(*PHOTON, 135), width=3)
    draw_hp2026_identity_footer(img)

    png = OUTPUT / "hp2026_slide14_phenix_systematics_closing.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        14,
        "RHIC comparison and uncertainty budget",
        "Now I want to put the result in context without turning this into the final summary yet. The left plot compares the sPHENIX isolated prompt-photon cross section to the earlier PHENIX direct-photon measurement at the same collision energy. The blue points are sPHENIX in the wider |eta| less than 0.7 acceptance. The pink points are the published PHENIX result in |eta| less than 0.25 with no isolation requirement. The purple points are the PHENIX result after the bin-width and acceptance corrections described in the paper, and the lower ratio panel is corrected PHENIX divided by sPHENIX.\n\nThe compact comparison message is that the two measurements agree within uncertainties, while sPHENIX extends the rapidity acceptance and transverse-energy reach.\n\nThe lower-right plot is then the uncertainty-budget check. Each colored line is one systematic source propagated through the corrected cross section. The black envelope is the total systematic uncertainty, including the luminosity normalization. At low transverse energy, the uncertainty is driven by purity closure, which matches the region where residual background matters most. At high transverse energy, energy scale and resolution become the important terms, which is what we expect for a steeply falling spectrum.\n\nSo this slide has two jobs: it shows that the result is consistent with RHIC legacy photons, and it shows that the uncertainty budget has a physics-readable structure.",
        stem="hp2026_slide14_phenix_systematics_closing_script",
    )
    return png, script


SLIDE_RENDERERS = [
    (5, slide05_focus_core_compactness),
    (5, slide05_focus_shoulders),
    (5, slide04),
    (6, slide06_bdt_focus_manual_cuts),
    (6, slide06_bdt_focus_tree_build),
    (6, slide06_bdt_conceptual),
    (7, slide07_npb_preselection_species_control),
    (8, slide05),
    (9, slide06),
    (10, slide07),
    (11, slide08),
    (12, slide09),
    (13, slide13_focus_spectrum),
    (13, slide13_focus_theory_data),
    (13, slide10),
    (14, slide11),
]


def write_manifest(outputs: list[tuple[int, Path, Path]], contact_sheet: Path) -> Path:
    figure_entries = []
    for spec in FIGURES.values():
        figure_entries.append(
            {
                "key": spec.key,
                "label": spec.label,
                "source_type": spec.source,
                "public_status": "current-paper-source_pending_public_label_check",
                "source_pdf": str(PAPER.relative_to(ROOT)),
                "pdf_page": spec.page,
                "rendered_page": str((PAGE_ASSETS / f"paper_page_{spec.page:02d}.png").relative_to(ROOT)),
                "crop_box_px_at_3000_long_edge": list(spec.box),
                "asset": str(figure_path(spec.key).relative_to(ROOT)),
            }
        )
    ian_figure_entries = []
    for spec in IAN_FIGURES.values():
        ian_figure_entries.append(
            {
                "key": spec.key,
                "label": spec.label,
                "source_type": spec.source,
                "public_status": "placeholder_only_not_public_safe",
                "source_pdf": str(IAN.relative_to(ROOT)),
                "pdf_page": spec.page,
                "rendered_page": str((IAN_PAGE_ASSETS / f"ian_page_{spec.page:02d}.png").relative_to(ROOT)),
                "crop_box_px_at_3000_long_edge": list(spec.box),
                "asset": str(ian_figure_path(spec.key).relative_to(ROOT)),
            }
        )
    slide08_flow_entries = []
    for key, asset in slide08_e11e33_flow_paths().items():
        slide08_flow_entries.append(
            {
                "key": f"slide08_{key}_e11e33_main_panel",
                "label": f"Slide 8 E11/E33 selection-flow crop: {key.replace('_', ' ')}",
                "source_type": "existing backup-slide screenshot crop",
                "public_status": "placeholder_only_not_public_safe",
                "source_screenshot": str(SLIDE08_E11E33_FLOW_SCREENSHOT),
                "asset": str(asset.relative_to(ROOT)),
            }
        )
    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "google_slides_mutation": False,
        "deck_context": {
            "target_deck": "HPslides_v1",
            "current_existing_slides": 3,
            "generated_remaining_slides": [n for n, _, _ in outputs],
            "talk_timing": "15 minutes speaking plus 5 minutes questions",
        },
        "source_policy": {
            "main_plot_posture": "paper-first",
            "main_physics_plot_source": str(PAPER.relative_to(ROOT)),
            "private_method_detail_reference": str(IAN.relative_to(ROOT)),
            "public_hp_plot_rule": "Final public HP PPG12 plot images must come from the current paper PDF or another explicitly public/approved source recorded in this manifest.",
            "ian_and_screenshot_plots_public_status": "placeholder_only_not_public_safe",
            "requires_preliminary_or_public_labels": True,
            "no_internal_root_or_data_plots": True,
            "generated_graphics_role": "explanatory scaffolding only; not physics result plots",
        },
        "outputs": [
            {
                "slide": n,
                "png": str(png.relative_to(ROOT)),
                **({"pdf": str(png.with_suffix(".pdf").relative_to(ROOT))} if png.with_suffix(".pdf").exists() else {}),
                "speaker_script": str(script.relative_to(ROOT)),
                "size": [W, H],
                "mode": "RGB",
            }
            for n, png, script in outputs
        ],
        "qa_contact_sheet": str(contact_sheet.relative_to(ROOT)),
        "embedded_assets": figure_entries
        + ian_figure_entries
        + slide08_flow_entries
        + [
            {
                "source_type": "generated explanatory graphic",
                "description": "PIL-drawn shower-shape panels, conceptual BDT gates/tree/ensemble, NPB species-control diagram, BDT flow, isolation cone, correction flow, and closing baseline diagram.",
            }
        ],
        "notes": [
            "No slide numbers or provenance footers are baked into the PNGs.",
            "All main physics/result plots are cropped from the current PPG12 paper draft, not regenerated from internal ROOT/data outputs.",
            "IAN-derived figures and screenshot crops are placeholders only and must be removed, replaced by current-paper/public-approved figures, or explicitly approved before public HP use.",
            "Slides containing paper plots with sPHENIX Internal labels are incomplete for public HP use until the plots carry Preliminary/approved public labels or documented collaboration approval.",
        ],
    }
    path = OUTPUT / "manifest.json"
    with path.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return path


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    ASSETS.mkdir(parents=True, exist_ok=True)
    prepare_figures()
    prepare_ian_figures()
    outputs = []
    for slide_no, renderer in SLIDE_RENDERERS:
        png, script = renderer()
        outputs.append((slide_no, png, script))
    contact_sheet = write_contact_sheet(outputs)
    manifest = write_manifest(outputs, contact_sheet)
    print(manifest)
    for _, png, _ in outputs:
        print(png)
    print(contact_sheet)


if __name__ == "__main__":
    main()

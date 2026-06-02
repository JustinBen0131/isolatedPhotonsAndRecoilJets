#!/usr/bin/env python3
"""Render HP2026 Slides 4-11 full-slide PNG candidates.

Main physics/result plots are cropped from the current PPG12 paper draft.
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

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "AGENTS.md").exists())
WORKSPACE = ROOT / "outputs/manual-20260601-hp2026-fulltalk/presentations/hp2026-fulltalk"
OUTPUT = WORKSPACE / "output"
ASSETS = WORKSPACE / "assets"
PAGE_ASSETS = ASSETS / "paper_pages"
FIGURE_ASSETS = ASSETS / "paper_figures"
SCRIPT_DIR = OUTPUT / "speaker_scripts"
CONTACT_SHEET = OUTPUT / "hp2026_slides04_11_contact_sheet.png"

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


@dataclass(frozen=True)
class FigureSpec:
    key: str
    label: str
    page: int
    box: tuple[int, int, int, int]
    source: str = "paper plot"


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


def font(path: Path, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(path), size)


def text_box(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.ImageFont) -> tuple[int, int]:
    box = draw.textbbox((0, 0), text, font=fnt)
    return box[2] - box[0], box[3] - box[1]


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


def base_slide(title: str, subtitle: str | None = None) -> Image.Image:
    img = Image.new("RGBA", (W, H), (*SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    header(draw, title, subtitle)
    return img


def rounded_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], fill=(255, 255, 255), radius=10) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=(*fill, 255), outline=(*PANEL_EDGE, 255), width=2)


def shadow(base: Image.Image, box: tuple[int, int, int, int], radius: int = 10) -> None:
    layer = Image.new("RGBA", base.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(layer, "RGBA")
    d.rounded_rectangle((box[0] + 8, box[1] + 10, box[2] + 8, box[3] + 10), radius=radius, fill=(30, 42, 58, 30))
    layer = layer.filter(ImageFilter.GaussianBlur(12))
    base.alpha_composite(layer)


def figure_path(key: str) -> Path:
    return FIGURE_ASSETS / f"{key}.png"


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


def prepare_figures() -> dict[str, Path]:
    return {key: crop_figure(spec) for key, spec in FIGURES.items()}


def place_figure(
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
    img = Image.open(figure_path(key)).convert("RGBA")
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


def save_script(slide_no: int, title: str, body: str) -> Path:
    SCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    path = SCRIPT_DIR / f"hp2026_slide{slide_no:02d}_script.md"
    path.write_text(f"# HP2026 Slide {slide_no} Speaker Script\n\n{body.strip()}\n", encoding="utf-8")
    return path


def write_contact_sheet(outputs: list[tuple[int, Path, Path]]) -> Path:
    thumb_w, thumb_h = 640, 360
    label_h = 52
    sheet = Image.new("RGB", (2 * thumb_w, 4 * (thumb_h + label_h)), "white")
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


def slide04() -> tuple[Path, Path]:
    img = base_slide("From EMCal clusters to prompt-photon candidates", "The first handle is shower shape: prompt photons and decay backgrounds look different in the EMCal.")
    draw = ImageDraw.Draw(img, "RGBA")
    place_figure(img, "fig1_shower_shape", (132, 332, 1542, 1118))
    claim_box(draw, (1600, 336, 2390, 548), "A measurable shape handle", "Prompt photons deposit a narrow, single-core shower. Neutral-meson decays tend to leave wider or merged electromagnetic clusters.", accent=PHOTON)
    draw_shower_pair(img, (1600, 594, 2390, 890))
    claim_box(draw, (1600, 938, 2390, 1118), "Talk logic", "Before showing the cross section, we show how the analysis turns an EMCal cluster into a photon candidate.", accent=TEAL)
    bridge(draw, "Next: compress the shower-shape information into the photon-identification BDT.")
    png = OUTPUT / "hp2026_slide04_emcal_clusters_to_prompt_candidates.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        4,
        "From EMCal clusters to prompt-photon candidates",
        "The analysis starts from EMCal clusters. The important experimental fact is that prompt photons and decay backgrounds do not populate shower-shape space in the same way. A prompt photon tends to look like a narrow single electromagnetic shower, while neutral-meson decay photons can merge into a wider or multi-core cluster. This gives us a measurable handle before we even talk about the final cross section. The paper figure shows representative shower-shape variables after the non-collision background and pre-selection cuts. The key transition is that these variables become the input to a quantitative photon-identification discriminator.",
    )
    return png, script


def slide05() -> tuple[Path, Path]:
    img = base_slide("BDT photon ID turns shower shape into a quantitative cut", "The BDT is not the result; it is the compact accounting tool that makes the photon selection reproducible.")
    draw = ImageDraw.Draw(img, "RGBA")
    place_figure(img, "fig2_bdt_score", (132, 324, 1390, 1138))
    claim_box(draw, (1460, 324, 2390, 540), "High score is photon-like", "The BDT score separates prompt-photon signal MC from inclusive-jet background while keeping the data comparison visible.", accent=PHOTON)
    draw_bdt_flow(img, (1460, 590, 2390, 838))
    claim_box(draw, (1460, 890, 2390, 1138), "Why this is stronger than a cut list", "It lets the talk say: many shower-shape observables are combined into one calibrated ID axis, then paired with isolation for a data-driven purity estimate.", accent=SPHENIX_BLUE)
    bridge(draw, "Next: require the photon candidate to be isolated, so the selected object is physics-clean.")
    png = OUTPUT / "hp2026_slide05_bdt_photon_id_quantitative_cut.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        5,
        "BDT photon ID turns shower shape into a quantitative cut",
        "Instead of presenting the BDT as a black-box machine-learning slide, I want to frame it as accounting. The EMCal shower variables from the previous slide are compressed into a single score. High BDT score means the cluster is more photon-like, while lower scores are more consistent with the inclusive-jet background sample. This is the part of the analysis where we can be more quantitative than a simple visual shower-shape argument. The important message is that the BDT defines the photon-ID axis used later with isolation, rather than being the final physics result by itself.",
    )
    return png, script


def slide06() -> tuple[Path, Path]:
    img = base_slide("Isolation makes the photon sample physics-clean", "A quiet cone suppresses fragmentation and decay-rich activity around the candidate.")
    draw = ImageDraw.Draw(img, "RGBA")
    place_figure(img, "fig3_isolation", (132, 324, 1424, 1134))
    claim_box(draw, (1490, 324, 2390, 546), "The second analysis axis", "Photon ID asks whether the cluster is photon-like. Isolation asks whether the nearby event activity is quiet enough.", accent=SPHENIX_BLUE)
    draw_isolation_cone(img, (1490, 594, 2390, 902))
    claim_box(draw, (1490, 950, 2390, 1134), "Why it matters", "The isolated tight sample is cleaner, but still needs a measured residual-background correction.", accent=PHOTON)
    bridge(draw, "Next: use ID and isolation together to measure purity from sidebands.")
    png = OUTPUT / "hp2026_slide06_isolation_physics_clean.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        6,
        "Isolation makes the photon sample physics-clean",
        "The BDT gives us the photon-ID axis. Isolation gives us the second axis. The isolation requirement asks whether the calorimeter energy around the candidate is quiet. That suppresses fragmentation photons and high-pT neutral meson backgrounds, which tend to come with nearby activity. The paper figure shows the reconstructed isolation energy distributions for the tight-ID data, non-tight-ID background-enriched data, and signal MC. The point for the audience is simple: ID and isolation together make the selected sample much cleaner, but not perfectly pure. That is why the next step is a data-driven purity estimate.",
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
    png = OUTPUT / "hp2026_slide07_data_driven_purity_abcd.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        7,
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
    png = OUTPUT / "hp2026_slide08_corrections_purity_efficiency_unfolding.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        8,
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
    png = OUTPUT / "hp2026_slide09_systematics_controlled_chain.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        9,
        "Systematics are controlled across the analysis chain",
        "Before the result, I want one slide that shows the uncertainty budget. This figure breaks down the relative systematic uncertainties across photon ET. The useful way to narrate it is by region. At lower ET, the sample has more residual background, so purity and closure-related terms matter strongly. At higher ET, the spectrum is steep and the photon energy scale and resolution become more important. This slide lets the audience see that the final error bars are not opaque. They come from identifiable analysis steps that we have already walked through.",
    )
    return png, script


def slide10() -> tuple[Path, Path]:
    img = base_slide("Main result: isolated prompt-photon cross section", "The corrected p+p spectrum agrees with NLO pQCD within uncertainties.")
    draw = ImageDraw.Draw(img, "RGBA")
    place_figure(img, "fig8_cross_section", (132, 314, 1468, 1258))
    claim_box(draw, (1540, 324, 2390, 544), "What is measured", "differential cross section for isolated prompt photons in p+p at 200 GeV, with |eta| < 0.7 and 12 < ET < 32 GeV.", accent=PHOTON)
    claim_box(draw, (1540, 610, 2390, 830), "What it means", "The result is consistent with NLO pQCD calculations within the quoted experimental and theory uncertainties.", accent=SPHENIX_BLUE)
    claim_box(draw, (1540, 896, 2390, 1158), "Why it matters", "This is the p+p hard-scattering baseline needed before future sPHENIX gamma-jet and heavy-ion photon measurements.", accent=TEAL)
    png = OUTPUT / "hp2026_slide10_main_cross_section_result.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        10,
        "Main result: isolated prompt-photon cross section",
        "This is the central result slide. The plot shows the isolated prompt-photon cross section in p+p collisions at 200 GeV. The audience has now seen enough of the analysis chain to interpret the points: photon ID, isolation, data-driven purity, efficiencies, unfolding, and systematics. The main physics statement is that the corrected sPHENIX spectrum agrees with NLO pQCD predictions within uncertainties. This result is valuable by itself, and it is also the baseline measurement that future gamma-jet and heavy-ion photon analyses will build on.",
    )
    return png, script


def slide11() -> tuple[Path, Path]:
    img = base_slide("What this establishes for sPHENIX hard probes", "The photon cross section becomes the p+p baseline for the next measurements.")
    draw = ImageDraw.Draw(img, "RGBA")
    place_figure(img, "fig9_phenix", (1240, 326, 2390, 1188))
    claim_box(draw, (132, 326, 1138, 532), "First p+p photon baseline", "sPHENIX now has an isolated prompt-photon cross-section measurement at RHIC energy.", accent=PHOTON)
    claim_box(draw, (132, 590, 1138, 796), "Analysis chain is visible", "BDT photon ID, isolation, data-driven purity, efficiency corrections, and unfolding form a coherent measurement path.", accent=SPHENIX_BLUE)
    claim_box(draw, (132, 854, 1138, 1058), "Future hard-probes payoff", "The same baseline supports future gamma-jet and heavy-ion photon measurements without making this talk a photon+jet talk.", accent=TEAL)
    draw_closing_visual(img, (132, 1100, 1138, 1288))
    png = OUTPUT / "hp2026_slide11_sphenix_hard_probes_baseline.png"
    img.convert("RGB").save(png, "PNG")
    script = save_script(
        11,
        "What this establishes for sPHENIX hard probes",
        "I would close by returning to the role of the measurement. This is not just a method exercise; it establishes a p+p photon baseline for sPHENIX at RHIC. The comparison with PHENIX places the result in the RHIC context, while the analysis chain we walked through explains why the sPHENIX measurement is interpretable: photon ID, isolation, data-driven purity, corrections, and unfolding. The final forward-looking message is that this baseline is what future gamma-jet and heavy-ion photon measurements need, but the talk itself stays focused on the isolated prompt-photon cross section.",
    )
    return png, script


SLIDE_RENDERERS = [
    (4, slide04),
    (5, slide05),
    (6, slide06),
    (7, slide07),
    (8, slide08),
    (9, slide09),
    (10, slide10),
    (11, slide11),
]


def write_manifest(outputs: list[tuple[int, Path, Path]], contact_sheet: Path) -> Path:
    figure_entries = []
    for spec in FIGURES.values():
        figure_entries.append(
            {
                "key": spec.key,
                "label": spec.label,
                "source_type": spec.source,
                "source_pdf": str(PAPER.relative_to(ROOT)),
                "pdf_page": spec.page,
                "rendered_page": str((PAGE_ASSETS / f"paper_page_{spec.page:02d}.png").relative_to(ROOT)),
                "crop_box_px_at_3000_long_edge": list(spec.box),
                "asset": str(figure_path(spec.key).relative_to(ROOT)),
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
            "method_detail_reference": str(IAN.relative_to(ROOT)),
            "no_internal_root_or_data_plots": True,
            "generated_graphics_role": "explanatory scaffolding only; not physics result plots",
        },
        "outputs": [
            {
                "slide": n,
                "png": str(png.relative_to(ROOT)),
                "speaker_script": str(script.relative_to(ROOT)),
                "size": [W, H],
                "mode": "RGB",
            }
            for n, png, script in outputs
        ],
        "qa_contact_sheet": str(contact_sheet.relative_to(ROOT)),
        "embedded_assets": figure_entries
        + [
            {
                "source_type": "generated explanatory graphic",
                "description": "PIL-drawn shower sketches, BDT flow, isolation cone, correction flow, and closing baseline diagram.",
            }
        ],
        "notes": [
            "No slide numbers or provenance footers are baked into the PNGs.",
            "All main physics/result plots are cropped from the current PPG12 paper draft, not regenerated from internal ROOT/data outputs.",
            "IAN is used only as method/detail reference for wording and slide-story choices in this batch.",
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

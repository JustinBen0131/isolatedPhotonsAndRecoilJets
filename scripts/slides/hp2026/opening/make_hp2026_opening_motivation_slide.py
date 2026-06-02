#!/usr/bin/env python3
"""Render local HP2026 opening-sequence PNG candidates.

The slides are deterministic: text and schematic geometry are drawn with
Pillow, with no AI-generated logos or copied physics figures. They are intended
as local candidates only; Google Slides insertion is a separate approval step.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw, ImageEnhance, ImageFilter, ImageFont


W, H = 2560, 1440

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "AGENTS.md").exists())
DEFAULT_WORKSPACE = ROOT / "outputs/manual-20260601-hp2026-opening-slide/presentations/hp2026-opening-slide"
DEFAULT_OUTPUT = DEFAULT_WORKSPACE / "output"
ASSET_DIR = DEFAULT_WORKSPACE / "assets"
TITLE_ASSET_DIR = ROOT / "outputs/manual-20260601-hp2026-title/presentations/hp2026-title-slide/assets"

FONT_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES = FONT_DIR / "Times New Roman.ttf"
TIMES_BOLD = FONT_DIR / "Times New Roman Bold.ttf"
TIMES_ITALIC = FONT_DIR / "Times New Roman Italic.ttf"
ARIAL = FONT_DIR / "Arial.ttf"
ARIAL_BOLD = FONT_DIR / "Arial Bold.ttf"

INK = (18, 22, 28)
MUTED = (71, 79, 91)
LIGHT_MUTED = (117, 126, 138)
BLUE = (19, 41, 75)
SPHENIX_BLUE = (30, 143, 214)
TEAL = (18, 112, 142)
TEAL_SOFT = (102, 173, 190)
PHOTON = (245, 181, 34)
PHOTON_DARK = (196, 124, 10)
SOFT_BG = (252, 253, 254)
PANEL = (246, 249, 252)
PANEL_EDGE = (218, 226, 235)
CARD = (255, 255, 255)
CARD_EDGE = (223, 229, 236)

TITLE = "Why isolated prompt photons?"
SUBTITLE = "A color-neutral hard probe for the p+p baseline at RHIC"
BRIDGE = "Next: how sPHENIX turns EMCal clusters into isolated prompt-photon candidates."

REASONS = [
    (
        "Calibrated hard scale",
        "The photon tags the hard scattering with minimal final-state interaction.",
    ),
    (
        "Isolation sharpens the sample",
        "A quiet cone suppresses decay and fragmentation photons before purity correction.",
    ),
    (
        "p+p baseline for future hard probes",
        "The cross section anchors RHIC comparisons and later gamma-jet measurements.",
    ),
]

CONTEXT_TITLE = "A p+p photon baseline for the sPHENIX hard-probes program"
CONTEXT_SUBTITLE = "RHIC p+p collisions at 200 GeV are the reference system for calibrated hard scattering."
CONTEXT_BRIDGE = "Next: define the prompt-photon object and the isolation requirement."

CONTEXT_REASONS = [
    (
        "Hard probes need a reference",
        "p+p establishes the baseline before interpreting nuclear modifications in heavy-ion data.",
    ),
    (
        "Photons preserve the hard scale",
        "Prompt photons escape color final-state interactions; isolation suppresses decay and fragmentation backgrounds.",
    ),
    (
        "sPHENIX makes it measurable",
        "EMCal shower shape, isolation, and purity plus correction steps turn candidates into a cross section.",
    ),
]

REAL_DETECTOR_PHOTO = ASSET_DIR / "bnl_sphenix_banner_detector.jpg"
REAL_DETECTOR_RENDERING = ASSET_DIR / "bnl_sphenix_detector_rendering.jpg"
REAL_DETECTOR_SOURCE = {
    "primary_image": "https://www.bnl.gov/rhic/images/banner-sphenix-2.jpg",
    "primary_page": "https://www.bnl.gov/rhic/sphenix.php",
    "backup_photo": "https://www.bnl.gov/today/body_pics/2023/04/sphenix-assembly-hr.jpg",
    "backup_page": "https://www.bnl.gov/newsroom/news.php?a=221191",
    "rendering": "https://www.bnl.gov/rhic/images/spheix-rendering.jpg",
}

REAL_DETECTOR_TITLE = "A p+p photon baseline for sPHENIX hard probes"
REAL_DETECTOR_SUBTITLE = "Start from the detector and the hard-probes program, then narrow to isolated prompt photons."
REAL_DETECTOR_BRIDGE = "Talk path: detector context -> photon isolation -> purity/corrections -> p+p cross section."
REAL_DETECTOR_REASONS = [
    (
        "Hard-probes program",
        "sPHENIX was built for high-rate RHIC measurements of jets, photons, and heavy flavor.",
    ),
    (
        "Photon channel",
        "An isolated prompt photon tags the hard scattering without color final-state interactions.",
    ),
    (
        "This measurement",
        "p+p at 200 GeV gives the corrected cross-section baseline for later gamma-jet and heavy-ion comparisons.",
    ),
]

DETECTOR_DATA_TITLE = "sPHENIX at RHIC: detector and p+p data for photons"
DETECTOR_DATA_SUBTITLE = "Experiment context first: the detector handles the photon, the p+p sample anchors the baseline."
DETECTOR_DATA_BRIDGE = "With the detector and p+p reference sample established, the talk narrows to isolated prompt photons."
DETECTOR_DATA_CARDS = [
    (
        "New RHIC detector",
        "Full calorimetry and modern tracking give sPHENIX the ingredients for hard-probe measurements.",
    ),
    (
        "Photon-relevant handles",
        "EMCal shower shape plus EMCal/HCal isolation separate photon candidates from nearby activity.",
    ),
    (
        "p+p reference sample",
        "The 2024 and 2026 p+p running provide the data context for the prompt-photon cross section.",
    ),
]
DETECTOR_CALLOUTS = [
    ("EMCal", "energy + shower shape", PHOTON),
    ("HCal", "isolation energy", TEAL),
    ("Tracking", "charged-particle context", SPHENIX_BLUE),
    ("MBD / trigger", "trigger / luminosity", (134, 111, 186)),
]
PP_DATA_ROWS = [
    ("2024 p+p", "107/pb calorimeter-only", "13/pb all subsystems", SPHENIX_BLUE),
    ("2026 p+p", "17/pb all subsystems", "current baseline extension", PHOTON_DARK),
]
MUTED_PROGRAM_CONTEXT = "Broader sPHENIX running: 2023 commissioning; 2025/26 heavy-ion program kept as context, not the focus here."
DATA_SOURCE_EVIDENCE = {
    "values_on_slide": [
        "2024 p+p: calorimeter-only 107/pb; all subsystem 13/pb",
        "2026 p+p: all subsystem 17/pb",
    ],
    "source_files": [
        "usefulDocs/20260506_DIS_YeonjuGo.pdf",
        "usefulDocs/sPHENIX_AUM2026_jet_and_photon_measurement.pdf",
    ],
    "verification_note": (
        "Text extraction from Yeonju DIS2026 and Hanpu AUM2026 reference decks "
        "shows the same p+p values; Au+Au/O+O numbers were deliberately omitted "
        "from the main slide to keep Slide 2 focused on the p+p prompt-photon baseline."
    ),
}

PROGRESSIVE_TITLE = "sPHENIX at RHIC: detector and data for photons"
PROGRESSIVE_SUBTITLE_A = "First, the detector handles: tracking, calorimetry, trigger and luminosity context."
PROGRESSIVE_SUBTITLE_B = "Then, the data context: p+p is the measurement system; heavy ions motivate the next step."
PROGRESSIVE_BRIDGE_A = "Click: add the data-taking context without changing the detector anchor."
PROGRESSIVE_BRIDGE_B = "Next: with experiment and data context established, define the isolated prompt photon."

PROGRESSIVE_SYSTEM_CALLOUTS = [
    ("Tracking", "MVTX / INTT / TPC", SPHENIX_BLUE, (0.46, 0.47), (210, 390)),
    ("EMCal", "energy + shower shape", PHOTON, (0.51, 0.36), (210, 520)),
    ("HCal", "iHCal / oHCal isolation", TEAL, (0.70, 0.44), (210, 650)),
    ("Magnet", "1.4 T solenoid", BLUE, (0.72, 0.57), (950, 390)),
    ("MBD / trigger", "luminosity + event timing", (134, 111, 186), (0.67, 0.12), (950, 650)),
]

PROGRESSIVE_DATA_EVIDENCE = {
    "values_on_slide": [
        "2024 p+p: 107 pb^-1 calorimeter-only; 13 pb^-1 all subsystems",
        "2026 p+p: 17 pb^-1 all subsystems",
        "PPG12 Run 24 analysis sample: L = 64.4 pb^-1",
        "2025 Au+Au: 6.6 nb^-1 all subsystems",
        "2026 O+O: 52.6 nb^-1 calorimeter-only; 23.6 nb^-1 all subsystems",
    ],
    "source_files": [
        "usefulDocs/20260506_DIS_YeonjuGo.pdf",
        "usefulDocs/sPHENIX_AUM2026_jet_and_photon_measurement.pdf",
        "usefulDocs/sPHENIX_PPG12_Paper_2026-05-21_current_draft.pdf",
        "usefulDocs/PPG12_analysis_note_2026-05-21_v4_current_IAN.pdf",
    ],
    "verification_note": (
        "Yeonju DIS2026 and Hanpu AUM2026 reference decks provide the broad "
        "sPHENIX data-taking values. The current PPG12 paper draft and IAN "
        "separately quote the analyzed Run 24 prompt-photon luminosity as "
        "about 64.4 pb^-1."
    ),
}

STANDALONE_EXPERIMENT_TITLE = "sPHENIX at RHIC: experiment for photons"
STANDALONE_EXPERIMENT_SUBTITLE = (
    "A new RHIC detector with tracking, calorimetry, trigger and luminosity systems; these are the components used later in the photon measurement."
)
STANDALONE_EXPERIMENT_BRIDGE = "The data slide separates broad p+p availability from the Run 24 analysis sample."
STANDALONE_EXPERIMENT_TAKEAWAYS = [
    ("EMCal measures the photon candidate", "energy and shower shape are the first experimental handles", PHOTON),
    ("HCal helps define isolation", "nearby calorimeter activity tests whether the photon is experimentally quiet", TEAL),
    ("Tracking and MBD anchor the event", "vertex, trigger, timing and luminosity make a cross-section measurement possible", SPHENIX_BLUE),
]
EXPERIMENT_INTRO_CARDS = [
    (
        "Tracking",
        "MVTX / INTT / TPC / TPOT",
        "Vertex detectors, silicon tracking, the Time Projection Chamber, and TPOT reconstruct charged tracks and the collision vertex.",
        SPHENIX_BLUE,
        (0.47, 0.50),
        (132, 354),
        "left",
    ),
    (
        "EMCal",
        "Electromagnetic Calorimeter",
        "Measures photon and electron showers; provides the cluster energy and shower-shape information used for photon identification.",
        PHOTON,
        (0.42, 0.36),
        (132, 608),
        "left",
    ),
    (
        "HCal",
        "iHCal / oHCal",
        "Inner and outer Hadronic Calorimeters sample hadronic activity around the photon candidate, which helps define isolation.",
        TEAL,
        (0.73, 0.55),
        (132, 862),
        "left",
    ),
    (
        "Magnet",
        "1.4 T solenoid",
        "The solenoidal field bends charged particles, enabling momentum reconstruction in the tracking system.",
        BLUE,
        (0.64, 0.53),
        (2024, 354),
        "right",
    ),
    (
        "MBD",
        "Minimum Bias Detector",
        "Forward detector near the beam axis used for triggering, event timing, vertex context, and luminosity normalization.",
        (134, 111, 186),
        (0.60, 0.33),
        (2024, 608),
        "right",
    ),
    (
        "Barrel coverage",
        "Large-acceptance hard-probe detector",
        "The central tracking-plus-calorimetry barrel gives the large-acceptance context for the p+p photon measurement.",
        (88, 128, 118),
        None,
        (2024, 862),
        "right",
    ),
]
EXPERIMENT_INTRO_CARD_HEIGHTS = {
    "EMCal": 204,
    "Barrel coverage": 204,
}

SUBSYSTEM_GROUPS = [
    (
        "tracking",
        "Tracking system",
        "MVTX, INTT, TPC, TPOT",
        "Silicon vertex detectors and time-projection chamber inside a 1.4 T solenoid.",
    ),
    (
        "calorimetry",
        "Calorimetry",
        "EMCal, HCal",
        "Electromagnetic and hadronic calorimeters, including inner and outer HCal layers.",
    ),
    (
        "forward",
        "Forward detectors",
        "MBD, sEPD, ZDC",
        "Minimum-bias triggers, centrality context, event-plane information and luminosity.",
    ),
]

SUBSYSTEM_LABEL_ARROWS = [
    ("outer HCal", (1588, 392), (0.46, 0.26), "center"),
    ("inner HCal", (1348, 586), (0.48, 0.39), "left"),
    ("MVTX & INTT", (1236, 1124), (0.47, 0.54), "center"),
    ("TPC", (1850, 1128), (0.54, 0.54), "center"),
    ("EMCal", (2240, 1000), (0.56, 0.42), "center"),
    ("solenoid", (2248, 726), (0.70, 0.50), "center"),
]

STANDALONE_DATASET_TITLE = "Dataset context: p+p anchors this measurement"
STANDALONE_DATASET_SUBTITLE = (
    "Keep the total available p+p running separate from the smaller, defined PPG12 analysis sample."
)
STANDALONE_DATASET_BRIDGE = "With the p+p sample established, the talk narrows to the isolated prompt-photon object."
STANDALONE_PP_AVAILABILITY = [
    ("2024 p+p", "107", "pb^-1", "calorimeter-only", SPHENIX_BLUE),
    ("2024 p+p", "13", "pb^-1", "all subsystems", PHOTON),
    ("2026 p+p", "17", "pb^-1", "all subsystems", PHOTON_DARK),
]
STANDALONE_HEAVY_ION_CONTEXT = [
    ("2025 Au+Au", "6.6 nb^-1", "all subsystems"),
    ("2026 O+O", "52.6 nb^-1", "calorimeter-only"),
    ("2026 O+O", "23.6 nb^-1", "all subsystems"),
]
STANDALONE_DATASET_EVIDENCE = {
    "values_on_slide": [
        "2024 p+p: 107 pb^-1 calorimeter-only; 13 pb^-1 all subsystems",
        "2026 p+p: 17 pb^-1 all subsystems",
        "PPG12 Run 24 analysis sample: L = 64.4 pb^-1",
        "2025 Au+Au: 6.6 nb^-1 all subsystems",
        "2026 O+O: 52.6 nb^-1 calorimeter-only; 23.6 nb^-1 all subsystems",
    ],
    "source_files": [
        "usefulDocs/20260506_DIS_YeonjuGo.pdf",
        "usefulDocs/sPHENIX_AUM2026_jet_and_photon_measurement.pdf",
        "usefulDocs/sPHENIX_PPG12_Paper_2026-05-21_current_draft.pdf",
        "usefulDocs/PPG12_analysis_note_2026-05-21_v4_current_IAN.pdf",
    ],
    "verification_note": (
        "Yeonju DIS2026 and Hanpu AUM2026 agree on the broad sPHENIX p+p, Au+Au, and O+O data-taking values. "
        "The current PPG12 paper draft and IAN quote the Run 24 isolated prompt-photon analysis luminosity as about 64.4 pb^-1."
    ),
}

PHOTON_MOTIVATION_TITLE = "Why isolated prompt photons?"
PHOTON_MOTIVATION_SUBTITLE = "A calibrated hard-scattering probe for the p+p baseline at RHIC."
PHOTON_MOTIVATION_BRIDGE = "Next: reconstruct the photon candidate, estimate purity, correct to particle level, and compare the cross section."
PHOTON_MOTIVATION_REASONS = [
    (
        "Prompt photon",
        "Produced in the hard scattering and measured without strong final-state interaction.",
    ),
    (
        "Isolation",
        "A quiet cone suppresses decay photons and fragmentation-rich backgrounds before purity correction.",
    ),
    (
        "p+p baseline",
        "The corrected cross section anchors future gamma-jet and heavy-ion comparisons.",
    ),
]


def font(path: Path, size: int) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(path), size)


def text_box(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.ImageFont) -> tuple[int, int]:
    box = draw.textbbox((0, 0), text, font=fnt)
    return box[2] - box[0], box[3] - box[1]


def _font_size(fnt: ImageFont.ImageFont) -> int:
    return int(getattr(fnt, "size", 24))


def _script_font(fnt: ImageFont.ImageFont, scale: float = 0.58) -> ImageFont.FreeTypeFont:
    size = max(10, round(_font_size(fnt) * scale))
    return font(TIMES_BOLD if getattr(fnt, "path", "") == str(TIMES_BOLD) else TIMES, size)


def _next_script_run(text: str, idx: int) -> tuple[str, int]:
    idx += 1
    if idx >= len(text):
        return "", idx
    if text[idx] == "{":
        end = text.find("}", idx + 1)
        if end != -1:
            return text[idx + 1 : end], end + 1
    start = idx
    if text[idx] in "+-":
        idx += 1
    while idx < len(text) and text[idx].isalnum():
        idx += 1
    if idx == start:
        idx += 1
    return text[start:idx], idx


def rich_text_box(
    draw: ImageDraw.ImageDraw,
    text: str,
    fnt: ImageFont.ImageFont,
    sup_font: ImageFont.ImageFont | None = None,
    sub_font: ImageFont.ImageFont | None = None,
) -> tuple[int, int]:
    sup_font = sup_font or _script_font(fnt)
    sub_font = sub_font or _script_font(fnt)
    x = 0
    max_above = 0
    max_below = text_box(draw, "Ag", fnt)[1]
    i = 0
    while i < len(text):
        marker = text[i]
        if marker in "^_":
            run, i = _next_script_run(text, i)
            run_font = sup_font if marker == "^" else sub_font
            rw, rh = text_box(draw, run, run_font)
            x += rw
            if marker == "^":
                max_above = max(max_above, round(_font_size(fnt) * 0.20))
            else:
                max_below = max(max_below, text_box(draw, "Ag", fnt)[1] + round(_font_size(fnt) * 0.34) + rh // 3)
            continue
        start = i
        while i < len(text) and text[i] not in "^_":
            i += 1
        rw, _ = text_box(draw, text[start:i], fnt)
        x += rw
    return x, max_above + max_below


def draw_rich_text(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    fnt: ImageFont.ImageFont,
    fill: tuple[int, int, int],
    sup_font: ImageFont.ImageFont | None = None,
    sub_font: ImageFont.ImageFont | None = None,
) -> tuple[int, int]:
    sup_font = sup_font or _script_font(fnt)
    sub_font = sub_font or _script_font(fnt)
    x, y = xy
    i = 0
    while i < len(text):
        marker = text[i]
        if marker in "^_":
            run, i = _next_script_run(text, i)
            run_font = sup_font if marker == "^" else sub_font
            offset_y = -round(_font_size(fnt) * 0.22) if marker == "^" else round(_font_size(fnt) * 0.36)
            draw.text((x, y + offset_y), run, font=run_font, fill=fill)
            x += text_box(draw, run, run_font)[0]
            continue
        start = i
        while i < len(text) and text[i] not in "^_":
            i += 1
        run = text[start:i]
        draw.text((x, y), run, font=fnt, fill=fill)
        x += text_box(draw, run, fnt)[0]
    return x, y + rich_text_box(draw, text, fnt, sup_font, sub_font)[1]


def open_rgba(path: Path) -> Image.Image:
    return Image.open(path).convert("RGBA")


def crop_visible(img: Image.Image, white_threshold: int = 248) -> Image.Image:
    rgba = img.convert("RGBA")
    px = rgba.load()
    min_x, min_y, max_x, max_y = rgba.width, rgba.height, -1, -1
    for y in range(rgba.height):
        for x in range(rgba.width):
            r, g, b, a = px[x, y]
            visible = a > 8 and not (r >= white_threshold and g >= white_threshold and b >= white_threshold)
            if visible:
                min_x = min(min_x, x)
                min_y = min(min_y, y)
                max_x = max(max_x, x)
                max_y = max(max_y, y)
    if max_x < min_x:
        return rgba
    pad = 8
    return rgba.crop(
        (
            max(0, min_x - pad),
            max(0, min_y - pad),
            min(rgba.width, max_x + pad + 1),
            min(rgba.height, max_y + pad + 1),
        )
    )


def fit(img: Image.Image, max_w: int, max_h: int) -> Image.Image:
    scale = min(max_w / img.width, max_h / img.height)
    size = (max(1, int(img.width * scale)), max(1, int(img.height * scale)))
    return img.resize(size, Image.Resampling.LANCZOS)


def cover_crop(img: Image.Image, size: tuple[int, int], focus: tuple[float, float] = (0.5, 0.5)) -> Image.Image:
    target_w, target_h = size
    scale = max(target_w / img.width, target_h / img.height)
    resized = img.resize((math.ceil(img.width * scale), math.ceil(img.height * scale)), Image.Resampling.LANCZOS)
    max_x = resized.width - target_w
    max_y = resized.height - target_h
    left = int(max(0, min(max_x, max_x * focus[0])))
    top = int(max(0, min(max_y, max_y * focus[1])))
    return resized.crop((left, top, left + target_w, top + target_h))


def rounded_mask(size: tuple[int, int], radius: int) -> Image.Image:
    mask = Image.new("L", size, 0)
    draw = ImageDraw.Draw(mask)
    draw.rounded_rectangle((0, 0, size[0], size[1]), radius=radius, fill=255)
    return mask


def paste_rounded(base: Image.Image, img: Image.Image, box: tuple[int, int, int, int], radius: int = 12) -> None:
    x0, y0, x1, y1 = box
    crop = cover_crop(img.convert("RGBA"), (x1 - x0, y1 - y0))
    mask = rounded_mask(crop.size, radius)
    base.paste(crop, (x0, y0), mask)


def paste_fit(base: Image.Image, img: Image.Image, box: tuple[int, int, int, int], anchor: str = "center") -> None:
    x0, y0, x1, y1 = box
    fitted = fit(img, x1 - x0, y1 - y0)
    if anchor == "left":
        x = x0
    elif anchor == "right":
        x = x1 - fitted.width
    else:
        x = x0 + ((x1 - x0) - fitted.width) // 2
    y = y0 + ((y1 - y0) - fitted.height) // 2
    base.alpha_composite(fitted, (x, y))


def paste_fit_return_box(
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


def load_sphenix_logo() -> Image.Image | None:
    logo = TITLE_ASSET_DIR / "sphenix-logo-white-bg_0.png"
    if not logo.exists():
        return None
    return crop_visible(open_rgba(logo), white_threshold=252)


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    text: str,
    xy: tuple[int, int],
    max_width: int,
    fnt: ImageFont.ImageFont,
    fill: tuple[int, int, int],
    line_gap: int,
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


def alpha_layer(base: Image.Image, callback) -> None:
    layer = Image.new("RGBA", base.size, (0, 0, 0, 0))
    callback(ImageDraw.Draw(layer, "RGBA"))
    base.alpha_composite(layer)


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


def draw_polyline(
    draw: ImageDraw.ImageDraw,
    pts: list[tuple[float, float]],
    fill: tuple[int, int, int, int] | tuple[int, int, int],
    width: int,
) -> None:
    draw.line([(round(x), round(y)) for x, y in pts], fill=fill, width=width, joint="curve")


def interp(a: tuple[float, float], b: tuple[float, float], t: float) -> tuple[float, float]:
    return (a[0] * (1 - t) + b[0] * t, a[1] * (1 - t) + b[1] * t)


def draw_detector_face(base: Image.Image, impact: tuple[int, int]) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    tl = (1148, 392)
    tr = (1358, 454)
    br = (1304, 664)
    bl = (1094, 604)
    draw.polygon([tl, tr, br, bl], fill=(237, 244, 249, 255), outline=(172, 194, 212, 255))

    for i in range(1, 5):
        t = i / 5
        p1 = interp(tl, bl, t)
        p2 = interp(tr, br, t)
        draw.line((p1, p2), fill=(197, 211, 223, 255), width=3)
    for i in range(1, 6):
        t = i / 6
        p1 = interp(tl, tr, t)
        p2 = interp(bl, br, t)
        draw.line((p1, p2), fill=(197, 211, 223, 255), width=3)

    glow = Image.new("RGBA", base.size, (0, 0, 0, 0))
    gdraw = ImageDraw.Draw(glow, "RGBA")
    x, y = impact
    for r, alpha in ((74, 42), (50, 72), (28, 110)):
        gdraw.ellipse((x - r, y - r, x + r, y + r), fill=(*PHOTON, alpha))
    glow = glow.filter(ImageFilter.GaussianBlur(8))
    base.alpha_composite(glow)
    draw = ImageDraw.Draw(base, "RGBA")
    draw.ellipse((x - 22, y - 22, x + 22, y + 22), fill=(*PHOTON, 230), outline=(*PHOTON_DARK, 230), width=3)
    draw.text((1138, 686), "EMCal", font=font(TIMES_ITALIC, 31), fill=LIGHT_MUTED)


def draw_icon(draw: ImageDraw.ImageDraw, kind: int, center: tuple[int, int]) -> None:
    cx, cy = center
    if kind == 0:
        pts = feynman_points((cx - 34, cy + 12), (cx + 34, cy - 12), 7, 3.6, 80)
        draw_polyline(draw, pts, (*PHOTON_DARK, 255), 4)
        draw.ellipse((cx - 44, cy - 44, cx + 44, cy + 44), outline=(195, 207, 220, 255), width=3)
        draw.line((cx - 52, cy, cx - 34, cy), fill=(195, 207, 220, 255), width=3)
        draw.line((cx + 34, cy, cx + 52, cy), fill=(195, 207, 220, 255), width=3)
    elif kind == 1:
        draw.pieslice((cx - 47, cy - 47, cx + 47, cy + 47), start=320, end=40, fill=(250, 223, 154, 120), outline=(*PHOTON_DARK, 200), width=3)
        draw.line((cx, cy, cx + 46, cy - 17), fill=(*PHOTON_DARK, 220), width=3)
        draw.line((cx, cy, cx + 46, cy + 17), fill=(*PHOTON_DARK, 220), width=3)
        draw.ellipse((cx - 10, cy - 10, cx + 10, cy + 10), fill=(*PHOTON, 230))
    else:
        draw.line((cx - 42, cy + 34, cx + 44, cy + 34), fill=(73, 104, 128, 255), width=4)
        draw.line((cx - 42, cy + 34, cx - 42, cy - 34), fill=(73, 104, 128, 255), width=4)
        draw.line((cx - 30, cy + 22, cx + 32, cy - 18), fill=(*TEAL, 230), width=5)
        draw.ellipse((cx + 22, cy - 26, cx + 42, cy - 6), fill=(*SPHENIX_BLUE, 220))


def draw_context_icon(draw: ImageDraw.ImageDraw, kind: int, center: tuple[int, int]) -> None:
    cx, cy = center
    if kind == 0:
        draw.line((cx - 46, cy + 36, cx + 48, cy + 36), fill=(72, 101, 122, 255), width=4)
        draw.line((cx - 46, cy + 36, cx - 46, cy - 36), fill=(72, 101, 122, 255), width=4)
        draw.line((cx - 35, cy + 18, cx - 5, cy - 8, cx + 38, cy - 24), fill=(*SPHENIX_BLUE, 230), width=5)
        draw.ellipse((cx + 28, cy - 34, cx + 50, cy - 12), fill=(*SPHENIX_BLUE, 230))
    elif kind == 1:
        pts = feynman_points((cx - 50, cy + 4), (cx + 50, cy - 4), 8, 5.0, 100)
        draw_polyline(draw, pts, (*PHOTON_DARK, 255), 5)
        draw.pieslice((cx - 56, cy - 56, cx + 56, cy + 56), start=330, end=30, outline=(*PHOTON_DARK, 190), width=4)
        draw.line((cx, cy, cx + 54, cy - 17), fill=(*PHOTON_DARK, 175), width=3)
        draw.line((cx, cy, cx + 54, cy + 17), fill=(*PHOTON_DARK, 175), width=3)
    else:
        for r, color in ((49, TEAL), (34, PHOTON), (18, SPHENIX_BLUE)):
            draw.ellipse((cx - r, cy - r, cx + r, cy + r), outline=(*color, 230), width=5)
        draw.line((cx, cy, cx + 48, cy - 22), fill=(*PHOTON_DARK, 240), width=5)
        draw.ellipse((cx + 40, cy - 30, cx + 58, cy - 12), fill=(*PHOTON, 230))


def draw_segmented_ring(
    draw: ImageDraw.ImageDraw,
    center: tuple[int, int],
    r_inner: int,
    r_outer: int,
    fill: tuple[int, int, int],
    outline: tuple[int, int, int],
    start_offset: int = 0,
    highlight: tuple[int, int] | None = None,
) -> None:
    cx, cy = center
    for i in range(16):
        start = start_offset + i * 22
        end = start + 16
        color = fill
        alpha = 72
        if highlight is not None and highlight[0] <= i <= highlight[1]:
            color = PHOTON
            alpha = 150
        box_outer = (cx - r_outer, cy - r_outer, cx + r_outer, cy + r_outer)
        draw.pieslice(box_outer, start=start, end=end, fill=(*color, alpha))
    draw.ellipse((cx - r_outer, cy - r_outer, cx + r_outer, cy + r_outer), outline=(*outline, 210), width=4)
    draw.ellipse((cx - r_inner, cy - r_inner, cx + r_inner, cy + r_inner), outline=(*outline, 155), width=3)


def draw_context_visual(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    panel = (126, 332, 1434, 1128)
    draw.rounded_rectangle(panel, radius=10, fill=(*PANEL, 255), outline=(*PANEL_EDGE, 255), width=2)

    # Subtle RHIC program cue.
    draw.ellipse((230, 430, 1240, 1014), outline=(200, 214, 226, 150), width=6)
    draw.arc((230, 430, 1240, 1014), start=202, end=338, fill=(*SPHENIX_BLUE, 155), width=8)
    draw.arc((230, 430, 1240, 1014), start=24, end=158, fill=(*PHOTON_DARK, 140), width=8)
    draw.text((240, 388), "RHIC p+p at 200 GeV", font=font(TIMES_ITALIC, 34), fill=MUTED)
    draw.text((1018, 976), "hard-probes reference", font=font(TIMES_ITALIC, 31), fill=LIGHT_MUTED)

    center = (704, 690)
    draw.line((312, center[1], 1096, center[1]), fill=(166, 191, 211, 255), width=6)
    draw.ellipse((center[0] - 19, center[1] - 19, center[0] + 19, center[1] + 19), fill=(*SPHENIX_BLUE, 230))
    draw.text((336, center[1] - 52), "p", font=font(TIMES_BOLD, 38), fill=BLUE)
    draw.text((1036, center[1] - 52), "p", font=font(TIMES_BOLD, 38), fill=BLUE)

    draw_segmented_ring(draw, center, 82, 120, (126, 169, 198), (83, 124, 153), start_offset=10)
    draw_segmented_ring(draw, center, 148, 210, (245, 181, 34), (183, 137, 39), start_offset=2, highlight=(1, 3))
    draw_segmented_ring(draw, center, 230, 302, (18, 112, 142), (18, 96, 124), start_offset=12)
    draw.ellipse((center[0] - 330, center[1] - 330, center[0] + 330, center[1] + 330), outline=(57, 83, 112, 130), width=4)
    draw.text((center[0] - 98, center[1] - 22), "sPHENIX", font=font(TIMES_BOLD, 39), fill=BLUE)

    photon_start = (center[0] + 28, center[1] - 22)
    photon_end = (center[0] + 205, center[1] - 112)
    cone_left = (center[0] + 297, center[1] - 190)
    cone_right = (center[0] + 294, center[1] - 42)
    alpha_layer(base, lambda d: d.polygon([photon_start, cone_left, cone_right], fill=(245, 181, 34, 26)))
    draw.line((photon_start, cone_left), fill=(219, 157, 32, 90), width=3)
    draw.line((photon_start, cone_right), fill=(219, 157, 32, 90), width=3)
    wave = feynman_points(photon_start, photon_end, 13, 4.8, 180)
    glow = Image.new("RGBA", base.size, (0, 0, 0, 0))
    gdraw = ImageDraw.Draw(glow, "RGBA")
    draw_polyline(gdraw, wave, (*PHOTON, 115), 16)
    glow = glow.filter(ImageFilter.GaussianBlur(6))
    base.alpha_composite(glow)
    draw = ImageDraw.Draw(base, "RGBA")
    draw_polyline(draw, wave, (*PHOTON, 255), 8)
    draw_polyline(draw, wave, (*PHOTON_DARK, 210), 3)
    draw.ellipse((photon_end[0] - 19, photon_end[1] - 19, photon_end[0] + 19, photon_end[1] + 19), fill=(*PHOTON, 240), outline=(*PHOTON_DARK, 230), width=3)

    for angle, length, width, alpha in ((38, 166, 7, 105), (318, 142, 6, 85), (285, 160, 5, 75)):
        ex = center[0] + math.cos(math.radians(angle)) * length
        ey = center[1] + math.sin(math.radians(angle)) * length
        draw.line((center[0], center[1], ex, ey), fill=(*TEAL_SOFT, alpha), width=width)
        draw.ellipse((ex - 9, ey - 9, ex + 9, ey + 9), fill=(*TEAL_SOFT, alpha + 25))

    legend_x, legend_y = 1048, 500
    legend = [
        ("Tracking", (126, 169, 198)),
        ("EMCal", PHOTON),
        ("HCal", TEAL),
    ]
    for i, (label, color) in enumerate(legend):
        yy = legend_y + i * 58
        draw.rounded_rectangle((legend_x, yy, legend_x + 34, yy + 22), radius=4, fill=(*color, 190))
        draw.text((legend_x + 50, yy - 6), label, font=font(TIMES, 34), fill=MUTED)
    draw.text((898, 454), "isolated photon candidate", font=font(TIMES_ITALIC, 34), fill=BLUE)

    # Measurement chain: the whole talk in one quiet line.
    chain_y = 1032
    steps = [
        ("p+p", "reference"),
        ("EMCal", "clusters"),
        ("isolated", "photons"),
        ("corrected", "cross section"),
    ]
    x = 202
    for i, (top, bottom) in enumerate(steps):
        w = (190, 220, 230, 252)[i]
        draw.rounded_rectangle((x, chain_y, x + w, chain_y + 72), radius=8, fill=(255, 255, 255, 255), outline=(216, 225, 234, 255), width=2)
        top_w, _ = text_box(draw, top, font(TIMES_BOLD, 29))
        bottom_w, _ = text_box(draw, bottom, font(TIMES, 25))
        draw.text((x + (w - top_w) / 2, chain_y + 9), top, font=font(TIMES_BOLD, 29), fill=INK)
        draw.text((x + (w - bottom_w) / 2, chain_y + 42), bottom, font=font(TIMES, 25), fill=MUTED)
        if i < len(steps) - 1:
            draw.line((x + w + 18, chain_y + 36, x + w + 72, chain_y + 36), fill=(158, 178, 194, 255), width=3)
            draw.polygon([(x + w + 72, chain_y + 36), (x + w + 58, chain_y + 27), (x + w + 58, chain_y + 45)], fill=(158, 178, 194, 255))
        x += w + 78


def draw_context_cards(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x, y = 1510, 338
    w, h, gap = 870, 224, 32
    for idx, (heading, body) in enumerate(CONTEXT_REASONS):
        top = y + idx * (h + gap)
        stripe = SPHENIX_BLUE if idx == 0 else PHOTON if idx == 1 else TEAL
        draw.rounded_rectangle((x, top, x + w, top + h), radius=10, fill=(*CARD, 255), outline=(*CARD_EDGE, 255), width=2)
        draw.rounded_rectangle((x, top, x + 14, top + h), radius=6, fill=(*stripe, 255))
        draw_context_icon(draw, idx, (x + 86, top + 112))
        draw.text((x + 162, top + 38), heading, font=font(TIMES_BOLD, 41), fill=INK)
        draw_wrapped(draw, body, (x + 164, top + 98), 632, font(TIMES, 32), fill=MUTED, line_gap=10)


def load_real_detector_asset() -> tuple[Image.Image, str]:
    if REAL_DETECTOR_PHOTO.exists():
        return open_rgba(REAL_DETECTOR_PHOTO), str(REAL_DETECTOR_PHOTO.relative_to(ROOT))
    if REAL_DETECTOR_RENDERING.exists():
        return open_rgba(REAL_DETECTOR_RENDERING), str(REAL_DETECTOR_RENDERING.relative_to(ROOT))
    raise FileNotFoundError(
        f"missing official BNL detector asset: expected {REAL_DETECTOR_PHOTO} or {REAL_DETECTOR_RENDERING}"
    )


def enhance_detector_photo(img: Image.Image) -> Image.Image:
    rgb = img.convert("RGB")
    rgb = ImageEnhance.Brightness(rgb).enhance(1.05)
    rgb = ImageEnhance.Color(rgb).enhance(0.98)
    rgb = ImageEnhance.Contrast(rgb).enhance(1.05)
    rgb = ImageEnhance.Sharpness(rgb).enhance(1.05)
    return rgb.convert("RGBA")


def overlay_photo_gradient(base: Image.Image, box: tuple[int, int, int, int], radius: int) -> None:
    x0, y0, x1, y1 = box
    w, h = x1 - x0, y1 - y0
    grad = Image.new("RGBA", (w, h), (0, 0, 0, 0))
    px = grad.load()
    for y in range(h):
        y_alpha = int(20 * (1 - abs((y / max(1, h - 1)) - 0.48)))
        bottom_alpha = int(38 * max(0, (y - h * 0.54) / (h * 0.46)))
        for x in range(w):
            left_alpha = int(28 * max(0, 1 - x / (w * 0.36)))
            px[x, y] = (10, 20, 34, max(y_alpha, bottom_alpha, left_alpha))
    alpha = ImageChops.multiply(grad.getchannel("A"), rounded_mask((w, h), radius))
    grad.putalpha(alpha)
    base.alpha_composite(grad, (x0, y0))


def draw_real_detector_photo_panel(base: Image.Image, detector_img: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    panel = (132, 328, 2428, 872)
    radius = 12

    shadow = Image.new("RGBA", base.size, (0, 0, 0, 0))
    sdraw = ImageDraw.Draw(shadow, "RGBA")
    sdraw.rounded_rectangle((panel[0] + 8, panel[1] + 10, panel[2] + 8, panel[3] + 10), radius=radius, fill=(30, 42, 58, 38))
    shadow = shadow.filter(ImageFilter.GaussianBlur(12))
    base.alpha_composite(shadow)

    photo = enhance_detector_photo(detector_img)
    crop = cover_crop(photo, (panel[2] - panel[0], panel[3] - panel[1]), focus=(0.54, 0.48))
    base.paste(crop, (panel[0], panel[1]), rounded_mask(crop.size, radius))
    overlay_photo_gradient(base, panel, radius)
    draw.rounded_rectangle(panel, radius=radius, outline=(216, 224, 232, 255), width=2)

    label = "sPHENIX detector at RHIC"
    label_font = font(TIMES_BOLD, 34)
    label_w, label_h = text_box(draw, label, label_font)
    label_box = (172, 360, 172 + label_w + 48, 360 + label_h + 34)
    draw.rounded_rectangle(label_box, radius=8, fill=(255, 255, 255, 218), outline=(222, 228, 235, 210), width=1)
    draw.text((196, 374), label, font=label_font, fill=BLUE)

    # These overlays are intentionally non-technical: they cue "photon through
    # detector" without asserting an exact subsystem boundary on the photo.
    highlight = Image.new("RGBA", base.size, (0, 0, 0, 0))
    hdraw = ImageDraw.Draw(highlight, "RGBA")
    hdraw.ellipse((1222, 406, 1966, 784), outline=(*PHOTON, 90), width=22)
    hdraw.ellipse((1234, 418, 1954, 772), outline=(*SPHENIX_BLUE, 66), width=7)
    highlight = highlight.filter(ImageFilter.GaussianBlur(4))
    base.alpha_composite(highlight)
    draw = ImageDraw.Draw(base, "RGBA")
    draw.ellipse((1230, 414, 1958, 776), outline=(*PHOTON, 172), width=4)

    start = (1088, 624)
    impact = (1568, 566)
    wave = feynman_points(start, impact, amplitude=12, cycles=3.6, steps=180)
    glow = Image.new("RGBA", base.size, (0, 0, 0, 0))
    gdraw = ImageDraw.Draw(glow, "RGBA")
    draw_polyline(gdraw, wave, (*PHOTON, 125), 18)
    glow = glow.filter(ImageFilter.GaussianBlur(8))
    base.alpha_composite(glow)
    draw = ImageDraw.Draw(base, "RGBA")
    draw_polyline(draw, wave, (*PHOTON, 240), 8)
    draw_polyline(draw, wave, (*PHOTON_DARK, 215), 3)
    draw.ellipse((impact[0] - 18, impact[1] - 18, impact[0] + 18, impact[1] + 18), fill=(*PHOTON, 225), outline=(255, 240, 190, 230), width=2)

    note = "isolated photon measurement enters as the p+p reference"
    note_font = font(TIMES_ITALIC, 34)
    note_w, note_h = text_box(draw, note, note_font)
    note_box = (2428 - note_w - 70, 798, 2428 - 28, 798 + note_h + 32)
    draw.rounded_rectangle(note_box, radius=8, fill=(5, 16, 29, 150), outline=(255, 255, 255, 70), width=1)
    draw.text((note_box[0] + 22, note_box[1] + 13), note, font=note_font, fill=(255, 255, 255, 238))


def draw_real_detector_cards(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x0, y0 = 132, 928
    gap = 40
    w = (2428 - 132 - 2 * gap) // 3
    h = 242
    accents = [SPHENIX_BLUE, PHOTON, TEAL]
    for idx, (heading, body) in enumerate(REAL_DETECTOR_REASONS):
        x = x0 + idx * (w + gap)
        draw.rounded_rectangle((x, y0, x + w, y0 + h), radius=8, fill=(*CARD, 255), outline=(220, 227, 235, 255), width=2)
        draw.rounded_rectangle((x, y0, x + w, y0 + 12), radius=6, fill=(*accents[idx], 255))
        draw_context_icon(draw, idx, (x + 76, y0 + 126))
        draw.text((x + 146, y0 + 52), heading, font=font(TIMES_BOLD, 39), fill=INK)
        draw_wrapped(draw, body, (x + 146, y0 + 108), w - 186, font(TIMES, 30), fill=MUTED, line_gap=9)


def load_detector_rendering_asset() -> tuple[Image.Image, str]:
    if REAL_DETECTOR_RENDERING.exists():
        return open_rgba(REAL_DETECTOR_RENDERING), str(REAL_DETECTOR_RENDERING.relative_to(ROOT))
    if REAL_DETECTOR_PHOTO.exists():
        return open_rgba(REAL_DETECTOR_PHOTO), str(REAL_DETECTOR_PHOTO.relative_to(ROOT))
    raise FileNotFoundError(
        f"missing official BNL detector asset: expected {REAL_DETECTOR_RENDERING} or {REAL_DETECTOR_PHOTO}"
    )


def draw_callout_stack(base: Image.Image, x: int, y: int, w: int) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    for idx, (label, body, color) in enumerate(DETECTOR_CALLOUTS):
        top = y + idx * 104
        draw.rounded_rectangle((x, top, x + w, top + 78), radius=8, fill=(255, 255, 255, 245), outline=(218, 226, 235, 255), width=2)
        draw.rounded_rectangle((x, top, x + 12, top + 78), radius=5, fill=(*color, 245))
        draw.ellipse((x + 30, top + 23, x + 58, top + 51), fill=(*color, 170))
        draw.text((x + 78, top + 13), label, font=font(TIMES_BOLD, 29), fill=INK)
        draw_wrapped(draw, body, (x + 78, top + 45), w - 94, font(TIMES, 22), fill=MUTED, line_gap=2)


def draw_detector_data_image_panel(base: Image.Image, detector_img: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    panel = (132, 330, 1438, 1032)
    draw.rounded_rectangle(panel, radius=10, fill=(255, 255, 255, 255), outline=(216, 225, 234, 255), width=2)

    image_box = (170, 372, 1122, 948)
    draw.rounded_rectangle(image_box, radius=8, fill=(248, 250, 252, 255), outline=(226, 232, 238, 255), width=1)
    rendering = detector_img.convert("RGBA")
    rendering = ImageEnhance.Color(rendering).enhance(1.04)
    rendering = ImageEnhance.Contrast(rendering).enhance(1.03)
    rendering = ImageEnhance.Sharpness(rendering).enhance(1.08)
    paste_fit(base, rendering, image_box)

    draw = ImageDraw.Draw(base, "RGBA")
    label = "sPHENIX detector rendering"
    draw.text((186, 962), label, font=font(TIMES_ITALIC, 26), fill=LIGHT_MUTED)

    # Local visual emphasis only: a quiet prompt-photon path into the
    # calorimetry region without pretending to label exact CAD boundaries.
    start = (432, 640)
    impact = (943, 538)
    wave = feynman_points(start, impact, amplitude=11, cycles=5.4, steps=180)
    glow = Image.new("RGBA", base.size, (0, 0, 0, 0))
    gdraw = ImageDraw.Draw(glow, "RGBA")
    draw_polyline(gdraw, wave, (*PHOTON, 92), 16)
    glow = glow.filter(ImageFilter.GaussianBlur(7))
    base.alpha_composite(glow)
    draw = ImageDraw.Draw(base, "RGBA")
    draw_polyline(draw, wave, (*PHOTON, 218), 7)
    draw_polyline(draw, wave, (*PHOTON_DARK, 160), 2)
    for r, alpha in ((58, 34), (35, 64), (18, 150)):
        draw.ellipse((impact[0] - r, impact[1] - r, impact[0] + r, impact[1] + r), outline=(*PHOTON, alpha), width=4)

    draw_callout_stack(base, 1128, 398, 268)
    draw.text((170, 336), "Photon-relevant detector handles", font=font(TIMES_BOLD, 34), fill=BLUE)


def draw_detector_data_cards(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x, y = 1512, 330
    w, h, gap = 916, 166, 28
    accents = [SPHENIX_BLUE, PHOTON, TEAL]
    for idx, (heading, body) in enumerate(DETECTOR_DATA_CARDS):
        top = y + idx * (h + gap)
        draw.rounded_rectangle((x, top, x + w, top + h), radius=9, fill=(*CARD, 255), outline=(*CARD_EDGE, 255), width=2)
        draw.rounded_rectangle((x, top, x + 13, top + h), radius=5, fill=(*accents[idx], 255))
        draw_context_icon(draw, idx, (x + 76, top + 83))
        draw.text((x + 150, top + 34), heading, font=font(TIMES_BOLD, 38), fill=INK)
        draw_wrapped(draw, body, (x + 150, top + 86), 708, font(TIMES, 28), fill=MUTED, line_gap=8)


def draw_pp_data_panel(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    panel = (1512, 942, 2428, 1192)
    draw.rounded_rectangle(panel, radius=10, fill=(255, 255, 255, 255), outline=(216, 225, 234, 255), width=2)
    draw.text((1548, 974), "p+p data context for this talk", font=font(TIMES_BOLD, 35), fill=BLUE)

    chip_y = 1030
    chip_w = 398
    for idx, (year, line1, line2, color) in enumerate(PP_DATA_ROWS):
        x = 1548 + idx * (chip_w + 46)
        draw.rounded_rectangle((x, chip_y, x + chip_w, chip_y + 102), radius=8, fill=(247, 250, 252, 255), outline=(222, 229, 236, 255), width=2)
        draw.rounded_rectangle((x, chip_y, x + chip_w, chip_y + 32), radius=7, fill=(*color, 235))
        year_w, _ = text_box(draw, year, font(TIMES_BOLD, 25))
        draw.text((x + (chip_w - year_w) / 2, chip_y + 4), year, font=font(TIMES_BOLD, 25), fill=(255, 255, 255))
        draw.text((x + 24, chip_y + 47), line1, font=font(TIMES_BOLD, 25), fill=INK)
        draw.text((x + 24, chip_y + 76), line2, font=font(TIMES, 23), fill=MUTED)

    draw_wrapped(draw, MUTED_PROGRAM_CONTEXT, (1548, 1148), 820, font(TIMES_ITALIC, 23), fill=LIGHT_MUTED, line_gap=6)


def draw_progressive_header(base: Image.Image, subtitle: str, active_step: str) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*PHOTON, 255))

    logo = load_sphenix_logo()
    if logo is not None:
        paste_fit(base, logo, (2090, 78, 2388, 160), anchor="right")

    draw.text((132, 84), PROGRESSIVE_TITLE, font=font(TIMES_BOLD, 72), fill=INK)
    draw.text((136, 182), subtitle, font=font(TIMES_ITALIC, 36), fill=BLUE)

    tabs = [
        ("Detector systems", "detector"),
        ("Data context", "data"),
    ]
    x = 136
    for label, key in tabs:
        active = key == active_step
        label_font = font(TIMES_BOLD if active else TIMES, 26)
        tw, th = text_box(draw, label, label_font)
        fill = (*SPHENIX_BLUE, 235) if active else (239, 246, 250, 255)
        text_fill = (255, 255, 255) if active else BLUE
        outline = (*SPHENIX_BLUE, 235) if active else (211, 224, 234, 255)
        draw.rounded_rectangle((x, 242, x + tw + 44, 284), radius=20, fill=fill, outline=outline, width=2)
        draw.text((x + 22, 250), label, font=label_font, fill=text_fill)
        x += tw + 64
    draw.line((132, 306, W - 132, 306), fill=(221, 226, 232), width=3)


def draw_leader(
    draw: ImageDraw.ImageDraw,
    start: tuple[int, int],
    target: tuple[int, int],
    color: tuple[int, int, int],
    side: str,
    approximate: bool = False,
) -> None:
    bend_x = start[0] + (94 if side == "left" else -94)
    pts = [start, (bend_x, start[1]), (bend_x, target[1]), target]
    line_alpha = 185 if approximate else 235
    dot_alpha = 225 if approximate else 255
    halo_alpha = 105 if approximate else 155
    line_w = 6 if approximate else 7
    halo_w = line_w + 8
    draw.line(pts, fill=(255, 255, 255, 235), width=halo_w, joint="curve")
    draw.line(pts, fill=(*color, line_alpha), width=line_w, joint="curve")

    arrow_poly = None
    arrow_outline_poly = None
    if len(pts) >= 2:
        prev = pts[-2]
        dx = target[0] - prev[0]
        dy = target[1] - prev[1]
        if dx or dy:
            angle = math.atan2(dy, dx)
            size = 25 if not approximate else 21
            spread = 0.58
            tip = (
                target[0] - math.cos(angle) * (38 if not approximate else 32),
                target[1] - math.sin(angle) * (38 if not approximate else 32),
            )
            base_x = tip[0] - math.cos(angle) * size
            base_y = tip[1] - math.sin(angle) * size
            p1 = (
                base_x + math.cos(angle + math.pi / 2) * size * spread,
                base_y + math.sin(angle + math.pi / 2) * size * spread,
            )
            p2 = (
                base_x + math.cos(angle - math.pi / 2) * size * spread,
                base_y + math.sin(angle - math.pi / 2) * size * spread,
            )
            arrow_poly = [tip, p1, p2]
            outline_size = size + 8
            outline_base_x = tip[0] - math.cos(angle) * outline_size
            outline_base_y = tip[1] - math.sin(angle) * outline_size
            arrow_outline_poly = [
                tip,
                (
                    outline_base_x + math.cos(angle + math.pi / 2) * outline_size * spread,
                    outline_base_y + math.sin(angle + math.pi / 2) * outline_size * spread,
                ),
                (
                    outline_base_x + math.cos(angle - math.pi / 2) * outline_size * spread,
                    outline_base_y + math.sin(angle - math.pi / 2) * outline_size * spread,
                ),
            ]

    halo_r = 42 if not approximate else 36
    ring_r = 26 if not approximate else 22
    dot_r = 12 if not approximate else 10
    draw.ellipse(
        (target[0] - halo_r, target[1] - halo_r, target[0] + halo_r, target[1] + halo_r),
        fill=(*color, 38 if not approximate else 28),
        outline=(*color, halo_alpha),
        width=4,
    )
    draw.ellipse(
        (target[0] - ring_r, target[1] - ring_r, target[0] + ring_r, target[1] + ring_r),
        fill=(255, 255, 255, 222),
        outline=(*color, dot_alpha),
        width=5,
    )
    draw.ellipse(
        (target[0] - dot_r, target[1] - dot_r, target[0] + dot_r, target[1] + dot_r),
        fill=(*color, dot_alpha),
    )
    if arrow_poly and arrow_outline_poly:
        draw.polygon(arrow_outline_poly, fill=(255, 255, 255, 238))
        draw.polygon(arrow_poly, fill=(*color, line_alpha))


def draw_component_callout(
    draw: ImageDraw.ImageDraw,
    image_bounds: tuple[int, int, int, int],
    label: str,
    detail: str,
    color: tuple[int, int, int],
    target_rel: tuple[float, float],
    box_xy: tuple[int, int],
) -> None:
    ix0, iy0, ix1, iy1 = image_bounds
    target = (
        round(ix0 + target_rel[0] * (ix1 - ix0)),
        round(iy0 + target_rel[1] * (iy1 - iy0)),
    )
    x, y = box_xy
    w, h = 362, 88
    side = "left" if x < target[0] else "right"
    start = (x + w, y + h // 2) if side == "left" else (x, y + h // 2)
    draw_leader(draw, start, target, color, side)

    draw.rounded_rectangle((x, y, x + w, y + h), radius=8, fill=(255, 255, 255, 236), outline=(214, 224, 234, 255), width=2)
    draw.rounded_rectangle((x, y, x + 12, y + h), radius=5, fill=(*color, 245))
    draw.ellipse((x + 28, y + 30, x + 56, y + 58), fill=(*color, 185))
    draw.text((x + 76, y + 14), label, font=font(TIMES_BOLD, 28), fill=INK)
    draw.text((x + 76, y + 48), detail, font=font(TIMES, 22), fill=MUTED)


def draw_progressive_detector_map(base: Image.Image, detector_img: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    panel = (132, 340, 1468, 1088)
    image_box = (164, 414, 1434, 1032)

    shadow = Image.new("RGBA", base.size, (0, 0, 0, 0))
    sdraw = ImageDraw.Draw(shadow, "RGBA")
    sdraw.rounded_rectangle((panel[0] + 6, panel[1] + 8, panel[2] + 6, panel[3] + 8), radius=12, fill=(30, 42, 58, 34))
    shadow = shadow.filter(ImageFilter.GaussianBlur(12))
    base.alpha_composite(shadow)

    draw.rounded_rectangle(panel, radius=12, fill=(255, 255, 255, 255), outline=(216, 225, 234, 255), width=2)
    draw.text((168, 360), "Detector systems used by the photon measurement", font=font(TIMES_BOLD, 34), fill=BLUE)

    draw.rounded_rectangle(image_box, radius=10, fill=(248, 250, 252, 255), outline=(224, 231, 238, 255), width=1)
    rendering = detector_img.convert("RGBA")
    rendering = ImageEnhance.Color(rendering).enhance(1.05)
    rendering = ImageEnhance.Contrast(rendering).enhance(1.04)
    rendering = ImageEnhance.Sharpness(rendering).enhance(1.10)
    image_bounds = paste_fit_return_box(base, rendering, image_box)

    overlay = Image.new("RGBA", base.size, (0, 0, 0, 0))
    odraw = ImageDraw.Draw(overlay, "RGBA")
    for label, detail, color, target, box_xy in PROGRESSIVE_SYSTEM_CALLOUTS:
        draw_component_callout(odraw, image_bounds, label, detail, color, target, box_xy)
    base.alpha_composite(overlay)

    draw = ImageDraw.Draw(base, "RGBA")

def draw_progressive_detector_summary(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    panel = (1512, 340, 2428, 1088)
    draw.rounded_rectangle(panel, radius=12, fill=(255, 255, 255, 255), outline=(216, 225, 234, 255), width=2)
    draw.text((1550, 374), "Detector handles for this measurement", font=font(TIMES_BOLD, 40), fill=BLUE)
    draw.line((1550, 436, 2388, 436), fill=(222, 228, 235), width=2)

    cards = [
        ("EMCal measures the object", "cluster energy and shower shape are the first photon handles", PHOTON),
        ("Calorimetry defines isolation", "nearby EMCal/HCal activity tests whether the candidate is quiet", TEAL),
        ("Tracking and MBD anchor the event", "vertex, trigger, timing, and luminosity make the cross section possible", SPHENIX_BLUE),
    ]
    y = 478
    for heading, body, color in cards:
        draw.rounded_rectangle((1550, y, 2388, y + 156), radius=9, fill=(247, 250, 252, 255), outline=(222, 229, 236, 255), width=2)
        draw.rounded_rectangle((1550, y, 1564, y + 156), radius=5, fill=(*color, 255))
        draw.text((1594, y + 26), heading, font=font(TIMES_BOLD, 34), fill=INK)
        draw_wrapped(draw, body, (1594, y + 78), 742, font(TIMES, 27), fill=MUTED, line_gap=7)
        y += 186

    note = "Experiment context first; photon selection is introduced after the data sample."
    draw.rounded_rectangle((1550, 1000, 2388, 1052), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=1)
    draw.text((1580, 1014), note, font=font(TIMES_ITALIC, 24), fill=BLUE)


def draw_number_badge(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], value: str, label: str, color: tuple[int, int, int]) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=10, fill=(247, 250, 252, 255), outline=(222, 229, 236, 255), width=2)
    draw.rounded_rectangle((x0, y0, x1, y0 + 12), radius=6, fill=(*color, 255))
    value_font = font(TIMES_BOLD, 43)
    vw, _ = text_box(draw, value, value_font)
    draw.text((x0 + (x1 - x0 - vw) / 2, y0 + 30), value, font=value_font, fill=INK)
    label_font = font(TIMES, 23)
    if label.startswith(("pb^-1", "nb^-1")) and " " in label:
        unit, rest = label.split(" ", 1)
        draw_rich_text(draw, (x0 + 22, y0 + 86), unit, label_font, MUTED)
        draw_wrapped(draw, rest, (x0 + 22, y0 + 112), x1 - x0 - 44, label_font, fill=MUTED, line_gap=4)
    else:
        draw_wrapped(draw, label, (x0 + 22, y0 + 86), x1 - x0 - 44, label_font, fill=MUTED, line_gap=4)


def draw_progressive_data_panel(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    panel = (1512, 340, 2428, 1088)
    draw.rounded_rectangle(panel, radius=12, fill=(255, 255, 255, 255), outline=(216, 225, 234, 255), width=2)
    draw.text((1550, 374), "Data context for this talk", font=font(TIMES_BOLD, 40), fill=BLUE)
    draw_wrapped(
        draw,
        "Available p+p data set the context; the PPG12 result uses a defined Run 24 analysis sample.",
        (1552, 426),
        790,
        font(TIMES_ITALIC, 23),
        fill=MUTED,
        line_gap=4,
    )

    draw.rounded_rectangle((1550, 496, 2388, 734), radius=10, fill=(250, 252, 254, 255), outline=(221, 228, 236, 255), width=2)
    draw.text((1582, 520), "p+p running", font=font(TIMES_BOLD, 34), fill=INK)
    draw_number_badge(draw, (1582, 568, 1818, 710), "107", "pb^-1 2024 calorimeter-only", SPHENIX_BLUE)
    draw_number_badge(draw, (1844, 568, 2100, 710), "13", "pb^-1 2024 all subsystems", PHOTON)
    draw_number_badge(draw, (2126, 568, 2356, 710), "17", "pb^-1 2026 all subsystems", PHOTON_DARK)

    draw.rounded_rectangle((1550, 760, 2388, 858), radius=10, fill=(239, 246, 250, 255), outline=(208, 222, 233, 255), width=2)
    draw.rounded_rectangle((1550, 760, 1566, 858), radius=6, fill=(*TEAL, 255))
    draw.text((1594, 780), "PPG12 analysis sample is separate", font=font(TIMES_BOLD, 30), fill=INK)
    draw_rich_text(draw, (1594, 820), "Run 24 isolated prompt-photon result: L = 64.4 pb^-1", font(TIMES, 29), BLUE)

    draw.text((1550, 900), "Muted future-program context", font=font(TIMES_BOLD, 30), fill=MUTED)
    side_cards = [
        ("2025 Au+Au", "6.6 nb^-1 all subsystems"),
        ("2026 O+O", "52.6 nb^-1 calo-only; 23.6 nb^-1 all subsystems"),
    ]
    y = 944
    for heading, body in side_cards:
        draw.rounded_rectangle((1550, y, 2388, y + 58), radius=8, fill=(249, 250, 251, 255), outline=(224, 230, 236, 255), width=1)
        draw.text((1580, y + 11), heading, font=font(TIMES_BOLD, 23), fill=LIGHT_MUTED)
        draw_rich_text(draw, (1772, y + 12), body, font(TIMES, 23), LIGHT_MUTED)
        y += 70


def draw_progressive_bridge(base: Image.Image, text: str) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    draw.rounded_rectangle((132, 1228, W - 132, 1318), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw.text((174, 1253), text, font=font(TIMES_ITALIC, 38), fill=BLUE)


def render_detector_data_progressive_frame(output_dir: Path, frame: str) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    detector_img, _detector_asset = load_detector_rendering_asset()
    img = Image.new("RGBA", (W, H), (*SOFT_BG, 255))
    if frame == "A":
        draw_progressive_header(img, PROGRESSIVE_SUBTITLE_A, "detector")
        draw_progressive_detector_map(img, detector_img)
        draw_progressive_detector_summary(img)
        draw_progressive_bridge(img, PROGRESSIVE_BRIDGE_A)
        png = output_dir / "hp2026_slide02A_detector_subsystems_context.png"
    elif frame == "B":
        draw_progressive_header(img, PROGRESSIVE_SUBTITLE_B, "data")
        draw_progressive_detector_map(img, detector_img)
        draw_progressive_data_panel(img)
        draw_progressive_bridge(img, PROGRESSIVE_BRIDGE_B)
        png = output_dir / "hp2026_slide02B_detector_data_context_build.png"
    else:
        raise ValueError(f"unknown progressive frame: {frame}")
    img.convert("RGB").save(png, "PNG")
    return png


def write_progressive_script(output_dir: Path) -> Path:
    script = """# HP2026 Slide 2 Progressive Build Speaker Script

Frame 2A, before the click:

I want to start with the experimental context, not the photon physics yet. sPHENIX is the new RHIC detector built around tracking and full calorimetry. For this measurement the important handles are the EMCal for energy and shower shape, the hadronic calorimeters for nearby activity and isolation, the tracking system for event and vertex context, and the MBD/trigger/luminosity systems that make a cross section possible.

Click to Frame 2B:

Now we add the data context. The broad p+p running available to sPHENIX includes the 2024 calorimeter-only sample and the all-subsystem p+p running from 2024 and 2026. I would keep the PPG12 analysis luminosity separate: the isolated prompt-photon result uses the defined Run 24 analysis sample, L = 64.4 pb^-1. The Au+Au and O+O numbers are shown only as muted context for the future hard-probes program, not as the focus of this talk.

Transition to Slide 3:

With the detector and data context established, the next slide can narrow cleanly to the physics object: why an isolated prompt photon is a useful p+p hard-scattering baseline.
"""
    path = output_dir / "hp2026_slide02_progressive_build_script.md"
    path.write_text(script, encoding="utf-8")
    return path


def write_progressive_contact_sheet(output_dir: Path, frame_a: Path, frame_b: Path) -> Path:
    contact = output_dir / "hp2026_slide02_progressive_build_contact_sheet.png"
    thumb_w, thumb_h = 960, 540
    label_h = 56
    sheet = Image.new("RGB", (2 * thumb_w, thumb_h + label_h), "white")
    draw = ImageDraw.Draw(sheet)
    for idx, path in enumerate((frame_a, frame_b)):
        x = idx * thumb_w
        draw.text((x + 18, 14), path.stem, font=font(TIMES, 27), fill=INK)
        slide = Image.open(path).convert("RGB").resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        sheet.paste(slide, (x, label_h))
    sheet.save(contact, "PNG")
    return contact


def write_progressive_manifest(output_dir: Path, frame_a: Path, frame_b: Path, script: Path, contact: Path) -> Path:
    _detector_img, detector_asset = load_detector_rendering_asset()
    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "google_slides_mutation": False,
        "purpose": "Two-frame click-through replacement candidate for HPslides_v1 Slide 2.",
        "outputs": {
            "frame_2a": str(frame_a.relative_to(ROOT)),
            "frame_2b": str(frame_b.relative_to(ROOT)),
            "speaker_script": str(script.relative_to(ROOT)),
            "contact_sheet": str(contact.relative_to(ROOT)),
        },
        "size": [W, H],
        "mode": "RGB",
        "source_assets": {
            "detector_asset": detector_asset,
            "official_source_urls": REAL_DETECTOR_SOURCE,
        },
        "data_source_evidence": PROGRESSIVE_DATA_EVIDENCE,
        "design_notes": [
            "Frame 2A and Frame 2B use the same detector anchor to behave like one spoken slide with one click-through build.",
            "The previous photon squiggle overlay is removed; detector annotation uses clean leader lines and callout boxes.",
            "Subsystem callout leaders are presentation-grade visual anchors, not CAD-precision subsystem boundaries.",
            "p+p data are visually primary; Au+Au and O+O are muted future-program context.",
            "No slide number or provenance footer is baked into either PNG.",
        ],
    }
    manifest_path = output_dir / "hp2026_slide02_progressive_build_manifest.json"
    with manifest_path.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return manifest_path


def render_detector_data_progressive_build(output_dir: Path) -> Path:
    frame_a = render_detector_data_progressive_frame(output_dir, "A")
    frame_b = render_detector_data_progressive_frame(output_dir, "B")
    script = write_progressive_script(output_dir)
    contact = write_progressive_contact_sheet(output_dir, frame_a, frame_b)
    write_progressive_manifest(output_dir, frame_a, frame_b, script, contact)
    return contact


def draw_standalone_header(base: Image.Image, title: str, subtitle: str) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*PHOTON, 255))
    logo = load_sphenix_logo()
    if logo is not None:
        paste_fit(base, logo, (2078, 78, 2388, 160), anchor="right")
    draw.text((132, 86), title, font=font(TIMES_BOLD, 76), fill=INK)
    draw_wrapped(draw, subtitle, (136, 184), 1660, font(TIMES_ITALIC, 35), fill=BLUE, line_gap=5)
    draw.line((132, 300, W - 132, 300), fill=(221, 226, 232), width=3)


def draw_standalone_bridge(base: Image.Image, text: str) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    draw.rounded_rectangle((132, 1228, W - 132, 1318), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw.text((174, 1253), text, font=font(TIMES_ITALIC, 38), fill=BLUE)


def draw_standalone_experiment_summary(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    panel = (1512, 340, 2428, 1088)
    draw.rounded_rectangle(panel, radius=12, fill=(255, 255, 255, 255), outline=(216, 225, 234, 255), width=2)
    draw.text((1550, 372), "What sPHENIX gives this analysis", font=font(TIMES_BOLD, 38), fill=BLUE)
    draw.line((1550, 430, 2388, 430), fill=(222, 228, 235), width=2)

    y = 472
    for heading, body, color in STANDALONE_EXPERIMENT_TAKEAWAYS:
        draw.rounded_rectangle((1550, y, 2388, y + 148), radius=9, fill=(247, 250, 252, 255), outline=(222, 229, 236, 255), width=2)
        draw.rounded_rectangle((1550, y, 1564, y + 148), radius=5, fill=(*color, 255))
        draw.text((1594, y + 24), heading, font=font(TIMES_BOLD, 32), fill=INK)
        draw_wrapped(draw, body, (1594, y + 76), 740, font(TIMES, 26), fill=MUTED, line_gap=7)
        y += 174


def draw_experiment_intro_card(
    draw: ImageDraw.ImageDraw,
    box_xy: tuple[int, int],
    label: str,
    expansion: str,
    body: str,
    color: tuple[int, int, int],
) -> tuple[int, int, int, int]:
    x, y = box_xy
    w = 404
    h = EXPERIMENT_INTRO_CARD_HEIGHTS.get(label, 174)
    if x < W / 2:
        w = 430
    draw.rounded_rectangle((x, y, x + w, y + h), radius=10, fill=(255, 255, 255, 248), outline=(194, 208, 222, 255), width=3)
    draw.rounded_rectangle((x, y, x + 18, y + h), radius=6, fill=(*color, 255))
    draw.ellipse((x + 32, y + 28, x + 66, y + 62), fill=(*color, 215))
    draw.text((x + 82, y + 16), label, font=font(TIMES_BOLD, 34), fill=INK)
    draw_wrapped(draw, expansion, (x + 82, y + 58), w - 112, font(TIMES_ITALIC, 24), fill=BLUE, line_gap=3)
    draw_wrapped(draw, body, (x + 30, y + 102), w - 60, font(TIMES, 21), fill=MUTED, line_gap=4)
    return (x, y, x + w, y + h)


def draw_subsystem_group_icon(draw: ImageDraw.ImageDraw, kind: str, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    cx = (x0 + x1) // 2
    cy = (y0 + y1) // 2
    if kind == "tracking":
        for offset, width in ((-22, 3), (0, 3), (22, 3)):
            draw.arc((x0 + 4, y0 + 6 + offset, x1 - 8, y1 - 8 + offset), start=205, end=330, fill=(70, 78, 88, 230), width=width)
            draw.arc((x0 + 18, y0 + 12 + offset, x1 - 22, y1 - 14 + offset), start=205, end=330, fill=(126, 139, 151, 220), width=width)
        for dx in (-30, -10, 14, 34):
            draw.line((cx + dx, y0 + 10, cx + dx - 10, y1 - 8), fill=(70, 78, 88, 210), width=3)
    elif kind == "calorimetry":
        colors = [(220, 226, 214), (230, 211, 177), (179, 212, 225)]
        for idx, color in enumerate(colors):
            y = y0 + 18 + idx * 22
            poly = [(x0 + 18, y), (x0 + 64, y + 16), (x0 + 42, y + 34), (x0 - 4, y + 16)]
            draw.polygon(poly, fill=(*color, 255), outline=(70, 78, 88, 230))
            draw.line((x0 + 64, y + 16, x0 + 64, y + 35, x0 + 42, y + 52), fill=(70, 78, 88, 180), width=2)
    else:
        draw.line((x0 + 6, cy, x1 - 16, cy - 34), fill=(70, 78, 88, 215), width=3)
        draw.line((x0 + 8, cy, x1 - 14, cy + 24), fill=(70, 78, 88, 215), width=3)
        draw.line((x0 + 8, cy, x1 - 44, y0 + 8), fill=(70, 78, 88, 215), width=3)
        draw.rectangle((x1 - 34, cy - 48, x1 - 8, cy - 24), fill=(207, 221, 228, 255), outline=(70, 78, 88, 230), width=2)
        draw.rectangle((x1 - 28, cy + 12, x1 - 4, cy + 38), fill=(224, 205, 178, 255), outline=(70, 78, 88, 230), width=2)
        draw.polygon([(x0 + 8, cy), (x0 + 24, cy - 9), (x0 + 24, cy + 9)], fill=(*PHOTON_DARK, 210))


def draw_subsystem_group_text(draw: ImageDraw.ImageDraw) -> None:
    panel = (132, 344, 920, 1190)
    draw.rounded_rectangle(panel, radius=12, fill=(255, 255, 255, 248), outline=(218, 226, 235, 255), width=2)
    draw.text((174, 384), "Subsystems used in this measurement", font=font(TIMES_BOLD, 37), fill=BLUE)
    draw.line((174, 440, 878, 440), fill=(220, 228, 236, 255), width=2)

    y_positions = [492, 704, 916]
    for (kind, heading, parenthetical, body), y in zip(SUBSYSTEM_GROUPS, y_positions):
        draw_subsystem_group_icon(draw, kind, (174, y + 8, 264, y + 104))
        heading_font = font(TIMES_BOLD, 34)
        heading_w, _ = text_box(draw, heading, heading_font)
        draw.text((300, y + 4), heading, font=heading_font, fill=INK)
        draw.text((300 + heading_w + 8, y + 8), f"({parenthetical})", font=font(TIMES, 27), fill=INK)
        draw_wrapped(draw, body, (300, y + 52), 560, font(TIMES, 28), fill=INK, line_gap=6)


def draw_detector_direct_arrow(
    draw: ImageDraw.ImageDraw,
    label_pos: tuple[int, int],
    target: tuple[int, int],
    label: str,
    align: str,
) -> None:
    lx, ly = label_pos
    label_font = font(TIMES_BOLD, 35)
    tw, th = text_box(draw, label, label_font)
    if align == "center":
        text_x = lx - tw // 2
    elif align == "right":
        text_x = lx - tw
    else:
        text_x = lx
    text_y = ly
    draw.text((text_x, text_y), label, font=label_font, fill=INK, stroke_width=3, stroke_fill=(255, 255, 255, 235))

    start = (lx, ly + th + 12)
    if ly > target[1]:
        start = (lx, ly - 8)
    dx = target[0] - start[0]
    dy = target[1] - start[1]
    if not (dx or dy):
        return
    angle = math.atan2(dy, dx)
    line_end = (
        target[0] - math.cos(angle) * 18,
        target[1] - math.sin(angle) * 18,
    )
    draw.line((start, line_end), fill=(255, 255, 255, 235), width=13)
    draw.line((start, line_end), fill=(*PHOTON_DARK, 230), width=7)
    size = 28
    spread = 0.55
    base = (
        target[0] - math.cos(angle) * size,
        target[1] - math.sin(angle) * size,
    )
    p1 = (
        base[0] + math.cos(angle + math.pi / 2) * size * spread,
        base[1] + math.sin(angle + math.pi / 2) * size * spread,
    )
    p2 = (
        base[0] + math.cos(angle - math.pi / 2) * size * spread,
        base[1] + math.sin(angle - math.pi / 2) * size * spread,
    )
    draw.polygon([target, p1, p2], fill=(255, 255, 255, 238))
    draw.polygon([target, p1, p2], fill=(*PHOTON_DARK, 235))


def draw_experiment_intro_layout(base: Image.Image, detector_img: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    draw_subsystem_group_text(draw)

    panel = (950, 326, 2428, 1220)
    image_box = (980, 342, 2400, 1190)

    shadow = Image.new("RGBA", base.size, (0, 0, 0, 0))
    sdraw = ImageDraw.Draw(shadow, "RGBA")
    sdraw.rounded_rectangle((panel[0] + 6, panel[1] + 8, panel[2] + 6, panel[3] + 8), radius=14, fill=(30, 42, 58, 32))
    shadow = shadow.filter(ImageFilter.GaussianBlur(12))
    base.alpha_composite(shadow)

    draw.rounded_rectangle(panel, radius=14, fill=(255, 255, 255, 255), outline=(216, 225, 234, 255), width=2)
    draw.rounded_rectangle(image_box, radius=10, fill=(248, 250, 252, 255), outline=(224, 231, 238, 255), width=1)

    rendering = detector_img.convert("RGBA")
    rendering = ImageEnhance.Color(rendering).enhance(1.08)
    rendering = ImageEnhance.Contrast(rendering).enhance(1.06)
    rendering = ImageEnhance.Sharpness(rendering).enhance(1.12)
    image_bounds = paste_fit_return_box(base, rendering, image_box)

    overlay = Image.new("RGBA", base.size, (0, 0, 0, 0))
    odraw = ImageDraw.Draw(overlay, "RGBA")
    ix0, iy0, ix1, iy1 = image_bounds
    for label, label_pos, target_rel, align in SUBSYSTEM_LABEL_ARROWS:
        target = (
            round(ix0 + target_rel[0] * (ix1 - ix0)),
            round(iy0 + target_rel[1] * (iy1 - iy0)),
        )
        draw_detector_direct_arrow(odraw, label_pos, target, label, align)
    base.alpha_composite(overlay)


def render_standalone_experiment_slide(output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    detector_img, _detector_asset = load_detector_rendering_asset()
    img = Image.new("RGBA", (W, H), (*SOFT_BG, 255))
    draw_standalone_header(img, STANDALONE_EXPERIMENT_TITLE, STANDALONE_EXPERIMENT_SUBTITLE)
    draw_experiment_intro_layout(img, detector_img)
    png = output_dir / "hp2026_slide02_sphenix_experiment_for_photons.png"
    img.convert("RGB").save(png, "PNG")
    write_standalone_experiment_script(output_dir)
    return png


def draw_dataset_primary_panel(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    panel = (132, 340, 1560, 926)
    draw.rounded_rectangle(panel, radius=12, fill=(255, 255, 255, 255), outline=(216, 225, 234, 255), width=2)
    draw.text((178, 382), "Available p+p data set the context", font=font(TIMES_BOLD, 43), fill=BLUE)
    draw_wrapped(
        draw,
        "These are the broad sPHENIX p+p running numbers to orient the audience before the analysis details.",
        (180, 438),
        1285,
        font(TIMES_ITALIC, 28),
        fill=MUTED,
        line_gap=5,
    )

    card_y = 548
    card_w = 394
    gap = 54
    for idx, (run, value, unit, detail, color) in enumerate(STANDALONE_PP_AVAILABILITY):
        x = 178 + idx * (card_w + gap)
        draw.rounded_rectangle((x, card_y, x + card_w, card_y + 260), radius=12, fill=(247, 250, 252, 255), outline=(222, 229, 236, 255), width=2)
        draw.rounded_rectangle((x, card_y, x + card_w, card_y + 18), radius=7, fill=(*color, 255))
        draw.text((x + 28, card_y + 38), run, font=font(TIMES_BOLD, 31), fill=INK)
        value_font = font(TIMES_BOLD, 76)
        vw, _ = text_box(draw, value, value_font)
        draw.text((x + (card_w - vw) / 2, card_y + 92), value, font=value_font, fill=INK)
        unit_font = font(TIMES, 34)
        uw, _ = rich_text_box(draw, unit, unit_font)
        draw_rich_text(draw, (round(x + (card_w - uw) / 2), card_y + 168), unit, unit_font, BLUE)
        detail_font = font(TIMES_ITALIC, 27)
        dw, _ = text_box(draw, detail, detail_font)
        draw.text((x + (card_w - dw) / 2, card_y + 214), detail, font=detail_font, fill=MUTED)

    draw.rounded_rectangle((178, 840, 1514, 884), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=1)
    draw.text((206, 852), "p+p is the baseline system; the analysis sample is defined separately below.", font=font(TIMES_ITALIC, 24), fill=BLUE)


def draw_dataset_analysis_sample(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    panel = (132, 972, 1560, 1166)
    draw.rounded_rectangle(panel, radius=12, fill=(239, 246, 250, 255), outline=(205, 221, 232, 255), width=2)
    draw.rounded_rectangle((132, 972, 150, 1166), radius=6, fill=(*TEAL, 255))
    draw.text((184, 1006), "PPG12 result uses a defined Run 24 analysis sample", font=font(TIMES_BOLD, 38), fill=INK)
    draw_rich_text(draw, (184, 1062), "L = 64.4 pb^-1", font(TIMES_BOLD, 62), BLUE)
    draw_wrapped(
        draw,
        "Luminosity for the isolated prompt-photon cross section; separate from the total p+p running shown above.",
        (720, 1048),
        720,
        font(TIMES, 27),
        fill=MUTED,
        line_gap=7,
    )


def draw_dataset_side_context(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    panel = (1620, 340, 2428, 1166)
    draw.rounded_rectangle(panel, radius=12, fill=(255, 255, 255, 255), outline=(216, 225, 234, 255), width=2)
    draw.text((1662, 382), "Broader program context", font=font(TIMES_BOLD, 39), fill=BLUE)
    draw_wrapped(
        draw,
        "Heavy-ion running motivates why the p+p photon result matters, but it should stay secondary on this opening data slide.",
        (1664, 438),
        700,
        font(TIMES_ITALIC, 27),
        fill=MUTED,
        line_gap=6,
    )

    y = 548
    for idx, (system, value, detail) in enumerate(STANDALONE_HEAVY_ION_CONTEXT):
        color = LIGHT_MUTED if idx else TEAL_SOFT
        draw.rounded_rectangle((1662, y, 2388, y + 112), radius=9, fill=(249, 250, 251, 255), outline=(224, 230, 236, 255), width=1)
        draw.rounded_rectangle((1662, y, 1676, y + 112), radius=5, fill=(*color, 180))
        draw.text((1704, y + 24), system, font=font(TIMES_BOLD, 29), fill=MUTED)
        draw_rich_text(draw, (1930, y + 24), value, font(TIMES_BOLD, 31), MUTED)
        draw.text((1704, y + 68), detail, font=font(TIMES, 24), fill=LIGHT_MUTED)
        y += 138

    draw.rounded_rectangle((1662, 970, 2388, 1118), radius=9, fill=(247, 250, 252, 255), outline=(222, 229, 236, 255), width=2)
    draw.text((1704, 994), "Why p+p matters", font=font(TIMES_BOLD, 30), fill=INK)
    draw_wrapped(
        draw,
        "p+p provides the calibrated photon baseline for later hard-probes measurements.",
        (1704, 1036),
        650,
        font(TIMES, 24),
        fill=MUTED,
        line_gap=5,
    )


def render_standalone_dataset_slide(output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    img = Image.new("RGBA", (W, H), (*SOFT_BG, 255))
    draw_standalone_header(img, STANDALONE_DATASET_TITLE, STANDALONE_DATASET_SUBTITLE)
    draw_dataset_primary_panel(img)
    draw_dataset_analysis_sample(img)
    draw_dataset_side_context(img)
    draw_standalone_bridge(img, STANDALONE_DATASET_BRIDGE)
    png = output_dir / "hp2026_slide03_pp_dataset_context.png"
    img.convert("RGB").save(png, "PNG")
    write_standalone_dataset_script(output_dir)
    return png


def write_standalone_experiment_script(output_dir: Path) -> Path:
    script = """# HP2026 Slide 2 Script - sPHENIX Experiment

At a high level, I want to start with the detector itself, because this measurement is built from the experimental handles that sPHENIX gives us.

sPHENIX is a large-acceptance detector at RHIC, designed around tracking, electromagnetic and hadronic calorimetry, a solenoidal magnetic field, and forward trigger and luminosity systems. For this talk, the important point is not every detail of the apparatus, but how these subsystems come together to make an isolated prompt-photon cross-section measurement possible.

Reading from the inside outward, the tracking system, MVTX, INTT, TPC, and TPOT, reconstructs charged tracks and the collision vertex. That gives the event context and helps distinguish charged activity from neutral electromagnetic energy.

Surrounding that is the calorimetry. The electromagnetic calorimeter, or EMCal, is the primary detector for the photon candidate: it measures the cluster energy and provides shower-shape information used later for photon identification. The inner and outer hadronic calorimeters, iHCal and oHCal, sample nearby hadronic activity, which matters because isolation is fundamentally asking whether the photon candidate is surrounded by additional event activity.

The detector sits inside a 1.4 T solenoid, so charged particles bend and can be reconstructed with momentum information. And the MBD, the Minimum Bias Detector, provides fast triggering, timing, vertex context, and the luminosity normalization needed to turn selected candidates into a cross section.

One subtle point is that the MBD sits forward along the beam direction, while this rendering is mainly useful for orienting the central tracking and calorimeter barrel. So overall, this slide is the experimental foundation: tracking gives the event and charged-particle context, calorimetry gives the photon and isolation handles, and MBD/luminosity make the measurement quantitative. With that detector context established, the next step is to separate the broad p+p running available to sPHENIX from the specific Run 24 sample used for the PPG12 result.
"""
    path = output_dir / "hp2026_slide02_sphenix_experiment_for_photons_script.md"
    path.write_text(script, encoding="utf-8")
    return path


def write_standalone_dataset_script(output_dir: Path) -> Path:
    script = """# HP2026 Slide 3 Script - p+p Dataset Context

Now that the detector context is set, I want to separate two pieces of the data story that can easily get mixed together.

First is the broad p+p data available to sPHENIX. In 2024, sPHENIX recorded 107 inverse picobarns with calorimeter information and 13 inverse picobarns with all subsystems. In 2026, there is an additional 17 inverse picobarns with all subsystems. Those numbers tell us the scale of the p+p program and why this is becoming a precision baseline environment.

Second is the actual analysis sample used for the isolated prompt-photon result. The current PPG12 measurement uses the defined Run 24 sample with an integrated luminosity of 64.4 inverse picobarns. So when I quote 64.4 inverse picobarns, that is not the total available p+p running; it is the luminosity for the cross-section result being shown in this talk.

I also want to show the heavy-ion context, but keep it secondary. The 2025 Au+Au and 2026 O+O data are what make this p+p measurement strategically important: p+p establishes the baseline that future hard-probes measurements will need.

So the takeaway from this slide is simple: sPHENIX has the detector and data context for a p+p isolated prompt-photon baseline, and the analysis uses a clearly defined Run 24 sample. With that in place, the next slide can define the physics object itself: the isolated prompt photon.
"""
    path = output_dir / "hp2026_slide03_pp_dataset_context_script.md"
    path.write_text(script, encoding="utf-8")
    return path


def write_standalone_contact_sheet(output_dir: Path, slide2: Path, slide3: Path) -> Path:
    contact = output_dir / "hp2026_experiment_dataset_standalone_contact_sheet.png"
    thumb_w, thumb_h = 960, 540
    label_h = 56
    sheet = Image.new("RGB", (2 * thumb_w, thumb_h + label_h), "white")
    draw = ImageDraw.Draw(sheet)
    for idx, path in enumerate((slide2, slide3)):
        x = idx * thumb_w
        draw.text((x + 18, 14), path.stem, font=font(TIMES, 27), fill=INK)
        slide = Image.open(path).convert("RGB").resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        sheet.paste(slide, (x, label_h))
    sheet.save(contact, "PNG")
    return contact


def write_standalone_experiment_dataset_manifest(output_dir: Path, slide2: Path, slide3: Path, contact: Path) -> Path:
    _detector_img, detector_asset = load_detector_rendering_asset()
    script2 = output_dir / "hp2026_slide02_sphenix_experiment_for_photons_script.md"
    script3 = output_dir / "hp2026_slide03_pp_dataset_context_script.md"
    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "google_slides_mutation": False,
        "recommendation": "Use these as two normal consecutive slides, not as a progressive click-through build.",
        "outputs": {
            "slide_2_experiment": str(slide2.relative_to(ROOT)),
            "slide_3_dataset": str(slide3.relative_to(ROOT)),
            "contact_sheet": str(contact.relative_to(ROOT)),
            "slide_2_script": str(script2.relative_to(ROOT)),
            "slide_3_script": str(script3.relative_to(ROOT)),
        },
        "size": [W, H],
        "mode": "RGB",
        "source_assets": {
            "detector_asset": detector_asset,
            "official_source_urls": REAL_DETECTOR_SOURCE,
        },
        "detector_subsystem_reference": {
            "source_file": "usefulDocs/20260506_DIS_YeonjuGo.pdf",
            "reference_slide": "Yeonju DIS2026 sPHENIX Experiment at RHIC detector overview",
            "subsystems_checked": [
                "MVTX / INTT / TPC / TPOT tracking",
                "Electromagnetic Calorimeter (EMCal)",
                "Inner and outer Hadronic Calorimeters (iHCal / oHCal)",
                "Superconducting magnet",
                "Minimum Bias Detector (MBD)",
                "Barrel/midrapidity tracking-plus-calorimetry coverage",
            ],
            "annotation_policy": "Leader lines anchor visible presentation-grade regions in the public detector rendering. MBD is marked as a forward beam-axis vicinity, and the barrel coverage card deliberately has no leader because it is not a separate detector component.",
        },
        "data_source_evidence": STANDALONE_DATASET_EVIDENCE,
        "design_notes": [
            "Standalone slides are preferred here because they are easier to present, review, screenshot, and insert than a hidden progressive build.",
            "Slide 2 introduces the experiment only: the detector is centered, subsystem acronyms are expanded, and surrounding cards explain tracking, calorimetry, magnet, trigger, and luminosity context.",
            "Slide 2 callout leaders were corrected to avoid false detector precision: EMCal, HCal, tracking, and magnet point to visible subsystem regions; MBD points only to the forward beam-axis vicinity; barrel coverage has no arrow.",
            "Slide 2 callout leaders use white-underlay strokes, stronger colored lines, arrowheads, and large endpoint rings so the subsystem mapping is readable at presentation distance.",
            "Slide 2 gives the EMCal and Barrel coverage cards extra vertical room so their final body-text lines are not cramped.",
            "Slide 3 introduces the dataset only: broad p+p availability, separate PPG12 Run 24 analysis luminosity, and muted heavy-ion context.",
            "Slide 4 can remain the isolated prompt-photon motivation slide, preserving a clean progression: experiment -> data -> object.",
            "No slide number or provenance footer is baked into either PNG.",
        ],
    }
    manifest_path = output_dir / "hp2026_experiment_dataset_standalone_manifest.json"
    with manifest_path.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return manifest_path


def render_experiment_dataset_standalone(output_dir: Path) -> Path:
    slide2 = render_standalone_experiment_slide(output_dir)
    slide3 = render_standalone_dataset_slide(output_dir)
    contact = write_standalone_contact_sheet(output_dir, slide2, slide3)
    write_standalone_experiment_dataset_manifest(output_dir, slide2, slide3, contact)
    return contact


def write_detector_data_script(output_dir: Path) -> Path:
    script = """# HP2026 Slide 2 Speaker Script

I want to start with just enough sPHENIX context to make the photon measurement legible. sPHENIX is the new RHIC detector built for high-rate hard-probe measurements, with tracking and full calorimetry that are directly relevant for photons.

For this talk the important detector handles are the EMCal shower and energy measurement, the nearby calorimeter energy used for isolation, tracking information that helps characterize the event, and the trigger/luminosity context from the minimum-bias systems.

The data context is p+p. The reference slides consistently quote 2024 p+p as 107/pb calorimeter-only and 13/pb with all subsystems, and 2026 p+p as 17/pb with all subsystems. I am not trying to summarize the entire sPHENIX run program here; I am setting the baseline system for this prompt-photon cross-section measurement.

Transition: with the detector and p+p sample established, the next slide defines why the isolated prompt photon is such a clean object.
"""
    path = output_dir / "hp2026_slide02_detector_data_context_script.md"
    path.write_text(script, encoding="utf-8")
    return path


def render_detector_data_context(output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    detector_img, _detector_asset = load_detector_rendering_asset()
    img = Image.new("RGBA", (W, H), (*SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")

    draw.rectangle((0, 0, W, H), fill=(*SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*PHOTON, 255))

    logo = load_sphenix_logo()
    if logo is not None:
        paste_fit(img, logo, (2090, 78, 2388, 160), anchor="right")

    draw.text((132, 86), DETECTOR_DATA_TITLE, font=font(TIMES_BOLD, 72), fill=INK)
    draw.text((136, 184), DETECTOR_DATA_SUBTITLE, font=font(TIMES_ITALIC, 38), fill=BLUE)
    draw.line((132, 286, W - 132, 286), fill=(221, 226, 232), width=3)

    draw_detector_data_image_panel(img, detector_img)
    draw_detector_data_cards(img)
    draw_pp_data_panel(img)

    draw.rounded_rectangle((132, 1228, W - 132, 1318), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw.text((174, 1253), DETECTOR_DATA_BRIDGE, font=font(TIMES_ITALIC, 38), fill=BLUE)

    png = output_dir / "hp2026_slide02_detector_data_context.png"
    img.convert("RGB").save(png, "PNG")
    write_detector_data_script(output_dir)
    return png


def draw_isolated_photon_visual(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    panel = (126, 333, 1422, 1112)
    draw.rounded_rectangle(panel, radius=10, fill=(*PANEL, 255), outline=(*PANEL_EDGE, 255), width=2)

    draw.text((184, 382), "p+p hard scattering", font=font(TIMES_ITALIC, 33), fill=MUTED)
    y0 = 666
    for x, label in ((300, "p"), (514, "p")):
        draw.ellipse((x - 48, y0 - 48, x + 48, y0 + 48), fill=(255, 255, 255, 255), outline=(98, 139, 169, 255), width=4)
        draw.text((x - 13, y0 - 23), label, font=font(TIMES_BOLD, 45), fill=BLUE)
    draw.line((348, y0, 628, y0), fill=(178, 199, 216, 255), width=6)
    draw.line((466, y0, 628, y0), fill=(178, 199, 216, 255), width=6)

    collision = (654, y0)
    for r, alpha in ((78, 28), (52, 46), (25, 96)):
        draw.ellipse((collision[0] - r, collision[1] - r, collision[0] + r, collision[1] + r), fill=(*SPHENIX_BLUE, alpha))
    draw.ellipse((collision[0] - 12, collision[1] - 12, collision[0] + 12, collision[1] + 12), fill=(*SPHENIX_BLUE, 220))

    # Quiet colored activity is drawn outside the isolation cone to explain
    # suppression without turning this into a photon+jet opening.
    for angle, length, width, alpha in ((112, 170, 5, 75), (248, 155, 5, 70), (308, 132, 4, 58)):
        ex = collision[0] + math.cos(math.radians(angle)) * length
        ey = collision[1] + math.sin(math.radians(angle)) * length
        draw.line((collision[0], collision[1], ex, ey), fill=(*TEAL_SOFT, alpha), width=width)
        draw.ellipse((ex - 9, ey - 9, ex + 9, ey + 9), fill=(*TEAL_SOFT, alpha + 25))

    photon_start = (686, 636)
    cone_tip = (1115, 546)
    detector = (1238, 454, 1324, 692)
    dx, dy = cone_tip[0] - photon_start[0], cone_tip[1] - photon_start[1]
    length = math.hypot(dx, dy)
    nx, ny = -dy / length, dx / length
    cone_left = (cone_tip[0] + nx * 124, cone_tip[1] + ny * 124)
    cone_right = (cone_tip[0] - nx * 124, cone_tip[1] - ny * 124)
    alpha_layer(base, lambda d: d.polygon([photon_start, cone_left, cone_right], fill=(245, 181, 34, 28)))
    draw.line((photon_start, cone_left), fill=(219, 157, 32, 80), width=3)
    draw.line((photon_start, cone_right), fill=(219, 157, 32, 80), width=3)
    draw.arc((cone_tip[0] - 130, cone_tip[1] - 130, cone_tip[0] + 130, cone_tip[1] + 130), start=346, end=44, fill=(*PHOTON_DARK, 160), width=4)

    wave = feynman_points(photon_start, cone_tip, amplitude=16, cycles=7.3, steps=240)
    glow = Image.new("RGBA", base.size, (0, 0, 0, 0))
    gdraw = ImageDraw.Draw(glow, "RGBA")
    draw_polyline(gdraw, wave, (*PHOTON, 118), 17)
    glow = glow.filter(ImageFilter.GaussianBlur(7))
    base.alpha_composite(glow)
    draw = ImageDraw.Draw(base, "RGBA")
    draw_polyline(draw, wave, (*PHOTON, 255), 9)
    draw_polyline(draw, wave, (*PHOTON_DARK, 212), 3)

    draw.rounded_rectangle(detector, radius=8, fill=(236, 244, 250, 255), outline=(172, 194, 212, 255), width=3)
    for i in range(1, 5):
        x = detector[0] + i * (detector[2] - detector[0]) / 5
        draw.line((x, detector[1], x, detector[3]), fill=(197, 211, 223, 255), width=2)
    for i in range(1, 7):
        y = detector[1] + i * (detector[3] - detector[1]) / 7
        draw.line((detector[0], y, detector[2], y), fill=(197, 211, 223, 255), width=2)
    draw.ellipse((cone_tip[0] - 21, cone_tip[1] - 21, cone_tip[0] + 21, cone_tip[1] + 21), fill=(*PHOTON, 235), outline=(*PHOTON_DARK, 210), width=3)

    draw.text((810, 472), "isolated prompt photon", font=font(TIMES_ITALIC, 35), fill=BLUE)
    draw.text((954, 690), "quiet cone", font=font(TIMES_ITALIC, 29), fill=(145, 111, 43))

    bg_box = (236, 870, 1230, 1018)
    draw.rounded_rectangle(bg_box, radius=8, fill=(255, 255, 255, 235), outline=(222, 229, 236, 255), width=2)
    y = draw_wrapped(
        draw,
        "Isolation turns the object from a crowded EMCal cluster into a cleaner hard-scattering tag.",
        (272, 900),
        900,
        font(TIMES, 29),
        fill=INK,
        line_gap=6,
    )
    draw_wrapped(
        draw,
        "The p+p cross section is the baseline result; gamma-jet and heavy-ion comparisons come later.",
        (272, y + 8),
        900,
        font(TIMES_ITALIC, 26),
        fill=MUTED,
        line_gap=5,
    )


def draw_photon_motivation_cards(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x, y = 1510, 345
    w, h, gap = 870, 216, 34
    accents = [PHOTON, SPHENIX_BLUE, TEAL]
    for idx, (heading, body) in enumerate(PHOTON_MOTIVATION_REASONS):
        top = y + idx * (h + gap)
        draw.rounded_rectangle((x, top, x + w, top + h), radius=10, fill=(*CARD, 255), outline=(*CARD_EDGE, 255), width=2)
        draw.rounded_rectangle((x, top, x + 14, top + h), radius=6, fill=(*accents[idx], 255))
        draw_icon(draw, idx, (x + 86, top + 108))
        draw.text((x + 162, top + 41), heading, font=font(TIMES_BOLD, 43), fill=INK)
        draw_wrapped(draw, body, (x + 164, top + 101), 625, font(TIMES, 33), fill=MUTED, line_gap=10)


def write_isolated_photon_script(output_dir: Path) -> Path:
    script = """# HP2026 Slide 3 Speaker Script

Now I can define the physics object. A prompt photon is tied to the short-distance hard scattering, and because it is color neutral it does not lose energy through the strong final-state interactions that affect colored probes.

Experimentally, the word isolated matters. We require the energy around the photon candidate to be quiet, which strongly reduces decay photons and fragmentation-rich backgrounds before the purity correction. That is what makes the candidate sample interpretable as a prompt-photon measurement.

The result I want the audience to hold onto is the p+p cross section. It is valuable on its own, and it also becomes the reference point for future gamma-jet and heavy-ion measurements at sPHENIX.

Transition: after this motivation, the talk can move into the actual reconstruction chain: clusters, identification, isolation, purity, correction, and then the cross-section comparison.
"""
    path = output_dir / "hp2026_slide03_isolated_photon_motivation_script.md"
    path.write_text(script, encoding="utf-8")
    return path


def render_isolated_photon_motivation(output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    img = Image.new("RGBA", (W, H), (*SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")

    draw.rectangle((0, 0, W, H), fill=(*SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*PHOTON, 255))

    draw.text((132, 98), PHOTON_MOTIVATION_TITLE, font=font(TIMES_BOLD, 90), fill=INK)
    draw.text((136, 205), PHOTON_MOTIVATION_SUBTITLE, font=font(TIMES_ITALIC, 43), fill=BLUE)
    draw.line((132, 292, W - 132, 292), fill=(221, 226, 232), width=3)

    draw_isolated_photon_visual(img)
    draw_photon_motivation_cards(img)

    draw.rounded_rectangle((132, 1192, W - 132, 1282), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw.text((174, 1216), PHOTON_MOTIVATION_BRIDGE, font=font(TIMES_ITALIC, 38), fill=BLUE)

    png = output_dir / "hp2026_slide03_isolated_photon_motivation.png"
    img.convert("RGB").save(png, "PNG")
    write_isolated_photon_script(output_dir)
    return png


def write_opening_sequence_manifest(output_dir: Path, slide2: Path, slide3: Path) -> Path:
    detector_script = output_dir / "hp2026_slide02_detector_data_context_script.md"
    photon_script = output_dir / "hp2026_slide03_isolated_photon_motivation_script.md"
    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "primary_outputs": [
            str(slide2.relative_to(ROOT)),
            str(slide3.relative_to(ROOT)),
        ],
        "speaker_scripts": [
            str(detector_script.relative_to(ROOT)),
            str(photon_script.relative_to(ROOT)),
        ],
        "size": [W, H],
        "mode": "RGB",
        "google_slides_mutation": False,
        "sequence": [
            {
                "slide": 2,
                "output": str(slide2.relative_to(ROOT)),
                "title": DETECTOR_DATA_TITLE,
                "role": "first content slide after title; establish sPHENIX detector and p+p data context",
                "detector_asset": str(REAL_DETECTOR_RENDERING.relative_to(ROOT)),
                "official_source_urls": REAL_DETECTOR_SOURCE,
                "data_taking_numbers": DATA_SOURCE_EVIDENCE["values_on_slide"],
                "speaker_script": str(detector_script.relative_to(ROOT)),
            },
            {
                "slide": 3,
                "output": str(slide3.relative_to(ROOT)),
                "title": PHOTON_MOTIVATION_TITLE,
                "role": "narrow from detector/data context to the isolated prompt-photon physics object",
                "speaker_script": str(photon_script.relative_to(ROOT)),
            },
        ],
        "data_source_evidence": DATA_SOURCE_EVIDENCE,
        "source_basis": [
            "Official BNL RHIC sPHENIX detector imagery for Slide 2.",
            "Yeonju DIS2026 and Hanpu AUM2026 reference decks for p+p data-taking numbers.",
            "Yeonju/Hanpu-style opening progression used as structure only: experiment/data first, photon object second.",
            "No copied Yeonju screenshots, no copied internal figures, no Google Slides mutation.",
        ],
        "comparison_candidates_retained": [
            str((output_dir / "hp2026_slide02_sphenix_real_detector_opener.png").relative_to(ROOT)),
            str((output_dir / "hp2026_slide02_sphenix_hard_probes_photon_baseline.png").relative_to(ROOT)),
            str((output_dir / "hp2026_slide02_why_isolated_prompt_photons.png").relative_to(ROOT)),
        ],
        "notes": [
            "Slide 2 uses only photon-relevant detector callouts: EMCal, HCal, tracking, and MBD/trigger.",
            "Au+Au/O+O detailed running numbers are omitted from the Slide 2 PNG to keep the talk opening focused on p+p.",
            "Slide 3 is photon-only and deliberately avoids a photon+jet narrative this early.",
            "No slide numbers or provenance footers are baked into the PNGs.",
        ],
    }
    manifest_path = output_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return manifest_path


def render_opening_sequence(output_dir: Path) -> Path:
    slide2 = render_detector_data_context(output_dir)
    slide3 = render_isolated_photon_motivation(output_dir)
    write_opening_sequence_manifest(output_dir, slide2, slide3)
    return slide2


def write_real_detector_script(output_dir: Path) -> Path:
    script = """# HP2026 Slide 2 Speaker Script

Before getting into the photon selection itself, I want to place the measurement in the sPHENIX program. sPHENIX is built to make high-rate hard-probe measurements at RHIC: jets, photons, and heavy flavor with modern tracking and calorimetry.

The reason isolated prompt photons are special is that the photon carries information from the hard scattering but does not undergo strong final-state interactions. The isolation requirement makes that statement experimentally useful by suppressing the large decay and fragmentation backgrounds before the purity and correction steps.

So the job of this talk is to connect the detector to one clean baseline result: the isolated prompt-photon cross section in p+p collisions at 200 GeV. Once that baseline is under control, it becomes a reference point for future gamma-jet and heavy-ion measurements.

Transition: I will first define the photon object and isolation requirement, then show how the analysis turns candidates into a corrected cross section.
"""
    path = output_dir / "hp2026_slide02_sphenix_real_detector_opener_script.md"
    path.write_text(script, encoding="utf-8")
    return path


def render_real_detector_opener(output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    detector_img, detector_asset = load_real_detector_asset()
    img = Image.new("RGBA", (W, H), (*SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")

    draw.rectangle((0, 0, W, H), fill=(*SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*PHOTON, 255))

    logo = load_sphenix_logo()
    if logo is not None:
        paste_fit(img, logo, (2070, 84, 2384, 166), anchor="right")

    draw.text((132, 86), REAL_DETECTOR_TITLE, font=font(TIMES_BOLD, 82), fill=INK)
    draw.text((136, 188), REAL_DETECTOR_SUBTITLE, font=font(TIMES_ITALIC, 39), fill=BLUE)
    draw.line((132, 284, W - 132, 284), fill=(221, 226, 232), width=3)

    draw_real_detector_photo_panel(img, detector_img)
    draw_real_detector_cards(img)

    draw.line((132, 1232, 2428, 1232), fill=(217, 225, 233, 255), width=2)
    draw.text((132, 1268), REAL_DETECTOR_BRIDGE, font=font(TIMES_ITALIC, 39), fill=BLUE)

    png = output_dir / "hp2026_slide02_sphenix_real_detector_opener.png"
    img.convert("RGB").save(png, "PNG")
    script_path = write_real_detector_script(output_dir)

    manifest_path = output_dir / "manifest.json"
    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "primary_output": str(png.relative_to(ROOT)),
        "speaker_script": str(script_path.relative_to(ROOT)),
        "size": [W, H],
        "mode": "RGB",
        "title": REAL_DETECTOR_TITLE,
        "slide_role": "first content slide after the title slide",
        "core_claim": "The p+p isolated prompt-photon result is the clean baseline entry point for the sPHENIX hard-probes program.",
        "detector_asset": detector_asset,
        "official_source_urls": REAL_DETECTOR_SOURCE,
        "source_basis": [
            "Official BNL RHIC sPHENIX detector page and BNL-hosted detector image asset.",
            "Reference inventory: HPslides_v1 and current PPG12 material are truth base; Yeonju/Hanpu/Stefan guide pacing.",
            "Deck linearization: establish sPHENIX/hard-probes context before photon definitions and analysis details.",
        ],
        "google_slides_mutation": False,
        "candidate_files": [
            str(png.relative_to(ROOT)),
            str(script_path.relative_to(ROOT)),
            str((output_dir / "hp2026_slide02_sphenix_hard_probes_photon_baseline.png").relative_to(ROOT)),
            str((output_dir / "hp2026_slide02_why_isolated_prompt_photons.png").relative_to(ROOT)),
        ],
        "notes": [
            "No slide number or provenance footer.",
            "Detector image is an official BNL-hosted asset, not AI-generated artwork.",
            "Photon wave and detector highlight are local visual emphasis only, not a precise subsystem annotation.",
            "BDT, purity, unfolding, and detailed detector claims are reserved for later slides.",
        ],
    }
    with manifest_path.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return png


def draw_conceptual_visual(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    panel = (126, 333, 1422, 1112)
    draw.rounded_rectangle(panel, radius=10, fill=(*PANEL, 255), outline=(*PANEL_EDGE, 255), width=2)

    # A quiet RHIC-beam cue leading into one hard scattering.
    y0 = 704
    draw.line((212, y0, 720, y0), fill=(178, 199, 216, 255), width=6)
    draw.line((326, y0 - 62, 435, y0 - 14), fill=(178, 199, 216, 150), width=4)
    draw.line((326, y0 + 62, 435, y0 + 14), fill=(178, 199, 216, 150), width=4)
    draw.ellipse((268, y0 - 48, 364, y0 + 48), fill=(255, 255, 255, 255), outline=(98, 139, 169, 255), width=4)
    draw.ellipse((468, y0 - 48, 564, y0 + 48), fill=(255, 255, 255, 255), outline=(98, 139, 169, 255), width=4)
    draw.text((290, y0 - 19), "p", font=font(TIMES_BOLD, 42), fill=BLUE)
    draw.text((491, y0 - 19), "p", font=font(TIMES_BOLD, 42), fill=BLUE)

    collision = (640, 704)
    for r, alpha in ((70, 36), (46, 58), (24, 96)):
        draw.ellipse((collision[0] - r, collision[1] - r, collision[0] + r, collision[1] + r), fill=(*SPHENIX_BLUE, alpha))
    draw.ellipse((collision[0] - 12, collision[1] - 12, collision[0] + 12, collision[1] + 12), fill=(*SPHENIX_BLUE, 210))
    draw.text((520, 806), "hard scattering", font=font(TIMES_ITALIC, 34), fill=MUTED)

    # Faint colored recoil: present enough for the hard process, not enough to
    # make the opening slide feel like a photon+jet talk.
    for angle, length, width, alpha in ((54, 230, 9, 125), (77, 185, 6, 90), (34, 170, 5, 80)):
        ex = collision[0] + math.cos(math.radians(angle)) * length
        ey = collision[1] + math.sin(math.radians(angle)) * length
        draw.line((collision[0], collision[1], ex, ey), fill=(*TEAL_SOFT, alpha), width=width)
        draw.ellipse((ex - 12, ey - 12, ex + 12, ey + 12), fill=(*TEAL_SOFT, alpha + 35))
    draw.text((771, 914), "recoil activity", font=font(TIMES_ITALIC, 30), fill=(102, 132, 148))

    photon_start = (682, 668)
    impact = (1210, 520)
    dx, dy = impact[0] - photon_start[0], impact[1] - photon_start[1]
    length = math.hypot(dx, dy)
    nx, ny = -dy / length, dx / length
    cone_left = (impact[0] + nx * 145, impact[1] + ny * 145)
    cone_right = (impact[0] - nx * 145, impact[1] - ny * 145)
    alpha_layer(
        base,
        lambda d: d.polygon([photon_start, cone_left, cone_right], fill=(245, 181, 34, 30)),
    )
    draw.line((photon_start, cone_left), fill=(219, 157, 32, 90), width=3)
    draw.line((photon_start, cone_right), fill=(219, 157, 32, 90), width=3)

    wave = feynman_points(photon_start, impact, 17, 8.8, 260)
    glow = Image.new("RGBA", base.size, (0, 0, 0, 0))
    gdraw = ImageDraw.Draw(glow, "RGBA")
    draw_polyline(gdraw, wave, (245, 181, 34, 120), 17)
    glow = glow.filter(ImageFilter.GaussianBlur(7))
    base.alpha_composite(glow)
    draw = ImageDraw.Draw(base, "RGBA")
    draw_polyline(draw, wave, (*PHOTON, 255), 9)
    draw_polyline(draw, wave, (*PHOTON_DARK, 210), 3)

    draw_detector_face(base, impact)
    draw = ImageDraw.Draw(base, "RGBA")
    draw.text((826, 475), "isolated prompt photon", font=font(TIMES_ITALIC, 36), fill=BLUE)
    draw.text((938, 591), "quiet cone", font=font(TIMES_ITALIC, 29), fill=(145, 111, 43))


def draw_reason_cards(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x, y = 1510, 345
    w, h, gap = 870, 216, 34
    for idx, (heading, body) in enumerate(REASONS):
        top = y + idx * (h + gap)
        stripe = SPHENIX_BLUE if idx == 0 else TEAL if idx == 2 else PHOTON
        draw.rounded_rectangle((x, top, x + w, top + h), radius=10, fill=(*CARD, 255), outline=(*CARD_EDGE, 255), width=2)
        draw.rounded_rectangle((x, top, x + 14, top + h), radius=6, fill=(*stripe, 255))
        draw_icon(draw, idx, (x + 86, top + 108))
        draw.text((x + 162, top + 41), heading, font=font(TIMES_BOLD, 43), fill=INK)
        draw_wrapped(draw, body, (x + 164, top + 101), 620, font(TIMES, 34), fill=MUTED, line_gap=11)


def render_why_photons(output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    img = Image.new("RGBA", (W, H), (*SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")

    draw.rectangle((0, 0, W, H), fill=(*SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*PHOTON, 255))

    draw.text((132, 98), TITLE, font=font(TIMES_BOLD, 90), fill=INK)
    draw.text((136, 205), SUBTITLE, font=font(TIMES_ITALIC, 43), fill=BLUE)
    draw.line((132, 292, W - 132, 292), fill=(221, 226, 232), width=3)

    draw_conceptual_visual(img)
    draw_reason_cards(img)

    draw.rounded_rectangle((132, 1192, W - 132, 1282), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw.text((174, 1216), BRIDGE, font=font(TIMES_ITALIC, 40), fill=BLUE)

    png = output_dir / "hp2026_slide02_why_isolated_prompt_photons.png"
    img.convert("RGB").save(png, "PNG")

    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "output": str(png.relative_to(ROOT)),
        "size": [W, H],
        "mode": "RGB",
        "title": TITLE,
        "slide_role": "first content slide after the title slide",
        "core_claim": "Isolated prompt photons provide a clean p+p hard-scattering baseline at RHIC.",
        "source_basis": [
            "HP2026 title-slide visual language: Times New Roman, restrained institutional palette",
            "HP2026 photon reference deck linearization: photon-only opening arc",
            "Conceptual motivation pattern from Yeonju/Hanpu/Stefan references without copied figures",
        ],
        "google_slides_mutation": False,
        "notes": [
            "No slide number or provenance footer.",
            "No copied public/internal plots or figures.",
            "Photon-only motivation; BDT, purity, unfolding, and detector details are left for later slides.",
        ],
    }
    with (output_dir / "manifest.json").open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return png


def render_sphenix_hard_probe_context(output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    img = Image.new("RGBA", (W, H), (*SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")

    draw.rectangle((0, 0, W, H), fill=(*SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*PHOTON, 255))

    logo = load_sphenix_logo()
    if logo is not None:
        paste_fit(img, logo, (2062, 82, 2380, 168), anchor="right")

    draw.text((132, 88), CONTEXT_TITLE, font=font(TIMES_BOLD, 76), fill=INK)
    draw.text((136, 192), CONTEXT_SUBTITLE, font=font(TIMES_ITALIC, 40), fill=BLUE)
    draw.line((132, 292, W - 132, 292), fill=(221, 226, 232), width=3)

    draw_context_visual(img)
    draw_context_cards(img)

    draw.rounded_rectangle((132, 1190, W - 132, 1284), radius=8, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw.text((174, 1217), CONTEXT_BRIDGE, font=font(TIMES_ITALIC, 39), fill=BLUE)

    png = output_dir / "hp2026_slide02_sphenix_hard_probes_photon_baseline.png"
    img.convert("RGB").save(png, "PNG")

    manifest_path = output_dir / "manifest.json"
    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "primary_output": str(png.relative_to(ROOT)),
        "size": [W, H],
        "mode": "RGB",
        "title": CONTEXT_TITLE,
        "slide_role": "first content slide after the title slide",
        "core_claim": "The p+p isolated prompt-photon result is the clean baseline entry point for the sPHENIX hard-probes program.",
        "source_basis": [
            "Reference inventory: HPslides_v1 and current PPG12 material are truth base; Yeonju/Hanpu/Stefan guide pacing.",
            "Deck linearization: establish sPHENIX/data context before photon definitions and analysis details.",
            "Conceptual sPHENIX/RHIC schematic drawn locally; no copied detector figure or result plot.",
        ],
        "google_slides_mutation": False,
        "candidate_files": [
            str(png.relative_to(ROOT)),
            str((output_dir / "hp2026_slide02_why_isolated_prompt_photons.png").relative_to(ROOT)),
        ],
        "notes": [
            "No slide number or provenance footer.",
            "Small sPHENIX logo uses cached official blue logo from the title-slide asset set.",
            "Detector/program visual is a conceptual schematic, not an official detector rendering.",
            "BDT, purity, unfolding, and detailed detector claims are reserved for later slides.",
        ],
    }
    with manifest_path.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return png


def render(output_dir: Path, variant: str) -> Path:
    if variant == "opening-sequence":
        return render_opening_sequence(output_dir)
    if variant == "experiment-dataset-standalone":
        return render_experiment_dataset_standalone(output_dir)
    if variant == "detector-data-progressive-build":
        return render_detector_data_progressive_build(output_dir)
    if variant == "detector-data-context":
        png = render_detector_data_context(output_dir)
        write_opening_sequence_manifest(output_dir, png, output_dir / "hp2026_slide03_isolated_photon_motivation.png")
        return png
    if variant == "isolated-photon-motivation":
        png = render_isolated_photon_motivation(output_dir)
        write_opening_sequence_manifest(output_dir, output_dir / "hp2026_slide02_detector_data_context.png", png)
        return png
    if variant == "real-detector-opener":
        return render_real_detector_opener(output_dir)
    if variant == "why-photons":
        return render_why_photons(output_dir)
    if variant == "sphenix-hard-probes":
        return render_sphenix_hard_probe_context(output_dir)
    raise ValueError(f"unknown variant: {variant}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--variant",
        choices=(
            "opening-sequence",
            "experiment-dataset-standalone",
            "detector-data-progressive-build",
            "detector-data-context",
            "isolated-photon-motivation",
            "real-detector-opener",
            "sphenix-hard-probes",
            "why-photons",
        ),
        default="opening-sequence",
        help="Opening slide concept or two-slide sequence to render.",
    )
    args = parser.parse_args()
    png = render(args.output_dir, args.variant)
    print(png)
    if args.variant == "detector-data-progressive-build":
        print(args.output_dir / "hp2026_slide02_progressive_build_manifest.json")
    elif args.variant == "experiment-dataset-standalone":
        print(args.output_dir / "hp2026_experiment_dataset_standalone_manifest.json")
    else:
        print(args.output_dir / "manifest.json")


if __name__ == "__main__":
    main()

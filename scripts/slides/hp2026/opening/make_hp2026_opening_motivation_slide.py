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
HP2026_FOOTER_RULE = (217, 225, 233)
HP2026_FOOTER_RULE_RGBA = (*HP2026_FOOTER_RULE, 255)

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
CURRENT_SLIDE2_SUBSYSTEM_SCREENSHOT = ASSET_DIR / "hp2026_slide02_current_subsystems_body_20260602.png"
USER_DETECTOR_CLEAN_REFERENCE = ASSET_DIR / "sphenix_detector_clean_reference_user_20260602.png"
USER_DETECTOR_ANNOTATED_REFERENCE = ASSET_DIR / "sphenix_detector_annotated_reference_user_20260602.png"
USER_DETECTOR_TARGET_LAYOUT = ASSET_DIR / "sphenix_detector_target_layout_user_20260602.png"
USER_DETECTOR_THREE_ARROW_TARGET = ASSET_DIR / "sphenix_detector_three_arrow_target_user_20260602.png"
USER_DETECTOR_EXACT_ARROW_LAYOUT = ASSET_DIR / "sphenix_detector_exact_arrow_layout_user_20260602.png"
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

RECREATED_SUBSYSTEM_ARROWS = [
    ("outer HCal", (1604, 178), (1860, 468), "center"),
    ("inner HCal", (1208, 562), (1800, 574), "left"),
    ("MVTX & INTT", (1340, 1204), (1822, 690), "center"),
    ("TPC", (2168, 1216), (2028, 690), "center"),
    ("EMCal", (2388, 1094), (2078, 684), "right"),
    ("solenoid", (2388, 788), (2198, 594), "right"),
]

RECREATED_SUBSYSTEM_ARROWS_USER_IMAGE = [
    ("outer HCal", (0.29, -0.035), (0.365, 0.305), "center"),
    ("inner HCal", (-0.205, 0.315), (0.365, 0.405), "left"),
    ("MVTX & INTT", (-0.090, 0.990), (0.455, 0.525), "center"),
    ("TPC", (0.700, 0.935), (0.575, 0.560), "center"),
    ("EMCal", (0.865, 0.855), (0.605, 0.515), "center"),
    ("solenoid", (0.840, 0.520), (0.705, 0.385), "center"),
]

EXACT_ARROW_VECTOR_OVERLAY = [
    ("outer HCal", (152, 18), (187, 50), (249, 168), "center"),
    ("inner HCal", (26, 188), (24, 207), (205, 207), "left"),
    ("MVTX & INTT", (9, 451), (105, 438), (256, 260), "left"),
    ("TPC", (392, 492), (382, 490), (329, 278), "center"),
    ("EMCal", (435, 452), (438, 438), (358, 274), "center"),
    ("solenoid", (443, 298), (423, 294), (356, 214), "center"),
]

STANDALONE_DATASET_TITLE = "Dataset context: p+p anchors this measurement"
STANDALONE_DATASET_SUBTITLE = (
    "Separate the broad all-subsystem p+p context from the defined Run 24 PPG12 cross-section sample."
)
STANDALONE_PP_AVAILABILITY = [
    ("2024 p+p", "13", "pb^-1", "all subsystems", PHOTON),
    ("2026 p+p", "17", "pb^-1", "all subsystems", PHOTON_DARK),
]
STANDALONE_HEAVY_ION_CONTEXT = [
    ("2025 Au+Au", "6.6 nb^-1", "all subsystems"),
    ("2026 O+O", "23.6 nb^-1", "all subsystems"),
]
STANDALONE_DATASET_EVIDENCE = {
    "values_on_slide": [
        "2024 p+p: 13 pb^-1 all subsystems",
        "2026 p+p: 17 pb^-1 all subsystems",
        "PPG12 Run 24 analysis sample: L = 64.4 pb^-1",
        "2025 Au+Au: 6.6 nb^-1 all subsystems",
        "2026 O+O: 23.6 nb^-1 all subsystems",
    ],
    "source_files": [
        "usefulDocs/20260506_DIS_YeonjuGo.pdf",
        "usefulDocs/sPHENIX_AUM2026_jet_and_photon_measurement.pdf",
        "usefulDocs/sPHENIX_PPG12_Paper_2026-05-21_current_draft.pdf",
        "usefulDocs/PPG12_analysis_note_2026-05-21_v4_current_IAN.pdf",
    ],
    "verification_note": (
        "Yeonju DIS2026 and Hanpu AUM2026 agree on the broad sPHENIX all-subsystem p+p, Au+Au, and O+O data-taking values. "
        "The current PPG12 paper draft and IAN quote the Run 24 isolated prompt-photon analysis luminosity as about 64.4 pb^-1."
    ),
}

PHOTON_MOTIVATION_TITLE = "Why isolated prompt photons?"
PHOTON_MOTIVATION_SUBTITLE = "A color-neutral hard-scattering tag for the p+p baseline at RHIC."
PHOTON_MOTIVATION_REASONS = [
    (
        "Hard scale",
        "The photon is produced in the short-distance scattering and leaves without strong final-state energy loss.",
    ),
    (
        "Isolation",
        "The quiet cone rejects decay-rich and fragmentation-rich activity around the candidate.",
    ),
    (
        "p+p baseline",
        "The corrected cross section anchors the RHIC reference before future heavy-ion photon measurements.",
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
            "slide_2_preserved_body": str(CURRENT_SLIDE2_SUBSYSTEM_SCREENSHOT.relative_to(ROOT)),
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


def draw_standalone_title_only_header(base: Image.Image, title: str) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*PHOTON, 255))
    logo = load_sphenix_logo()
    if logo is not None:
        paste_fit(base, logo, (2188, 58, 2432, 164), anchor="right")
    draw.text((132, 84), title, font=font(TIMES_BOLD, 90), fill=INK)
    draw.line((132, 246, W - 132, 246), fill=(221, 226, 232), width=3)


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


def draw_preserved_subsystems_slide_body(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*PHOTON, 255))

    logo = load_sphenix_logo()
    if logo is not None:
        paste_fit(base, logo, (2068, 54, 2388, 128), anchor="right")

    draw.line((96, 142, W - 96, 142), fill=(224, 229, 235, 255), width=2)

    if not CURRENT_SLIDE2_SUBSYSTEM_SCREENSHOT.exists():
        raise FileNotFoundError(f"Missing preserved Slide 2 screenshot: {CURRENT_SLIDE2_SUBSYSTEM_SCREENSHOT}")

    body = open_rgba(CURRENT_SLIDE2_SUBSYSTEM_SCREENSHOT)
    body = ImageEnhance.Sharpness(body).enhance(1.18)
    body_box = (72, 154, W - 72, H - 50)
    fitted_box = paste_fit_return_box(base, body, body_box)
    x0, y0, x1, y1 = fitted_box
    draw.rectangle((x0, y0, x1, y1), outline=(196, 203, 211, 255), width=2)


def render_standalone_experiment_slide(output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    img = Image.new("RGBA", (W, H), (*SOFT_BG, 255))
    draw_preserved_subsystems_slide_body(img)
    png = output_dir / "hp2026_slide02_sphenix_experiment_for_photons.png"
    img.convert("RGB").save(png, "PNG")
    write_standalone_experiment_script(output_dir)
    return png


def white_to_alpha(img: Image.Image, threshold: int = 246) -> Image.Image:
    rgba = img.convert("RGBA")
    px = rgba.load()
    for y in range(rgba.height):
        for x in range(rgba.width):
            r, g, b, a = px[x, y]
            if r >= threshold and g >= threshold and b >= threshold:
                px[x, y] = (r, g, b, 0)
    return crop_visible(rgba, white_threshold=252)


def white_to_alpha_preserve_canvas(img: Image.Image, threshold: int = 252) -> Image.Image:
    rgba = img.convert("RGBA")
    px = rgba.load()
    for y in range(rgba.height):
        for x in range(rgba.width):
            r, g, b, a = px[x, y]
            if r >= threshold and g >= threshold and b >= threshold:
                px[x, y] = (r, g, b, 0)
    return rgba


def draw_recreated_subsystem_icon(draw: ImageDraw.ImageDraw, kind: str, x: int, y: int) -> None:
    if kind == "tracking":
        for idx, scale in enumerate((0, 1, 2)):
            box = (x + 8 + idx * 8, y + 10 + idx * 28, x + 128 - idx * 8, y + 70 + idx * 28)
            draw.arc(box, start=205, end=330, fill=(64, 68, 74), width=5)
            draw.arc((box[0] + 22, box[1] + 8, box[2] - 22, box[3] - 8), start=205, end=330, fill=(116, 125, 136), width=4)
        for dx in (34, 56, 82, 104):
            draw.line((x + dx, y + 10, x + dx - 20, y + 132), fill=(65, 70, 78), width=4)
            draw.polygon([(x + dx, y + 6), (x + dx - 10, y + 24), (x + dx + 8, y + 22)], fill=(65, 70, 78))
    elif kind == "calorimetry":
        colors = [(211, 218, 203), (227, 205, 164), (173, 207, 222)]
        for idx, color in enumerate(colors):
            top = y + 14 + idx * 34
            poly = [(x + 20, top), (x + 118, top + 32), (x + 74, top + 66), (x - 18, top + 32)]
            draw.polygon(poly, fill=(*color, 255), outline=(45, 50, 56), width=3)
            draw.line((x + 118, top + 32, x + 118, top + 58, x + 74, top + 92), fill=(45, 50, 56), width=3)
            draw.line((x - 18, top + 32, x - 18, top + 58, x + 74, top + 92), fill=(45, 50, 56), width=3)
    else:
        cy = y + 78
        axis = (55, 61, 70)
        tan = (218, 200, 158)
        blue = (172, 204, 218)
        green = (184, 210, 198)
        soft = (*SOFT_BG, 255)

        # Forward detectors sit as annular stations along the beam direction.
        draw.line((x + 4, cy, x + 154, cy), fill=axis, width=5)
        draw.polygon([(x + 4, cy), (x + 22, cy - 10), (x + 22, cy + 10)], fill=axis)
        draw.polygon([(x + 154, cy), (x + 136, cy - 10), (x + 136, cy + 10)], fill=axis)
        draw.ellipse((x + 73, cy - 7, x + 87, cy + 7), fill=(*PHOTON, 255), outline=axis, width=2)

        def endcap(cx: int, fill: tuple[int, int, int], offset: int) -> None:
            body = (cx - 17, cy - 52, cx + 17, cy + 52)
            side = (cx + offset - 17, cy - 52, cx + offset + 17, cy + 52)
            draw.ellipse(side, fill=(*fill, 130), outline=axis, width=2)
            draw.line((cx + offset, cy - 52, cx, cy - 52), fill=axis, width=2)
            draw.line((cx + offset, cy + 52, cx, cy + 52), fill=axis, width=2)
            draw.ellipse(body, fill=(*fill, 230), outline=axis, width=3)
            # Radial sector lines make the face read like a compact segmented
            # forward detector disk without adding the long external cabling.
            rx, ry = 15, 48
            for deg in range(0, 360, 30):
                theta = math.radians(deg)
                ex = cx + rx * math.cos(theta)
                ey = cy + ry * math.sin(theta)
                draw.line((cx, cy, ex, ey), fill=(255, 255, 255, 190), width=2)
            draw.ellipse((cx - 11, cy - 32, cx + 11, cy + 32), outline=(255, 255, 255, 135), width=2)
            draw.ellipse((cx - 7, cy - 22, cx + 7, cy + 22), fill=soft, outline=axis, width=2)
            draw.ellipse((cx - 3, cy - 7, cx + 3, cy + 7), fill=(*PHOTON, 235), outline=axis, width=1)

        endcap(x + 48, blue, -8)
        endcap(x + 112, tan, 8)

        # Small far-forward neutral-energy stations, visually separated from MBD/sEPD.
        draw.rounded_rectangle((x + 4, cy - 24, x + 26, cy + 24), radius=4, fill=(*green, 245), outline=axis, width=2)
        draw.rounded_rectangle((x + 134, cy - 24, x + 156, cy + 24), radius=4, fill=(*green, 245), outline=axis, width=2)


def render_recreated_subsystem_icon(kind: str) -> Image.Image:
    canvas = Image.new("RGBA", (240, 200), (0, 0, 0, 0))
    icon_draw = ImageDraw.Draw(canvas, "RGBA")
    draw_recreated_subsystem_icon(icon_draw, kind, 48, 28)
    return crop_visible(canvas, white_threshold=255)


def draw_recreated_subsystem_text(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    rows = [
        (
            "tracking",
            342,
            SPHENIX_BLUE,
            (235, 247, 254),
            "Tracking system",
            "(MVTX, INTT, TPC, TPOT)",
            ["Silicon vertex detectors & time-projection chamber", "inside a 1.4 T solenoid"],
        ),
        (
            "calorimetry",
            634,
            PHOTON,
            (255, 248, 232),
            "Calorimetry",
            "(EMCal, HCal)",
            ["Electromagnetic & hadronic calorimeters", "(inner/outer HCal)"],
        ),
        (
            "forward",
            926,
            (49, 132, 97),
            (226, 246, 238),
            "Forward detectors",
            "(MBD, sEPD, ZDC)",
            ["Provide minimum-bias triggers, centrality, &", "event-plane information"],
        ),
    ]
    panel_left = 94
    panel_right = 1224
    panel_h = 244
    accent_w = 206
    text_x = panel_left + 258
    heading_font = font(TIMES_BOLD, 46)
    paren_font = font(TIMES, 39)
    body_font = font(TIMES, 40)

    for kind, panel_top, accent, accent_fill, heading, paren, lines in rows:
        panel = (panel_left, panel_top, panel_right, panel_top + panel_h)
        shadow_layer = Image.new("RGBA", base.size, (0, 0, 0, 0))
        shadow_draw = ImageDraw.Draw(shadow_layer, "RGBA")
        shadow_draw.rounded_rectangle(
            (panel[0] + 5, panel[1] + 7, panel[2] + 5, panel[3] + 7),
            radius=12,
            fill=(29, 45, 64, 34),
        )
        shadow_layer = shadow_layer.filter(ImageFilter.GaussianBlur(5))
        base.alpha_composite(shadow_layer)
        draw.rounded_rectangle(panel, radius=12, fill=(255, 255, 255, 246), outline=(210, 222, 232, 255), width=2)
        draw.rounded_rectangle((panel[0] + 1, panel[1] + 1, panel[0] + accent_w, panel[3] - 1), radius=12, fill=(*accent_fill, 235))
        draw.rounded_rectangle((panel[0] + 1, panel[1] + 1, panel[0] + 13, panel[3] - 1), radius=6, fill=(*accent, 235))
        draw.line((panel[0] + accent_w, panel[1] + 20, panel[0] + accent_w, panel[3] - 20), fill=(*PANEL_EDGE, 230), width=2)
        icon = render_recreated_subsystem_icon(kind)
        icon_box = (panel_left + 26, panel_top + 28, panel_left + accent_w - 24, panel_top + panel_h - 28)
        paste_fit(base, icon, icon_box)
        heading_w, _ = text_box(draw, heading, heading_font)
        heading_y = panel_top + 44
        body_y = panel_top + 112
        draw.text((text_x, heading_y), heading, font=heading_font, fill=INK)
        draw.text((text_x + heading_w + 12, heading_y + 5), paren, font=paren_font, fill=INK)
        for idx, line in enumerate(lines):
            draw.text((text_x, body_y + idx * 51), line, font=body_font, fill=INK)


def draw_recreated_detector_arrow(
    draw: ImageDraw.ImageDraw,
    label: str,
    label_xy: tuple[int, int],
    target_xy: tuple[int, int],
    align: str,
) -> None:
    label_font = font(TIMES_BOLD, 42)
    lx, ly = label_xy
    tw, th = text_box(draw, label, label_font)
    if align == "center":
        tx = lx - tw // 2
    elif align == "right":
        tx = lx - tw
    else:
        tx = lx
    draw.text((tx, ly), label, font=label_font, fill=(0, 0, 0), stroke_width=4, stroke_fill=(*SOFT_BG, 245))
    if align == "left":
        start = (tx + tw + 16, ly + th // 2)
    elif align == "right":
        start = (tx - 16, ly + th // 2)
    else:
        start = (lx, ly + th + 12)
        if ly > target_xy[1]:
            start = (lx, ly - 8)
    sx, sy = start
    ex, ey = target_xy
    angle = math.atan2(ey - sy, ex - sx)
    line_end = (ex - math.cos(angle) * 24, ey - math.sin(angle) * 24)
    arrow = (244, 167, 72, 245)
    draw.line((start, line_end), fill=(255, 255, 255, 220), width=15)
    draw.line((start, line_end), fill=arrow, width=8)
    size = 34
    spread = 0.58
    base = (ex - math.cos(angle) * size, ey - math.sin(angle) * size)
    p1 = (base[0] + math.cos(angle + math.pi / 2) * size * spread, base[1] + math.sin(angle + math.pi / 2) * size * spread)
    p2 = (base[0] + math.cos(angle - math.pi / 2) * size * spread, base[1] + math.sin(angle - math.pi / 2) * size * spread)
    draw.polygon([target_xy, p1, p2], fill=(255, 255, 255, 230))
    draw.polygon([target_xy, p1, p2], fill=arrow)


def draw_recreated_detector_arrow_on_image(
    draw: ImageDraw.ImageDraw,
    image_bounds: tuple[int, int, int, int],
    label: str,
    label_rel: tuple[float, float],
    target_rel: tuple[float, float],
    align: str,
) -> None:
    x0, y0, x1, y1 = image_bounds
    iw = x1 - x0
    ih = y1 - y0
    label_xy = (round(x0 + label_rel[0] * iw), round(y0 + label_rel[1] * ih))
    target_xy = (round(x0 + target_rel[0] * iw), round(y0 + target_rel[1] * ih))
    draw_recreated_detector_arrow(draw, label, label_xy, target_xy, align)


def draw_exact_reference_vector_overlay(
    draw: ImageDraw.ImageDraw,
    image_bounds: tuple[int, int, int, int],
    source_size: tuple[int, int],
) -> None:
    x0, y0, x1, y1 = image_bounds
    iw = x1 - x0
    ih = y1 - y0
    sw, sh = source_size

    def map_pt(pt: tuple[int, int]) -> tuple[int, int]:
        return (round(x0 + pt[0] / sw * iw), round(y0 + pt[1] / sh * ih))

    arrow = (244, 167, 72, 252)
    label_font = font(TIMES_BOLD, 46)
    for label, label_pos, start_pos, target_pos, align in EXACT_ARROW_VECTOR_OVERLAY:
        start = map_pt(start_pos)
        target = map_pt(target_pos)
        sx, sy = start
        ex, ey = target
        angle = math.atan2(ey - sy, ex - sx)
        line_end = (ex - math.cos(angle) * 25, ey - math.sin(angle) * 25)
        draw.line((start, line_end), fill=(255, 255, 255, 230), width=17)
        draw.line((start, line_end), fill=arrow, width=9)
        size = 35
        spread = 0.58
        base = (ex - math.cos(angle) * size, ey - math.sin(angle) * size)
        p1 = (base[0] + math.cos(angle + math.pi / 2) * size * spread, base[1] + math.sin(angle + math.pi / 2) * size * spread)
        p2 = (base[0] + math.cos(angle - math.pi / 2) * size * spread, base[1] + math.sin(angle - math.pi / 2) * size * spread)
        draw.polygon([target, p1, p2], fill=(255, 255, 255, 232))
        draw.polygon([target, p1, p2], fill=arrow)

        lx, ly = map_pt(label_pos)
        tw, th = text_box(draw, label, label_font)
        if align == "center":
            tx = lx - tw // 2
        else:
            tx = lx
        draw.text((tx, ly), label, font=label_font, fill=(0, 0, 0), stroke_width=5, stroke_fill=(*SOFT_BG, 245))


def draw_recreated_footer(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    footer_top = 1326
    draw.rectangle((0, footer_top, W, H), fill=(*SOFT_BG, 255))
    draw.line((0, footer_top, W, footer_top), fill=HP2026_FOOTER_RULE_RGBA, width=2)
    asset_dir = TITLE_ASSET_DIR
    illinois = crop_visible(open_rgba(asset_dir / "illinois_logo_fullcolor_rgb.png"), white_threshold=252)
    hp = white_to_alpha(open_rgba(asset_dir / "hp2026_indico_logo.png"), threshold=246)
    illinois = fit(illinois, 54, 62)
    hp = fit(hp, 116, 60)
    cy = 1384
    base.alpha_composite(illinois, (30, cy - illinois.height // 2))
    draw.text((104, cy - 17), "Justin Bennett", font=font(TIMES, 31), fill=(43, 49, 57))
    center = "Hard Probes 2026 / June 24, 2026"
    center_font = font(TIMES_BOLD, 31)
    cw, ch = text_box(draw, center, center_font)
    group_w = hp.width + 20 + cw
    gx = (W - group_w) // 2
    base.alpha_composite(hp, (gx, cy - hp.height // 2))
    draw.text((gx + hp.width + 20, cy - ch // 2 - 1), center, font=center_font, fill=(43, 49, 57))


def write_recreated_subsystems_script(output_dir: Path) -> Path:
    script = """# HP2026 Slide 2 Script - sPHENIX Subsystems

To start, I want to introduce the sPHENIX detector at the subsystem level needed for this photon measurement.

On the left, the detector is grouped into tracking, calorimetry, and forward detectors. The tracking system includes the MVTX, INTT, TPC, and TPOT, which give the charged-particle and vertex context inside the 1.4 T solenoid.

The calorimetry is central for this talk. The EMCal measures electromagnetic showers from photon candidates, while the inner and outer HCal layers measure hadronic activity around the event. Those calorimeter systems are what make photon energy, shower shape, and isolation experimentally accessible.

The forward detectors, including the MBD, sEPD, and ZDC, provide minimum-bias triggering, centrality context in nuclear running, event-plane information, and the normalization handles needed for cross-section measurements.

So the detector picture is simple: tracking anchors the event, EMCal measures the photon candidate, HCal helps characterize nearby activity, and the forward systems provide trigger and normalization context. With the detector established, the next step is the p+p data sample used for the measurement.
"""
    path = output_dir / "hp2026_slide02_sphenix_subsystems_recreated_script.md"
    path.write_text(script, encoding="utf-8")
    return path


def render_recreated_sphenix_subsystems(output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    detector_img = open_rgba(USER_DETECTOR_EXACT_ARROW_LAYOUT)
    detector_asset = str(USER_DETECTOR_EXACT_ARROW_LAYOUT.relative_to(ROOT))
    img = Image.new("RGBA", (W, H), (*SOFT_BG, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.rectangle((0, 0, W, H), fill=(*SOFT_BG, 255))
    draw.rectangle((0, 0, W, 22), fill=(*SPHENIX_BLUE, 255))
    draw.rectangle((0, 22, W, 30), fill=(*PHOTON, 255))
    logo = load_sphenix_logo()
    if logo is not None:
        paste_fit(img, logo, (2188, 58, 2432, 164), anchor="right")
    draw.text((112, 108), "sPHENIX Subsystems", font=font(TIMES_BOLD, 88), fill=(0, 0, 0))
    draw.text(
        (116, 214),
        "sPHENIX is a RHIC experiment at Brookhaven National Laboratory on Long Island, New York, "
        "with full azimuthal coverage and spans |η| < 1.1 in pseudorapidity.",
        font=font(TIMES_ITALIC, 34),
        fill=BLUE,
    )
    draw.line((112, 282, W - 112, 282), fill=(221, 226, 232, 255), width=3)
    draw_recreated_subsystem_text(img)

    detector_box = (1260, 306, 2428, 1320)
    rendering = detector_img.convert("RGBA")
    rendering = rendering.resize((rendering.width * 4, rendering.height * 4), Image.Resampling.LANCZOS)
    rendering = ImageEnhance.Color(rendering).enhance(1.06)
    rendering = ImageEnhance.Contrast(rendering).enhance(1.10)
    rendering = rendering.filter(ImageFilter.UnsharpMask(radius=1.4, percent=170, threshold=2))
    rendering = white_to_alpha_preserve_canvas(rendering, threshold=253)
    image_bounds = paste_fit_return_box(img, rendering, detector_box)
    draw_recreated_footer(img)
    png = output_dir / "hp2026_slide02_sphenix_subsystems_recreated_enhanced.png"
    img.convert("RGB").save(png, "PNG")
    script_path = write_recreated_subsystems_script(output_dir)
    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "google_slides_mutation": False,
        "output": str(png.relative_to(ROOT)),
        "speaker_script": str(script_path.relative_to(ROOT)),
        "size": [W, H],
        "mode": "RGB",
        "detector_asset": detector_asset,
        "style_notes": [
            "Recreated from scratch at 2560x1440 using the HP2026 soft background rather than pure white.",
            "Detector panel uses Justin's exact annotated screenshot for the detector diagram so the subsystem arrows match the requested reference layout.",
            "The exact annotated detector screenshot is upscaled and sharpened before placement to improve audience-facing readability without moving arrow targets.",
            "Forward-detector pictogram uses a beamline with annular endcap stations plus far-forward blocks to evoke MBD/sEPD/ZDC geometry.",
            "The standardized HP2026 identity footer is drawn after the original slide body so the title, subsystem text, icons, detector labels, arrows, and detector placement remain unchanged.",
        ],
        "reference_sources": [
            {
                "file": "usefulDocs/The_sPHENIX_Detector_-_sPHENIX_Russia_Workshop.pdf",
                "page": 3,
                "use": "detector subsystem label and arrow-target reference: outer HCal, inner HCal, MVTX+INTT, TPC, EMCal, solenoid",
            },
            {
                "file": str(USER_DETECTOR_CLEAN_REFERENCE.relative_to(ROOT)),
                "use": "user-provided clean detector screenshot saved as a stable local visual reference",
            },
            {
                "file": str(USER_DETECTOR_ANNOTATED_REFERENCE.relative_to(ROOT)),
                "use": "user-provided annotated detector screenshot used to retune oHCal/iHCal/solenoid/EMCal/TPC/MVTX+INTT arrow targets",
            },
            {
                "file": str(USER_DETECTOR_TARGET_LAYOUT.relative_to(ROOT)),
                "use": "user-provided target layout screenshot used to match final arrow style and placement",
            },
            {
                "file": str(USER_DETECTOR_THREE_ARROW_TARGET.relative_to(ROOT)),
                "use": "user-provided crop used to correct outer HCal, inner HCal, and solenoid arrow targets",
            },
            {
                "file": str(USER_DETECTOR_EXACT_ARROW_LAYOUT.relative_to(ROOT)),
                "use": "user-provided exact annotated detector diagram used directly in the final slide panel",
            }
        ],
    }
    manifest_path = output_dir / "hp2026_slide02_sphenix_subsystems_recreated_manifest.json"
    with manifest_path.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    return png


def draw_dataset_primary_panel(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    panel = (132, 292, 1518, 742)
    draw.rounded_rectangle(panel, radius=12, fill=(255, 255, 255, 255), outline=(216, 225, 234, 255), width=2)
    draw.ellipse((178, 326, 226, 374), fill=(*PHOTON, 255))
    draw.text((194, 334), "1", font=font(TIMES_BOLD, 30), fill=(255, 255, 255))
    draw.text((250, 324), "All-subsystem p+p running context", font=font(TIMES_BOLD, 47), fill=INK)
    draw_wrapped(
        draw,
        "Detector-complete p+p data establish the RHIC baseline context for this photon measurement.",
        (252, 388),
        1100,
        font(TIMES_ITALIC, 30),
        fill=MUTED,
        line_gap=7,
    )

    card_y = 508
    card_w = 392
    gap = 46
    for idx, (run, value, unit, detail, color) in enumerate(STANDALONE_PP_AVAILABILITY):
        x = 178 + idx * (card_w + gap)
        draw.rounded_rectangle((x, card_y, x + card_w, card_y + 154), radius=10, fill=(247, 250, 252, 255), outline=(222, 229, 236, 255), width=2)
        draw.rectangle((x, card_y, x + 12, card_y + 154), fill=(*color, 255))
        draw.text((x + 34, card_y + 28), run, font=font(TIMES_BOLD, 34), fill=INK)
        value_font = font(TIMES_BOLD, 69)
        draw.text((x + 202, card_y + 18), value, font=value_font, fill=INK)
        draw_rich_text(draw, (x + 294, card_y + 60), unit, font(TIMES, 30), BLUE)
        draw.text((x + 34, card_y + 102), detail, font=font(TIMES_ITALIC, 25), fill=MUTED)

    summary_x = 178 + 2 * (card_w + gap)
    draw.rounded_rectangle((summary_x, card_y, summary_x + 432, card_y + 154), radius=10, fill=(255, 251, 239, 255), outline=(239, 223, 184, 255), width=2)
    draw.text((summary_x + 34, card_y + 26), "combined context", font=font(TIMES_BOLD, 29), fill=MUTED)
    draw.text((summary_x + 34, card_y + 70), "30", font=font(TIMES_BOLD, 62), fill=INK)
    draw_rich_text(draw, (summary_x + 122, card_y + 105), "pb^-1", font(TIMES, 28), BLUE)
    draw.text((summary_x + 218, card_y + 88), "all-subsystem p+p", font=font(TIMES_ITALIC, 26), fill=MUTED)


def draw_dataset_analysis_sample(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    panel = (132, 802, 1518, 1266)
    draw.rounded_rectangle(panel, radius=12, fill=(239, 246, 250, 255), outline=(205, 221, 232, 255), width=2)
    draw.ellipse((178, 842, 226, 890), fill=(*TEAL, 255))
    draw.text((194, 850), "2", font=font(TIMES_BOLD, 30), fill=(255, 255, 255))
    draw.text((250, 838), "Defined PPG12 cross-section sample", font=font(TIMES_BOLD, 47), fill=INK)
    draw_wrapped(
        draw,
        "The result in this talk uses one defined Run 24 p+p analysis sample, not the full run-context number above.",
        (252, 902),
        1110,
        font(TIMES_ITALIC, 30),
        fill=MUTED,
        line_gap=7,
    )
    draw.rounded_rectangle((184, 1018, 1490, 1196), radius=12, fill=(255, 255, 255, 255), outline=(205, 221, 232, 255), width=2)
    draw.text((226, 1044), "analysis luminosity", font=font(TIMES_BOLD, 33), fill=LIGHT_MUTED)
    draw_rich_text(draw, (226, 1092), "L = 64.4 pb^-1", font(TIMES_BOLD, 78), BLUE)
    draw.text((820, 1099), "used for the isolated prompt-photon", font=font(TIMES_ITALIC, 31), fill=MUTED)
    draw.text((820, 1138), "cross-section measurement", font=font(TIMES_ITALIC, 31), fill=MUTED)
    draw.rounded_rectangle((184, 1220, 1490, 1252), radius=7, fill=(226, 240, 248, 255))
    draw.text((210, 1224), "Key distinction: run context and analysis luminosity are quoted for different purposes.", font=font(TIMES_ITALIC, 24), fill=BLUE)


def draw_dataset_side_context(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")

    def draw_centered_rich_line(
        runs: list[tuple[str, ImageFont.ImageFont, tuple[int, int, int]]],
        y: int,
        x0: int = 1646,
        x1: int = 2388,
    ) -> None:
        total = sum(text_box(draw, text, fnt)[0] for text, fnt, _color in runs)
        x = x0 + ((x1 - x0) - total) / 2
        for text, fnt, color in runs:
            draw.text((x, y), text, font=fnt, fill=color)
            x += text_box(draw, text, fnt)[0]

    panel = (1600, 292, 2428, 1266)
    draw.rounded_rectangle(panel, radius=12, fill=(255, 255, 255, 255), outline=(216, 225, 234, 255), width=2)
    draw.ellipse((1644, 326, 1692, 374), fill=(*BLUE, 255))
    draw.text((1660, 334), "3", font=font(TIMES_BOLD, 30), fill=(255, 255, 255))
    draw.text((1716, 324), "Nuclear-system data collected", font=font(TIMES_BOLD, 40), fill=INK)
    draw_wrapped(
        draw,
        "Luminosity context for the collected Au+Au and O+O samples that motivate the next stage of the photon program.",
        (1648, 390),
        714,
        font(TIMES_ITALIC, 32),
        fill=MUTED,
        line_gap=8,
    )

    y = 560
    for idx, (system, value, detail) in enumerate(STANDALONE_HEAVY_ION_CONTEXT):
        color = LIGHT_MUTED if idx else TEAL_SOFT
        draw.rounded_rectangle((1646, y, 2388, y + 174), radius=10, fill=(249, 250, 251, 255), outline=(224, 230, 236, 255), width=1)
        draw.rectangle((1646, y, 1662, y + 174), fill=(*color, 170))
        draw.text((1690, y + 34), system, font=font(TIMES_BOLD, 37), fill=INK)
        draw_rich_text(draw, (1996, y + 28), value, font(TIMES_BOLD, 45), BLUE if idx == 0 else MUTED)
        draw.text((1690, y + 104), detail, font=font(TIMES_ITALIC, 30), fill=LIGHT_MUTED)
        y += 216

    draw.rounded_rectangle((1646, 1004, 2388, 1210), radius=10, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    normal = font(TIMES, 31)
    bold = font(TIMES_BOLD, 31)
    draw_centered_rich_line([("p+p", bold, BLUE), (" serves as the reference baseline", normal, INK)], 1034)
    draw_centered_rich_line([("for future ", normal, INK), ("A+A", bold, BLUE), (" analyses", normal, INK)], 1076)
    draw_centered_rich_line([("using the photon as a ", normal, INK), ("color-neutral tag", bold, BLUE)], 1118)
    draw_centered_rich_line([("to study the ", normal, INK), ("QGP", bold, BLUE)], 1160)


def draw_dataset_flow_connectors(base: Image.Image) -> None:
    return


def render_standalone_dataset_slide(output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    img = Image.new("RGBA", (W, H), (*SOFT_BG, 255))
    draw_standalone_title_only_header(img, STANDALONE_DATASET_TITLE)
    draw_dataset_primary_panel(img)
    draw_dataset_analysis_sample(img)
    draw_dataset_side_context(img)
    draw_dataset_flow_connectors(img)
    draw_recreated_footer(img)
    png = output_dir / "hp2026_slide03_pp_dataset_context.png"
    img.convert("RGB").save(png, "PNG")
    write_standalone_dataset_script(output_dir)
    return png


def write_standalone_experiment_script(output_dir: Path) -> Path:
    script = """# HP2026 Slide 2 Script - sPHENIX Experiment

To start, I want to introduce the sPHENIX detector at the level needed for this photon measurement.

The main subsystems are grouped here into tracking, calorimetry, and forward detectors. The tracking system includes the MVTX, INTT, TPC, and TPOT. These are the silicon vertex detectors and time-projection chamber systems inside the 1.4 T solenoid, and they provide the charged-particle and vertex context for the event.

The calorimetry is the key part for this talk. The EMCal, or electromagnetic calorimeter, measures electromagnetic showers from photon candidates. The HCal, or hadronic calorimeter, sits around it in inner and outer layers and measures hadronic activity. Together, those calorimeters give both the photon energy measurement and the surrounding activity needed later for isolation.

The forward detectors, including the MBD, sEPD, and ZDC, provide the minimum-bias trigger, centrality information in heavy-ion running, event-plane information, and the luminosity handles needed to normalize measurements.

So the detector picture I want the audience to keep in mind is simple: sPHENIX gives us tracking for event context, EMCal for the photon candidate, HCal for the nearby hadronic activity, and forward systems for triggering and normalization. With the experiment oriented, the next step is to say what p+p data this measurement is built from.
"""
    path = output_dir / "hp2026_slide02_sphenix_experiment_for_photons_script.md"
    path.write_text(script, encoding="utf-8")
    return path


def write_standalone_dataset_script(output_dir: Path) -> Path:
    script = """# HP2026 Slide 3 Script - p+p Dataset Context

Now that the detector context is set, I want to make the data story very explicit, because there are two numbers here that mean different things.

First, the top-left panel is the broad all-subsystem p+p context. In 2024, sPHENIX recorded 13 inverse picobarns with all subsystems, and in 2026 there is another 17 inverse picobarns with all subsystems. Together, that gives 30 inverse picobarns of detector-complete p+p context around the photon program.

Second, the lower-left panel is the defined PPG12 sample used for the isolated prompt-photon cross section. That is the Run 24 analysis luminosity, 64.4 inverse picobarns. So when I quote 64.4 inverse picobarns later, that number belongs to the measured cross section, not to the broad all-subsystem run-context summary.

The right side is just a quick declaration of the nuclear-system data already collected. We have 2025 Au+Au data with 6.6 inverse nanobarns using all subsystems, and 2026 O+O data with 23.6 inverse nanobarns using all subsystems. I do not want to dwell on those samples here; they are included to show where this baseline goes next.

So the takeaway is simple: the analysis result I am presenting is p+p, and it uses one clearly defined Run 24 luminosity. That p+p cross section serves as the reference baseline for future A+A analyses using the photon as a color-neutral tag to study the quark-gluon plasma. With the detector and dataset established, the next step is to define the object we are measuring: the isolated prompt photon.
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
            "Slide 2 preserves the current HPslides_v1 sPHENIX Subsystems slide screenshot as the body content.",
            "Only the HP2026-style top stripe, sPHENIX mark, and standardized identity footer are added around the preserved Slide 2 body so it matches the rest of the deck without redesigning the subsystem content.",
            "The subsystem text, icons, detector labels, and arrows are inherited from the screenshot rather than redrawn; the old baked-in body slide number is visually suppressed to keep the footer consistent.",
            "Slide 3 introduces the dataset only: broad all-subsystem p+p availability, prominent PPG12 Run 24 cross-section luminosity, and muted all-subsystem heavy-ion context.",
            "Slide 3 deliberately shows only all-subsystem context and the PPG12 Run 24 result luminosity.",
            "Slide 3 uses a title-only header with larger type and expanded content panels for consistency with the preceding sPHENIX subsystem slide.",
            "Slide 3 uses the same standardized HP2026 identity footer as the rest of the generated talk slides; the dataset content is otherwise unchanged.",
            "Slide 4 can remain the isolated prompt-photon motivation slide, preserving a clean progression: experiment -> data -> object.",
            "No provenance footer or internal note is baked into either PNG.",
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
    panel = (126, 333, 1422, 1248)
    draw.rounded_rectangle(panel, radius=10, fill=(*PANEL, 255), outline=(*PANEL_EDGE, 255), width=2)

    def arrow_line(
        start: tuple[float, float],
        end: tuple[float, float],
        color: tuple[int, int, int],
        width: int = 4,
        alpha: int = 230,
        head: int = 14,
    ) -> None:
        sx, sy = start
        ex, ey = end
        draw.line((sx, sy, ex, ey), fill=(*color, alpha), width=width)
        angle = math.atan2(ey - sy, ex - sx)
        left = (ex - head * math.cos(angle - 0.55), ey - head * math.sin(angle - 0.55))
        right = (ex - head * math.cos(angle + 0.55), ey - head * math.sin(angle + 0.55))
        draw.polygon([(ex, ey), left, right], fill=(*color, alpha))

    def label(text: str, xy: tuple[int, int], size: int = 23, fill=MUTED, rich: bool = False) -> None:
        if rich:
            draw_rich_text(draw, xy, text, font(TIMES_ITALIC, size), fill)
        else:
            draw.text(xy, text, font=font(TIMES_ITALIC, size), fill=fill)

    def draw_small_feynman(box: tuple[int, int, int, int], title: str, equation: str, mode: str) -> None:
        x0, y0, x1, y1 = box
        draw.rounded_rectangle(box, radius=8, fill=(255, 255, 255, 238), outline=(222, 229, 236, 255), width=2)
        draw.text((x0 + 22, y0 + 18), title, font=font(TIMES_BOLD, 25), fill=INK)
        draw.text((x0 + 22, y0 + 54), equation, font=font(TIMES, 23), fill=BLUE)

        cx, cy = x0 + 160, y0 + 132
        if mode == "compton":
            arrow_line((cx - 95, cy - 64), (cx - 10, cy - 12), INK, width=3, alpha=210, head=10)
            pts = feynman_points((cx - 96, cy + 58), (cx - 8, cy + 14), 6, 5.0, 90)
            draw_polyline(draw, pts, (*TEAL, 210), 4)
            arrow_line((cx + 8, cy - 10), (cx + 104, cy - 58), PHOTON_DARK, width=4, alpha=235, head=11)
            arrow_line((cx + 6, cy + 12), (cx + 104, cy + 58), INK, width=3, alpha=205, head=10)
            labels = [("q", cx - 118, cy - 82), ("g", cx - 118, cy + 48), ("γ", cx + 112, cy - 76), ("q", cx + 112, cy + 48)]
        else:
            arrow_line((cx - 96, cy - 62), (cx - 10, cy - 12), INK, width=3, alpha=210, head=10)
            arrow_line((cx - 96, cy + 62), (cx - 10, cy + 12), INK, width=3, alpha=210, head=10)
            arrow_line((cx + 8, cy - 10), (cx + 104, cy - 58), PHOTON_DARK, width=4, alpha=235, head=11)
            pts = feynman_points((cx + 8, cy + 14), (cx + 104, cy + 58), 6, 5.0, 90)
            draw_polyline(draw, pts, (*TEAL, 210), 4)
            labels = [("q", cx - 118, cy - 82), ("q̄", cx - 126, cy + 48), ("γ", cx + 112, cy - 76), ("g", cx + 112, cy + 48)]
        draw.ellipse((cx - 12, cy - 12, cx + 12, cy + 12), fill=(*BLUE, 220))
        for txt, lx, ly in labels:
            draw.text((lx, ly), txt, font=font(TIMES_ITALIC, 21), fill=MUTED)

    def draw_direct_gamma_raa_inset(box: tuple[int, int, int, int]) -> None:
        x0, y0, x1, y1 = box
        draw.rounded_rectangle(box, radius=8, fill=(255, 255, 255, 238), outline=(222, 229, 236, 255), width=2)
        draw.text((x0 + 22, y0 + 16), "Color-neutral check", font=font(TIMES_BOLD, 25), fill=INK)
        draw.text((x0 + 22, y0 + 48), "PHENIX direct γ in Au+Au", font=font(TIMES_ITALIC, 21), fill=MUTED)

        px0, py0, px1, py1 = x0 + 78, y0 + 84, x1 - 38, y1 - 42
        draw.rectangle((px0, py0, px1, py1), fill=(255, 255, 255, 255), outline=(80, 84, 90, 255), width=2)
        for frac in (0.25, 0.5, 0.75):
            yy = py1 - frac * (py1 - py0)
            draw.line((px0, yy, px1, yy), fill=(230, 234, 238, 255), width=1)
        y_one = py1 - 0.5 * (py1 - py0)
        draw.line((px0, y_one, px1, y_one), fill=(50, 50, 50, 210), width=3)

        band = [
            (px0 + 78, y_one - 34),
            (px0 + 150, y_one - 46),
            (px0 + 226, y_one - 32),
            (px0 + 308, y_one - 44),
            (px0 + 394, y_one - 24),
            (px0 + 394, y_one + 28),
            (px0 + 308, y_one + 34),
            (px0 + 226, y_one + 42),
            (px0 + 150, y_one + 30),
            (px0 + 78, y_one + 44),
        ]
        draw.polygon(band, fill=(180, 186, 194, 92))

        points = [
            (0.20, 0.50), (0.30, 0.48), (0.39, 0.52), (0.48, 0.49),
            (0.57, 0.53), (0.67, 0.50), (0.77, 0.69), (0.88, 0.51),
        ]
        for fx, fy in points:
            x = px0 + fx * (px1 - px0)
            y = py1 - fy * (py1 - py0)
            draw.line((x, y - 18, x, y + 18), fill=(70, 70, 70, 185), width=2)
            draw.ellipse((x - 6, y - 6, x + 6, y + 6), fill=(35, 35, 35, 240))

        draw_rich_text(draw, (px0 - 52, py0 + 4), "R_AA", font(TIMES_BOLD, 21), INK)
        draw.text((px0 + 8, py0 + 8), "√sNN = 200 GeV", font=font(TIMES_BOLD, 18), fill=INK)
        draw.text((px0 + 8, py0 + 32), "Au+Au", font=font(TIMES, 17), fill=INK)
        draw_rich_text(draw, (px1 - 98, py1 + 8), "p_T", font(TIMES, 20), INK)
        draw_rich_text(draw, (px0 + 24, py1 + 8), "direct γ R_AA ≈ 1", font(TIMES_ITALIC, 20), BLUE)

    def draw_prompt_definition_box(box: tuple[int, int, int, int]) -> None:
        x0, y0, x1, y1 = box
        draw.rounded_rectangle(box, radius=8, fill=(255, 255, 255, 238), outline=(222, 229, 236, 255), width=2)
        draw.text((x0 + 22, y0 + 18), "Prompt photon definition", font=font(TIMES_BOLD, 24), fill=INK)

        y = y0 + 62
        draw.rounded_rectangle((x0 + 24, y, x0 + 118, y + 34), radius=6, fill=(*PHOTON, 235))
        draw.text((x0 + 42, y + 6), "direct", font=font(TIMES_BOLD, 20), fill=INK)
        draw.text((x0 + 128, y + 6), "+", font=font(TIMES_BOLD, 22), fill=MUTED)
        draw.rounded_rectangle((x0 + 154, y, x1 - 24, y + 34), radius=6, fill=(*TEAL, 220))
        draw.text((x0 + 170, y + 6), "fragmentation", font=font(TIMES_BOLD, 19), fill=(255, 255, 255))

        draw.text((x0 + 44, y + 56), "= prompt photons", font=font(TIMES_BOLD, 24), fill=BLUE)
        draw.rounded_rectangle((x0 + 24, y + 104, x0 + 48, y + 128), radius=5, fill=(*LIGHT_MUTED, 220))
        draw.text((x0 + 64, y + 101), "π⁰/η decays", font=font(TIMES_BOLD, 20), fill=MUTED)
        draw.text((x0 + 190, y + 103), "backgrounds", font=font(TIMES_ITALIC, 19), fill=MUTED)

    draw.text((184, 382), "RHIC p+p hard scattering", font=font(TIMES_ITALIC, 33), fill=MUTED)
    draw.text((184, 424), "√s = 200 GeV: direct photon tags the short-distance parton scattering", font=font(TIMES, 25), fill=LIGHT_MUTED)
    y0 = 646
    for x, label in ((300, "p"), (514, "p")):
        draw.ellipse((x - 48, y0 - 48, x + 48, y0 + 48), fill=(255, 255, 255, 255), outline=(98, 139, 169, 255), width=4)
        draw.text((x - 13, y0 - 23), label, font=font(TIMES_BOLD, 45), fill=BLUE)
    draw.line((348, y0, 628, y0), fill=(178, 199, 216, 255), width=6)
    draw.line((466, y0, 628, y0), fill=(178, 199, 216, 255), width=6)

    collision = (666, y0)
    for r, alpha in ((78, 26), (52, 44), (25, 96)):
        draw.ellipse((collision[0] - r, collision[1] - r, collision[0] + r, collision[1] + r), fill=(*SPHENIX_BLUE, alpha))
    draw.ellipse((collision[0] - 12, collision[1] - 12, collision[0] + 12, collision[1] + 12), fill=(*SPHENIX_BLUE, 220))
    draw.text((570, 714), "hard q/g vertex", font=font(TIMES_ITALIC, 24), fill=BLUE)

    # Muted colored recoil and decay-rich activity provide context without
    # making the opening a photon+jet slide.
    for angle, length, width, alpha in ((112, 164, 5, 66), (242, 148, 5, 60), (302, 142, 4, 50)):
        ex = collision[0] + math.cos(math.radians(angle)) * length
        ey = collision[1] + math.sin(math.radians(angle)) * length
        draw.line((collision[0], collision[1], ex, ey), fill=(*TEAL_SOFT, alpha), width=width)
        draw.ellipse((ex - 9, ey - 9, ex + 9, ey + 9), fill=(*TEAL_SOFT, alpha + 25))

    decay_center = (484, 518)
    draw.ellipse((decay_center[0] - 18, decay_center[1] - 18, decay_center[0] + 18, decay_center[1] + 18), fill=(170, 178, 188, 120), outline=(120, 130, 142, 150), width=2)
    draw_rich_text(draw, (decay_center[0] - 30, decay_center[1] - 56), "π^0, η decays", font(TIMES_ITALIC, 21), LIGHT_MUTED)
    for end in ((552, 486), (570, 540)):
        pts = feynman_points(decay_center, end, 4, 3.2, 80)
        draw_polyline(draw, pts, (150, 158, 168, 120), 3)

    photon_start = (698, 618)
    cone_tip = (1122, 540)
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
    draw.arc((cone_tip[0] - 178, cone_tip[1] - 178, cone_tip[0] + 178, cone_tip[1] + 178), start=346, end=44, fill=(*PHOTON_DARK, 55), width=3)

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

    draw.text((792, 462), "direct prompt photon", font=font(TIMES_ITALIC, 35), fill=BLUE)
    draw.text((924, 688), "isolation cone", font=font(TIMES_ITALIC, 29), fill=(145, 111, 43))
    draw.text((1212, 720), "EMCal", font=font(TIMES_ITALIC, 28), fill=LIGHT_MUTED)

    bg_box = (214, 868, 1284, 1178)
    draw.rounded_rectangle(bg_box, radius=8, fill=(255, 255, 255, 238), outline=(222, 229, 236, 255), width=2)
    draw.text((250, 900), "Why the photon is a clean tag", font=font(TIMES_BOLD, 32), fill=INK)
    draw.text((250, 938), "definition plus empirical color-neutral behavior at RHIC", font=font(TIMES_ITALIC, 23), fill=MUTED)
    draw_prompt_definition_box((246, 980, 566, 1158))
    draw_direct_gamma_raa_inset((596, 980, 1250, 1158))
    draw.rounded_rectangle((250, 1190, 1246, 1226), radius=5, fill=(232, 240, 247, 255))
    draw.text((270, 1198), "Direct photons are not strongly quenched; p+p provides the baseline that makes this comparison meaningful.", font=font(TIMES_ITALIC, 21), fill=BLUE)


def draw_photon_motivation_cards(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    x, y = 1510, 345
    w, h, gap = 870, 226, 32
    accents = [PHOTON, SPHENIX_BLUE, TEAL]
    for idx, (heading, body) in enumerate(PHOTON_MOTIVATION_REASONS):
        top = y + idx * (h + gap)
        draw.rounded_rectangle((x, top, x + w, top + h), radius=10, fill=(*CARD, 255), outline=(*CARD_EDGE, 255), width=2)
        draw.rounded_rectangle((x, top, x + 14, top + h), radius=6, fill=(*accents[idx], 255))
        draw_icon(draw, idx, (x + 86, top + 108))
        draw.text((x + 162, top + 38), heading, font=font(TIMES_BOLD, 43), fill=INK)
        draw_wrapped(draw, body, (x + 164, top + 98), 625, font(TIMES, 31), fill=MUTED, line_gap=8)


def draw_photon_motivation_payoff(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")
    box = (1510, 1108, 2380, 1290)
    draw.rounded_rectangle(box, radius=10, fill=(239, 246, 250, 255), outline=(213, 226, 235, 255), width=2)
    draw.rectangle((box[0], box[1], box[0] + 12, box[3]), fill=(*BLUE, 255))
    draw.text((box[0] + 42, box[1] + 28), "Measurement logic", font=font(TIMES_BOLD, 34), fill=INK)
    draw_wrapped(
        draw,
        "Color neutral gives the hard scale; isolation makes the candidate interpretable; p+p makes it the reference.",
        (box[0] + 42, box[1] + 76),
        box[2] - box[0] - 76,
        font(TIMES, 28),
        fill=MUTED,
        line_gap=6,
    )


def draw_prompt_photon_integrated_slide(base: Image.Image) -> None:
    draw = ImageDraw.Draw(base, "RGBA")

    def arrow_line(
        start: tuple[float, float],
        end: tuple[float, float],
        color: tuple[int, int, int],
        width: int = 4,
        alpha: int = 230,
        head: int = 13,
    ) -> None:
        sx, sy = start
        ex, ey = end
        draw.line((sx, sy, ex, ey), fill=(*color, alpha), width=width)
        angle = math.atan2(ey - sy, ex - sx)
        left = (ex - head * math.cos(angle - 0.55), ey - head * math.sin(angle - 0.55))
        right = (ex - head * math.cos(angle + 0.55), ey - head * math.sin(angle + 0.55))
        draw.polygon([(ex, ey), left, right], fill=(*color, alpha))

    def draw_clean_feynman(cx: int, cy: int, scale: float, mode: str) -> None:
        node_top = (cx, cy - 42 * scale)
        node_bot = (cx, cy + 64 * scale)
        if mode == "fragmentation":
            node_top = (cx, cy - 4 * scale)
            node_bot = (cx, cy + 86 * scale)
            arrow_line((cx - 108 * scale, cy - 88 * scale), node_top, INK, width=3, alpha=215, head=10)
            arrow_line(node_top, (cx + 112 * scale, cy - 82 * scale), INK, width=3, alpha=215, head=10)
            pts = feynman_points((cx + 46 * scale, cy - 38 * scale), (cx + 116 * scale, cy + 8 * scale), 5.0 * scale, 4.0, 80)
            draw_polyline(draw, pts, (*PHOTON_DARK, 235), max(3, round(4 * scale)))
            pts = feynman_points(node_bot, node_top, 5.0 * scale, 5.0, 90)
            draw_polyline(draw, pts, (*TEAL, 230), max(3, round(4 * scale)))
            arrow_line((cx - 108 * scale, cy + 150 * scale), node_bot, INK, width=3, alpha=215, head=10)
            arrow_line(node_bot, (cx + 112 * scale, cy + 150 * scale), INK, width=3, alpha=215, head=10)
            labels = [
                ("q", cx - 132 * scale, cy - 110 * scale),
                ("q", cx + 124 * scale, cy - 110 * scale),
                ("γ", cx + 122 * scale, cy + 8 * scale),
                ("q", cx - 132 * scale, cy + 116 * scale),
                ("q", cx + 124 * scale, cy + 116 * scale),
            ]
        else:
            arrow_line((cx - 108 * scale, cy - 112 * scale), node_top, INK, width=3, alpha=215, head=10)
            pts = feynman_points(node_top, (cx + 112 * scale, cy - 112 * scale), 5.5 * scale, 4.3, 80)
            draw_polyline(draw, pts, (*PHOTON_DARK, 235), max(3, round(4 * scale)))
            arrow_line(node_top, node_bot, INK, width=3, alpha=210, head=9)
            if mode == "compton":
                pts = feynman_points((cx - 112 * scale, cy + 128 * scale), node_bot, 5.5 * scale, 5.0, 90)
                draw_polyline(draw, pts, (*TEAL, 230), max(3, round(4 * scale)))
                labels = [
                    ("q", cx - 132 * scale, cy - 134 * scale),
                    ("γ", cx + 116 * scale, cy - 132 * scale),
                    ("g", cx - 128 * scale, cy + 112 * scale),
                    ("q", cx + 124 * scale, cy + 112 * scale),
                ]
            else:
                arrow_line((cx - 112 * scale, cy + 128 * scale), node_bot, INK, width=3, alpha=215, head=10)
                labels = [
                    ("q", cx - 132 * scale, cy - 134 * scale),
                    ("γ", cx + 116 * scale, cy - 132 * scale),
                    ("q", cx - 132 * scale, cy + 112 * scale),
                    ("g", cx + 116 * scale, cy + 112 * scale),
                ]
                pts = feynman_points(node_bot, (cx + 112 * scale, cy + 128 * scale), 5.5 * scale, 5.0, 90)
                draw_polyline(draw, pts, (*TEAL, 230), max(3, round(4 * scale)))
            arrow_line(node_bot, (cx + 112 * scale, cy + 128 * scale), INK, width=3, alpha=215, head=10) if mode == "compton" else None

        for node in (node_top, node_bot):
            draw.ellipse((node[0] - 7 * scale, node[1] - 7 * scale, node[0] + 7 * scale, node[1] + 7 * scale), fill=INK)
        label_font = font(TIMES_ITALIC, max(28, round(42 * scale)))
        for txt, lx, ly in labels:
            label_color = PHOTON_DARK if txt == "γ" else TEAL if txt == "g" else INK
            # Particle labels need a small halo to survive projection at the back of the room.
            for dx, dy in ((-2, 0), (2, 0), (0, -2), (0, 2), (-1, -1), (1, 1)):
                draw.text((lx + dx, ly + dy), txt, font=label_font, fill=(255, 255, 255, 235))
            draw.text((lx, ly), txt, font=label_font, fill=label_color)

    def draw_raa_fallback(box: tuple[int, int, int, int]) -> None:
        x0, y0, x1, y1 = box
        draw.rectangle(box, fill=(255, 255, 255, 255), outline=(55, 55, 55, 255), width=3)
        px0, py0, px1, py1 = x0 + 100, y0 + 36, x1 - 34, y1 - 58
        draw.line((px0, py1, px1, py1), fill=INK, width=3)
        draw.line((px0, py0, px0, py1), fill=INK, width=3)
        for t in range(11):
            x = px0 + t * (px1 - px0) / 10
            draw.line((x, py1, x, py1 + 12), fill=INK, width=2)
            draw.line((x, py0, x, py0 - 8), fill=INK, width=2)
        for t in range(10):
            y = py0 + t * (py1 - py0) / 9
            draw.line((px0, y, px0 - 14, y), fill=INK, width=2)
            draw.line((px1, y, px1 + 10, y), fill=INK, width=2)
        y_one = py0 + 0.47 * (py1 - py0)
        draw.line((px0, y_one, px1, y_one), fill=(35, 35, 35, 255), width=3)
        band = [
            (px0 + 280, y_one - 74), (px0 + 410, y_one - 94), (px0 + 560, y_one - 70),
            (px0 + 690, y_one - 46), (px0 + 690, y_one + 52), (px0 + 560, y_one + 76),
            (px0 + 410, y_one + 52), (px0 + 280, y_one + 78),
        ]
        draw.polygon(band, fill=(180, 184, 190, 100))
        pts = [
            (0.30, 0.44), (0.34, 0.49), (0.38, 0.51), (0.42, 0.48), (0.46, 0.52),
            (0.50, 0.50), (0.54, 0.46), (0.60, 0.50), (0.68, 0.58), (0.78, 0.49),
            (0.86, 0.69), (0.95, 0.50),
        ]
        for fx, fy in pts:
            x = px0 + fx * (px1 - px0)
            y = py1 - fy * (py1 - py0)
            draw.line((x, y - 26, x, y + 26), fill=(70, 70, 70, 180), width=2)
            draw.ellipse((x - 8, y - 8, x + 8, y + 8), fill=(36, 36, 36, 245))
        draw.rectangle((px0 - 30, y_one - 28, px0 + 6, y_one + 28), fill=(55, 55, 55, 255))
        draw.rectangle((px1 - 20, y_one - 28, px1 + 16, y_one + 28), fill=(95, 95, 95, 255))
        draw.text((x0 + 20, y0 + 34), "direct γ R_AA", font=font(TIMES_BOLD, 35), fill=INK)
        draw.text((x0 + 140, y0 + 46), "√sNN=200 GeV\nAu+Au, 0-92%", font=font(TIMES_BOLD, 29), fill=INK, spacing=3)
        draw.text((x0 + 615, y1 - 42), "p_T (GeV/c)", font=font(TIMES_BOLD, 31), fill=INK)

    def draw_exact_or_fallback_raa(box: tuple[int, int, int, int]) -> None:
        candidates = [
            ASSET_DIR / "direct_gamma_raa_user_constructed_prl109_fig3_backup_slide14.png",
            ASSET_DIR / "phenix_direct_gamma_raa.png",
            ASSET_DIR / "direct_gamma_raa_attached.png",
            ASSET_DIR / "direct_photon_raa.png",
        ]
        for candidate in candidates:
            if candidate.exists():
                img = open_rgba(candidate)
                paste_fit(base, img, box)
                return
        draw_raa_fallback(box)

    def draw_subset_label(
        box: tuple[int, int, int, int],
        text: str,
        color: tuple[int, int, int],
        fill: tuple[int, int, int],
        *,
        align: str = "left",
        y_offset: int = -24,
        size: int = 30,
    ) -> None:
        x0, y0, x1, _ = box
        label_font = font(TIMES_BOLD, size)
        tw, _ = text_box(draw, text, label_font)
        if align == "center":
            lx0 = round((x0 + x1 - tw) / 2) - 22
        else:
            lx0 = x0 + 24
        label_box = (lx0, y0 + y_offset, lx0 + 44 + tw, y0 + y_offset + 54)
        draw.rounded_rectangle(label_box, radius=20, fill=(*fill, 255), outline=(*color, 210), width=2)
        draw.text((label_box[0] + 18, label_box[1] + 7), text, font=label_font, fill=color)

    draw.rounded_rectangle((132, 320, 2390, 1278), radius=10, fill=(*PANEL, 255), outline=(*PANEL_EDGE, 255), width=2)

    draw.text((178, 360), "Prompt photons: production and color-neutral behavior", font=font(TIMES_BOLD, 48), fill=INK)
    draw.text((180, 418), "Direct and fragmentation photons are prompt; decay photons are backgrounds. Isolation selects the quiet prompt-photon subset.", font=font(TIMES_ITALIC, 29), fill=BLUE)

    # Production hierarchy: explicit subset contours instead of a loose bracket diagram.
    prod_card = (180, 496, 1410, 1216)
    draw.rounded_rectangle(prod_card, radius=8, fill=(255, 255, 255, 248), outline=(222, 229, 236, 255), width=2)
    draw.text((216, 528), "Production channels", font=font(TIMES_BOLD, 42), fill=INK)
    draw.text((218, 582), "Prompt photons include direct and fragmentation; decay photons are backgrounds.", font=font(TIMES_ITALIC, 26), fill=MUTED)

    prompt_box = (220, 674, 1368, 1060)
    direct_box = (252, 752, 844, 994)
    frag_box = (908, 752, 1336, 994)
    draw.rounded_rectangle(prompt_box, radius=30, fill=(255, 247, 244, 120), outline=(197, 64, 48, 210), width=3)
    draw_subset_label(
        prompt_box,
        "Prompt photon production",
        (197, 64, 48),
        (255, 247, 244),
        align="center",
        y_offset=-37,
        size=34,
    )
    draw.rounded_rectangle(direct_box, radius=20, fill=(237, 247, 254, 230), outline=(*SPHENIX_BLUE, 220), width=3)
    draw.rounded_rectangle(frag_box, radius=20, fill=(236, 247, 243, 230), outline=(*TEAL, 220), width=3)

    direct_label_box = (388, 700, 706, 748)
    frag_label_box = (930, 700, 1308, 748)
    draw.rounded_rectangle(direct_label_box, radius=19, fill=(237, 247, 254, 255), outline=(*SPHENIX_BLUE, 190), width=2)
    direct_label_font = font(TIMES_BOLD, 32)
    direct_tw, direct_th = text_box(draw, "Direct photons", direct_label_font)
    draw.text(
        (
            (direct_label_box[0] + direct_label_box[2] - direct_tw) / 2,
            (direct_label_box[1] + direct_label_box[3] - direct_th) / 2 - 3,
        ),
        "Direct photons",
        font=direct_label_font,
        fill=SPHENIX_BLUE,
    )
    draw.rounded_rectangle(frag_label_box, radius=19, fill=(236, 247, 243, 255), outline=(*TEAL, 190), width=2)
    frag_label_font = font(TIMES_BOLD, 32)
    frag_tw, frag_th = text_box(draw, "Fragmentation photons", frag_label_font)
    draw.text(
        (
            (frag_label_box[0] + frag_label_box[2] - frag_tw) / 2,
            (frag_label_box[1] + frag_label_box[3] - frag_th) / 2 - 3,
        ),
        "Fragmentation photons",
        font=frag_label_font,
        fill=TEAL,
    )

    centers = [(406, 880), (678, 880), (1118, 880)]
    labels = [("Compton scattering", INK), ("Annihilation", INK), ("Fragmentation radiation", TEAL)]
    for (cx, cy), (label_text, color), mode in zip(centers, labels, ("compton", "annihilation", "fragmentation")):
        draw_clean_feynman(cx, cy - 18, 0.64, mode)
        channel_font = font(TIMES_BOLD, 32)
        tw, _ = text_box(draw, label_text, channel_font)
        draw.text((cx - tw / 2, 1012), label_text, font=channel_font, fill=color)

    strip = (234, 1084, 1386, 1174)
    draw.rounded_rectangle(strip, radius=8, fill=(255, 251, 239, 255), outline=(238, 224, 190, 255), width=2)
    draw.text((260, 1100), "prompt photons", font=font(TIMES_BOLD, 34), fill=(197, 64, 48))
    draw.text((506, 1100), "=", font=font(TIMES_BOLD, 34), fill=INK)
    draw.text((548, 1100), "direct photons", font=font(TIMES_BOLD, 34), fill=SPHENIX_BLUE)
    draw.text((806, 1100), "+", font=font(TIMES_BOLD, 34), fill=INK)
    draw.text((852, 1100), "fragmentation photons", font=font(TIMES_BOLD, 34), fill=TEAL)
    draw.text((260, 1146), "decay photons are backgrounds", font=font(TIMES_ITALIC, 24), fill=MUTED)

    # Exact plot slot. The generator pastes a provided PNG if it exists.
    proof_card = (1458, 496, 2328, 968)
    draw.rounded_rectangle(proof_card, radius=8, fill=(255, 255, 255, 248), outline=(222, 229, 236, 255), width=2)
    draw.text((1496, 528), "Color-neutral behavior", font=font(TIMES_BOLD, 36), fill=INK)
    draw.text((1498, 576), "A RHIC direct-photon reference anchors the intuition.", font=font(TIMES_ITALIC, 26), fill=MUTED)
    plot_box = (1488, 620, 2298, 842)
    draw_exact_or_fallback_raa(plot_box)
    draw.rounded_rectangle((1496, 870, 2298, 930), radius=8, fill=(255, 251, 239, 255), outline=(238, 224, 190, 255), width=2)
    draw.rectangle((1496, 870, 1508, 930), fill=(*PHOTON, 255))
    draw.text((1534, 887), "Direct photons stay near unity in Au+Au.", font=font(TIMES_BOLD, 29), fill=INK)

    draw.rounded_rectangle((1458, 988, 2328, 1216), radius=8, fill=(255, 255, 255, 245), outline=(222, 229, 236, 255), width=2)
    draw.text((1496, 1016), "Isolation definition", font=font(TIMES_BOLD, 36), fill=INK)

    def draw_isolation_equation(x: int, y: int) -> None:
        eq_font = font(TIMES_BOLD, 28)
        script_font = font(TIMES_BOLD, 15)
        x_next, _ = draw_rich_text(draw, (x, y + 7), "E_T^iso =", eq_font, INK)
        sigma_x = x_next + 18
        draw.text((sigma_x, y - 2), "Σ", font=font(TIMES_BOLD, 50), fill=INK)
        sigma_w, _ = text_box(draw, "Σ", font(TIMES_BOLD, 50))
        sub = "cone towers"
        sub_w, _ = text_box(draw, sub, script_font)
        draw.text((sigma_x + sigma_w / 2 - sub_w / 2, y + 43), sub, font=script_font, fill=MUTED)
        rest_x = round(sigma_x + sigma_w + 16)
        draw_rich_text(draw, (rest_x, y + 7), "E_T^tower - E_T^candidate", eq_font, INK)

    def draw_isolation_cartoon(box: tuple[int, int, int, int], *, busy: bool) -> None:
        x0, y0, x1, y1 = box
        cx = (x0 + x1) // 2
        outline = (197, 64, 48) if busy else PHOTON_DARK
        title_color = (197, 64, 48) if busy else BLUE
        card_fill = (255, 248, 246, 255) if busy else (244, 249, 252, 255)
        cone_fill = (255, 247, 244, 178) if busy else (244, 250, 253, 208)
        label = "non-isolated" if busy else "isolated"

        draw.rounded_rectangle(box, radius=7, fill=card_fill, outline=(222, 229, 236, 255), width=2)
        title_font = font(TIMES_BOLD, 22)
        tw, _ = text_box(draw, label, title_font)
        draw.text((cx - tw / 2, y0 + 7), label, font=title_font, fill=title_color)

        rim_y = y0 + 58
        apex = (cx, y1 - 16)
        rim_w = min(58, (x1 - x0) // 2 - 16)
        rim_left = (cx - rim_w, rim_y)
        rim_right = (cx + rim_w, rim_y)
        draw.polygon([apex, rim_left, rim_right], fill=cone_fill)
        draw.line((apex, rim_left), fill=(*outline, 220), width=3)
        draw.line((apex, rim_right), fill=(*outline, 220), width=3)
        draw.arc((cx - rim_w, rim_y - 12, cx + rim_w, rim_y + 12), 0, 180, fill=(*outline, 224), width=3)
        draw.arc((cx - rim_w, rim_y - 12, cx + rim_w, rim_y + 12), 180, 360, fill=(*outline, 92), width=2)

        def activity_ray(start: tuple[int, int], end: tuple[int, int]) -> None:
            sx, sy = start
            ex, ey = end
            draw.line((sx, sy, ex, ey), fill=(255, 255, 255, 214), width=5)
            draw.line((sx, sy, ex, ey), fill=(*INK, 218), width=2)
            angle = math.atan2(ey - sy, ex - sx)
            head = 6
            left = (ex - head * math.cos(angle - 0.48), ey - head * math.sin(angle - 0.48))
            right = (ex - head * math.cos(angle + 0.48), ey - head * math.sin(angle + 0.48))
            draw.polygon([(ex, ey), left, right], fill=(*INK, 218))

        if busy:
            for start, end in (
                ((cx - 5, apex[1] - 3), (cx - 29, rim_y + 34)),
                ((cx - 1, apex[1] - 8), (cx - 12, rim_y + 58)),
                ((cx + 5, apex[1] - 7), (cx + 8, rim_y + 46)),
                ((cx + 10, apex[1] - 4), (cx + 26, rim_y + 36)),
                ((cx - 9, apex[1] - 12), (cx - 3, rim_y + 26)),
            ):
                activity_ray(start, end)

        photon_x = cx - 4
        wave = feynman_points((photon_x, rim_y - 15), (photon_x, apex[1] - 9), 2.2, 5.5, 108)
        draw_polyline(draw, wave, (*PHOTON_DARK, 228), 2)
        gamma_font = font(TIMES_ITALIC, 21)
        draw.text(
            (photon_x + 10, rim_y - 34),
            "γ",
            font=gamma_font,
            fill=INK,
            stroke_width=1,
            stroke_fill=(255, 255, 255, 235),
        )

    formula_box = (1496, 1066, 1986, 1138)
    draw.rounded_rectangle(formula_box, radius=8, fill=(255, 251, 239, 255), outline=(238, 224, 190, 255), width=2)
    draw.rectangle((1496, 1066, 1508, 1138), fill=(*PHOTON, 255))
    draw_isolation_equation(1536, 1078)
    draw_rich_text(
        draw,
        (1496, 1150),
        "small E_T^iso → low nearby activity",
        font(TIMES, 24),
        MUTED,
    )
    draw.text(
        (1496, 1186),
        "Suppresses fragmentation-rich nearby activity.",
        font=font(TIMES_BOLD, 22),
        fill=TEAL,
    )
    draw.line((1994, 1024, 1994, 1204), fill=(222, 229, 236, 255), width=2)
    draw_isolation_cartoon((2016, 1024, 2162, 1206), busy=False)
    draw_isolation_cartoon((2178, 1024, 2314, 1206), busy=True)


def write_isolated_photon_script(output_dir: Path) -> Path:
    script = """# HP2026 Slide 4 Speaker Script

Now I want to define the actual physics object of the talk. The important point is that prompt photons are photons associated with the short-distance parton scattering, not photons from neutral-meson decays.

On the left, I am showing that hierarchy explicitly. The first two diagrams are the direct photon channels: Compton scattering and annihilation. The fragmentation channel is also prompt, because the photon is still associated with the hard scattering, but experimentally it tends to come with more nearby activity. Decay photons are outside this prompt-photon box and are the background this analysis has to suppress.

The reason photons are powerful is shown on the right. Once the photon is produced, it is color neutral, so it does not undergo the same strong final-state energy loss as a colored parton. The direct-photon reference sits near unity, which is the qualitative behavior we want from a calibrated electromagnetic tag.

Isolation is the experimental step that makes this object clean enough to measure. Operationally, we calculate an isolation energy: the transverse energy from towers in a cone around the candidate, minus the candidate's own transverse energy. The two small graphics show the intuition directly: the isolated case has little activity around the photon candidate, while the non-isolated case has additional nearby particles or calorimeter energy in that same cone.

That is why isolation is especially useful for the fragmentation component. Fragmentation photons tend to come with nearby activity from the parent parton, so the isolation requirement preferentially suppresses fragmentation-rich and decay-rich candidates. So the object for this talk is not just any photon; it is an isolated prompt photon in p+p, which gives the baseline for future heavy-ion photon measurements.
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

    logo = load_sphenix_logo()
    if logo is not None:
        paste_fit(img, logo, (2188, 58, 2432, 164), anchor="right")

    draw.text((132, 98), PHOTON_MOTIVATION_TITLE, font=font(TIMES_BOLD, 90), fill=INK)
    draw.text((136, 205), PHOTON_MOTIVATION_SUBTITLE, font=font(TIMES_ITALIC, 43), fill=BLUE)
    draw.line((132, 292, W - 132, 292), fill=(221, 226, 232), width=3)

    draw_prompt_photon_integrated_slide(img)
    draw_recreated_footer(img)

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
                "color_neutrality_plot_asset": str(
                    (ASSET_DIR / "direct_gamma_raa_user_constructed_prl109_fig3_backup_slide14.png").relative_to(ROOT)
                ),
                "color_neutrality_plot_source": {
                    "paper_pdf": str((ROOT / "usefulDocs" / "PHENIX_direct_photon_RAA_PRL109_152302.pdf").relative_to(ROOT)),
                    "source_note": "User-constructed white-background screenshot crop from PRL 109.152302 Fig. 3 / backup slide 14.",
                },
            },
        ],
        "data_source_evidence": DATA_SOURCE_EVIDENCE,
        "source_basis": [
            "Official BNL RHIC sPHENIX detector imagery for Slide 2.",
            "Yeonju DIS2026 and Hanpu AUM2026 reference decks for p+p data-taking numbers.",
            "PHENIX PRL 109.152302 direct-photon R_AA plot used as a color-neutrality reference via Justin's constructed backup-slide crop.",
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
    if variant == "sphenix-subsystems-recreated":
        return render_recreated_sphenix_subsystems(output_dir)
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
            "sphenix-subsystems-recreated",
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

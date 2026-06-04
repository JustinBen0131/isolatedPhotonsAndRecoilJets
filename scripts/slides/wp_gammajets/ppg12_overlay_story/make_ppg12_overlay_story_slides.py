#!/usr/bin/env python3
"""Build full-slide PNG candidates for the pp/PPG12 BDT overlay story."""

from __future__ import annotations

import json
import math
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from PIL import Image, ImageDraw, ImageFilter, ImageFont

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ModuleNotFoundError:
    plt = None
    HAS_MATPLOTLIB = False


REPO = Path(__file__).resolve().parents[4]
OUT = REPO / "dataOutput/slides/wp_gammajets_6_3_26/ppg12_overlay_story_20260603"
RAW = OUT / "raw_assets"
SUMMARY = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_truthWindowOverlayFix_20260528_1305"
    / "validation/fullsim_shuhang_overlay_truthwindow_ppg12_weighted"
    / "pp_currentian_basev3e_bdt_score_overlay_vs_ppg12_pp_noCent_bdt_22_28_noNPB_eta0-pt3-cut0_summary.json"
)

W, H = 2560, 1440
FONT_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES = FONT_DIR / "Times New Roman.ttf"
TIMES_BOLD = FONT_DIR / "Times New Roman Bold.ttf"
TIMES_ITALIC = FONT_DIR / "Times New Roman Italic.ttf"

INK = (18, 24, 38)
MUTED = (76, 83, 101)
BLUE = (37, 99, 235)
RED = (220, 38, 38)
GREEN = (5, 122, 85)
AMBER = (180, 83, 9)
LIGHT_BLUE = (238, 246, 255)
LIGHT_GREEN = (234, 251, 242)
LIGHT_AMBER = (255, 247, 237)
LIGHT_RED = (255, 241, 242)
LINE = (211, 220, 230)


def font(size: int, *, bold: bool = False, italic: bool = False) -> ImageFont.FreeTypeFont:
    if bold:
        return ImageFont.truetype(str(TIMES_BOLD), size)
    if italic:
        return ImageFont.truetype(str(TIMES_ITALIC), size)
    return ImageFont.truetype(str(TIMES), size)


F = {
    "title": font(76, bold=True),
    "subtitle": font(37),
    "kicker": font(36, bold=True),
    "body": font(38),
    "body_small": font(34),
    "body_bold": font(38, bold=True),
    "label": font(42, bold=True),
    "tiny": font(31),
    "tiny_bold": font(31, bold=True),
    "panel": font(34, bold=True),
}


def canvas() -> Image.Image:
    return Image.new("RGB", (W, H), "white")


def rounded(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    *,
    fill: tuple[int, int, int] = (255, 255, 255),
    outline: tuple[int, int, int] = LINE,
    radius: int = 22,
    width: int = 3,
) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def wrap_text(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.FreeTypeFont, max_w: int) -> list[str]:
    lines: list[str] = []
    for paragraph in text.split("\n"):
        words = paragraph.split()
        if not words:
            lines.append("")
            continue
        line = words[0]
        for word in words[1:]:
            test = f"{line} {word}"
            if draw.textlength(test, font=fnt) <= max_w:
                line = test
            else:
                lines.append(line)
                line = word
        lines.append(line)
    return lines


def draw_textbox(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    fnt: ImageFont.FreeTypeFont,
    *,
    max_w: int,
    fill: tuple[int, int, int] = INK,
    line_gap: int = 8,
) -> int:
    x, y = xy
    for line in wrap_text(draw, text, fnt, max_w):
        draw.text((x, y), line, font=fnt, fill=fill)
        y += fnt.size + line_gap
    return y


def draw_title(draw: ImageDraw.ImageDraw, title: str, subtitle: str | None = None) -> None:
    draw.text((105, 52), title, font=F["title"], fill=INK)
    if subtitle:
        draw.text((109, 140), subtitle, font=F["subtitle"], fill=MUTED)


def shadow_paste(base: Image.Image, img: Image.Image, box: tuple[int, int, int, int], *, radius: int = 18) -> None:
    x1, y1, x2, y2 = box
    bw, bh = x2 - x1, y2 - y1
    # Approximate Google Slides shadow baseline: opacity 18%, angle 60,
    # distance 4 px, blur radius 14 px.
    panel = Image.new("RGBA", (bw + 36, bh + 36), (0, 0, 0, 0))
    sh = Image.new("RGBA", (bw, bh), (0, 0, 0, 46))
    sh = sh.filter(ImageFilter.GaussianBlur(14))
    panel.alpha_composite(sh, (2, 3))
    resized = fit_image(img, (bw, bh), fill=(255, 255, 255), mode="contain")
    panel.alpha_composite(resized.convert("RGBA"), (0, 0))
    base.alpha_composite(panel, (x1, y1))
    draw = ImageDraw.Draw(base)
    draw.rounded_rectangle((x1, y1, x2, y2), radius=radius, outline=LINE, width=2)


def fit_image(
    img: Image.Image,
    size: tuple[int, int],
    *,
    fill: tuple[int, int, int] = (255, 255, 255),
    mode: str = "contain",
) -> Image.Image:
    tw, th = size
    iw, ih = img.size
    scale = max(tw / iw, th / ih) if mode == "cover" else min(tw / iw, th / ih)
    nw, nh = max(1, int(iw * scale)), max(1, int(ih * scale))
    resample = Image.Resampling.LANCZOS
    resized = img.resize((nw, nh), resample)
    if mode == "cover":
        left = max(0, (nw - tw) // 2)
        top = max(0, (nh - th) // 2)
        return resized.crop((left, top, left + tw, top + th))
    out = Image.new("RGBA", (tw, th), fill + (255,))
    out.alpha_composite(resized.convert("RGBA"), ((tw - nw) // 2, (th - nh) // 2))
    return out


def padded_image(
    img: Image.Image,
    *,
    left: int = 18,
    top: int = 14,
    right: int = 18,
    bottom: int = 22,
    fill: tuple[int, int, int] = (255, 255, 255),
) -> Image.Image:
    out = Image.new("RGB", (img.width + left + right, img.height + top + bottom), fill)
    out.paste(img.convert("RGB"), (left, top))
    return out


def label_bar(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], text: str, fill: tuple[int, int, int]) -> None:
    rounded(draw, box, fill=fill, outline=fill, radius=14, width=1)
    tw = draw.textlength(text, font=F["panel"])
    draw.text(((box[0] + box[2] - tw) / 2, box[1] + 11), text, font=F["panel"], fill=INK)


def crop_old_slide42_reference() -> Image.Image:
    img = Image.open(RAW / "old_deck_slide42_stale.png").convert("RGB")
    return img.crop((130, 145, 815, 875))


def crop_old_slide42_rhs() -> Image.Image:
    img = Image.open(RAW / "old_deck_slide42_stale.png").convert("RGB")
    # Plot-only crop from the old slide-42 RHS. Do not include the old
    # blue slide header; slide-level labels are redrawn by this builder.
    return img.crop((845, 183, 1420, 875))


def clean_slide_thumb(name: str) -> Image.Image:
    img = Image.open(RAW / name).convert("RGB")
    draw = ImageDraw.Draw(img)
    draw.rectangle((1500, 825, 1600, 900), fill="white")
    return img


def make_fixed_overlay_with_ratio() -> Path:
    out = OUT / "fixed_overlay_contract_with_ratio.png"
    if not HAS_MATPLOTLIB:
        if not out.exists():
            raise RuntimeError("matplotlib is required to build the fixed overlay plot from scratch")
        return out
    data = json.loads(SUMMARY.read_text())
    bins = np.array(data["bins"], dtype=float)
    x = bins[:-1]
    y_sig = np.array(data["this_analysis_signal_hist"], dtype=float)
    y_inc = np.array(data["this_analysis_inclusive_hist"], dtype=float)
    y_ref_sig = np.array(data["shuhang_signal_hist"], dtype=float)
    y_ref_inc = np.array(data["shuhang_inclusive_hist"], dtype=float)
    diff = y_ref_inc - y_inc

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.25,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig = plt.figure(figsize=(7.3, 9.0), dpi=240)
    gs = fig.add_gridspec(2, 1, height_ratios=[4.2, 1.15], hspace=0.03)
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)
    kwargs = dict(where="post", linewidth=2.2)
    ax.step(x, y_sig, color="#dc2626", label="This analysis signal", **kwargs)
    ax.step(x, y_inc, color="#2563eb", label="This analysis inclusive", **kwargs)
    ax.step(x, y_ref_sig, color="#991b1b", linestyle="--", label="PPG12/Shuhang signal", **kwargs)
    ax.step(x, y_ref_inc, color="#1e40af", linestyle="--", label="PPG12/Shuhang inclusive", **kwargs)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, max(0.22, float(max(y_sig.max(), y_inc.max(), y_ref_sig.max(), y_ref_inc.max())) * 1.18))
    ax.set_ylabel("unit-normalized counts", fontsize=26)
    ax.text(
        0.055,
        0.82,
        r"$\it{\bf{sPHENIX}}$ Internal" + "\n"
        + r"$p$+$p$ $\sqrt{s}=200$ GeV" + "\n"
        + r"$|\eta|<0.7$" + "\n"
        + r"$22<E_T<28$ GeV, w/o NPB cut",
        transform=ax.transAxes,
        fontsize=19,
        va="top",
    )
    ax.legend(frameon=False, fontsize=17, loc="upper right", handlelength=2.6)
    ax.tick_params(labelbottom=False, labelsize=20, length=8)
    rax.scatter(x + 0.01, diff, s=12, color="black")
    rax.axhline(0, color="black", linestyle="--", linewidth=1)
    rax.set_ylim(-0.06, 0.06)
    rax.set_ylabel("PPG12 -\nthis inc.", fontsize=22)
    rax.set_xlabel("BDT score", fontsize=30)
    rax.tick_params(labelsize=20, length=8)
    fig.subplots_adjust(left=0.18, right=0.98, top=0.98, bottom=0.11)
    fig.savefig(out)
    plt.close(fig)
    return out


def bullet_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    title: str,
    body: str,
    *,
    fill: tuple[int, int, int],
    accent: tuple[int, int, int],
    title_font: ImageFont.FreeTypeFont = F["label"],
    body_font: ImageFont.FreeTypeFont = F["body_small"],
) -> None:
    rounded(draw, box, fill=fill, outline=accent, radius=24, width=3)
    x, y = box[0] + 30, box[1] + 26
    draw.text((x, y), title, font=title_font, fill=accent)
    draw_textbox(draw, (x, y + title_font.size + 16), body, body_font, max_w=box[2] - box[0] - 60, fill=INK, line_gap=8)


def compact_key_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    title: str,
    rows: list[tuple[str, str]],
    *,
    fill: tuple[int, int, int],
    accent: tuple[int, int, int],
) -> None:
    rounded(draw, box, fill=fill, outline=accent, radius=24, width=3)
    x, y = box[0] + 30, box[1] + 22
    draw.text((x, y), title, font=F["label"], fill=accent)
    y += F["label"].size + 18
    max_w = box[2] - box[0] - 60
    for label, body in rows:
        draw.text((x, y), label, font=F["tiny_bold"], fill=INK)
        y = draw_textbox(
            draw,
            (x, y + F["tiny_bold"].size + 3),
            body,
            F["tiny"],
            max_w=max_w,
            fill=INK,
            line_gap=5,
        )
        y += 12


def labeled_line_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    title: str,
    rows: list[tuple[str, str]],
    *,
    fill: tuple[int, int, int],
    accent: tuple[int, int, int],
) -> None:
    rounded(draw, box, fill=fill, outline=accent, radius=24, width=3)
    x, y = box[0] + 30, box[1] + 24
    draw.text((x, y), title, font=F["label"], fill=accent)
    y += F["label"].size + 22
    for label, body in rows:
        draw.text((x, y), label, font=F["tiny_bold"], fill=INK)
        body_x = x + int(draw.textlength(label, font=F["tiny_bold"])) + 10
        draw.text((body_x, y), body, font=F["tiny"], fill=INK)
        y += F["tiny"].size + 18


def diagnostic_definition_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    title: str,
    rows: list[tuple[str, str]],
    *,
    fill: tuple[int, int, int],
    accent: tuple[int, int, int],
    note: str | None = None,
) -> None:
    rounded(draw, box, fill=fill, outline=accent, radius=24, width=3)
    x, y = box[0] + 30, box[1] + 22
    draw.text((x, y), title, font=F["label"], fill=accent)
    y += F["label"].size + 18
    row_label_font = F["body_bold"]
    row_body_font = F["body_small"]
    label_w = max(int(draw.textlength(label, font=row_label_font)) for label, _ in rows) + 22
    body_x = x + label_w
    max_w = box[2] - body_x - 30
    for label, body in rows:
        draw.text((x, y), label, font=row_label_font, fill=INK)
        y_next = draw_textbox(draw, (body_x, y + 2), body, row_body_font, max_w=max_w, fill=INK, line_gap=5)
        y = max(y + row_body_font.size + 16, y_next + 12)
    if note:
        draw.line((x, box[3] - 70, box[2] - 30, box[3] - 70), fill=LINE, width=3)
        draw_textbox(draw, (x, box[3] - 52), note, F["body_small"], max_w=box[2] - box[0] - 60, fill=MUTED, line_gap=5)


def slide1_reference_object() -> Path:
    out = OUT / "slide01_define_ppg12_reference_object.png"
    im = canvas().convert("RGBA")
    draw = ImageDraw.Draw(im)
    draw_title(
        draw,
        "Define the PPG12 reference before comparing to it",
        "From the current PPG12 IAN Fig. 13: BDT score, 22 < ET < 28 GeV, no NPB cut.",
    )

    ref = crop_old_slide42_reference()
    shadow_paste(im, ref, (115, 265, 1045, 1265), radius=18)
    label_bar(draw, (140, 215, 1020, 276), "PPG12 IAN Fig. 13: 22 < ET < 28 GeV", (239, 242, 246))

    card_x, card_r = 1120, 2415
    bullet_card(
        draw,
        (card_x, 250, card_r, 435),
        "Reference object",
        "Unit-normalized BDT-score shapes for data, signal MC, inclusive MC, and NPB-tagged data.",
        fill=LIGHT_BLUE,
        accent=BLUE,
        title_font=F["label"],
        body_font=F["body_small"],
    )
    bullet_card(
        draw,
        (card_x, 470, card_r, 645),
        "Comparison question",
        "Signal-like BDT shape compared with the generic candidate population from inclusive jet MC.",
        fill=LIGHT_GREEN,
        accent=GREEN,
        title_font=F["label"],
        body_font=F["body_small"],
    )
    bullet_card(
        draw,
        (card_x, 680, card_r, 855),
        "Overlay must match the fill contract",
        "Raw-inclusive rows, PPG12 weights, truth-jet stitch windows, 22 < ET < 28 GeV, no NPB cut, unit normalization.",
        fill=LIGHT_AMBER,
        accent=AMBER,
        title_font=F["label"],
        body_font=F["body_small"],
    )

    rounded(draw, (card_x, 900, card_r, 1315), fill=(250, 250, 250), outline=LINE, radius=24)
    draw.text((1160, 935), "Key point: what the two MC curves mean", font=F["label"], fill=INK)
    draw.line((1160, 1000, 2375, 1000), fill=LINE, width=3)

    label_x = 1160
    body_x = 1370
    row_gap = 92
    y = 1040
    draw.text((label_x, y), "Signal =", font=F["body_bold"], fill=INK)
    draw.text((body_x, y), "truth prompt photons.", font=F["body_small"], fill=INK)
    y += row_gap
    draw.text((label_x, y), "Inclusive =", font=F["body_bold"], fill=INK)
    draw_textbox(
        draw,
        (body_x, y),
        "reconstructed candidates in inclusive jet MC after analysis selections.",
        F["body_small"],
        max_w=990,
        fill=INK,
        line_gap=6,
    )
    y += 126
    draw.line((1160, y - 34, 2375, y - 34), fill=LINE, width=3)
    draw.text((label_x, y), "This asks:", font=F["body_bold"], fill=INK)
    draw_textbox(
        draw,
        (body_x, y),
        "Signal-like output versus generic jet/inclusive candidates.",
        F["body_small"],
        max_w=990,
        fill=INK,
        line_gap=6,
    )
    im.convert("RGB").save(out, quality=95)
    return out


def slide2_overlay_exercise(fixed_path: Path) -> Path:
    out = OUT / "slide02_overlay_exercise_three_way_comparison.png"
    im = canvas().convert("RGBA")
    draw = ImageDraw.Draw(im)
    draw_title(
        draw,
        "The overlay exercise separated two valid BDT diagnostics",
        "The corrected inclusive diagnostic now closely matches PPG12 Fig. 13, with only minor residual inconsistency.",
    )

    panels = [
        ((70, 245, 760, 918), "Truth-labeled diagnostic", padded_image(crop_old_slide42_rhs(), bottom=34), LIGHT_RED, RED),
        ((935, 245, 1625, 918), "Fixed inclusive diagnostic", padded_image(Image.open(fixed_path).convert("RGB"), bottom=28), LIGHT_GREEN, GREEN),
        ((1800, 245, 2490, 918), "PPG12 IAN reference", padded_image(crop_old_slide42_reference(), bottom=30), (239, 242, 246), INK),
    ]
    for box, title, img, fill, accent in panels:
        label_bar(draw, (box[0], box[1] - 66, box[2], box[1] - 8), title, fill)
        shadow_paste(im, img, box, radius=18)
        draw.rounded_rectangle(box, radius=18, outline=accent, width=5)

    for x1, x2 in ((800, 895), (1665, 1760)):
        draw.line((x1, 582, x2, 582), fill=(115, 123, 140), width=6)
        draw.polygon([(x2, 582), (x2 - 34, 558), (x2 - 34, 606)], fill=(115, 123, 140))

    diagnostic_definition_card(
        draw,
        (120, 955, 1240, 1308),
        "Truth-labeled training diagnostic",
        [
            ("Signal =", "truth prompt photons."),
            ("Background =", "truth-labeled non-signal candidates."),
            ("Question =", "Can the model separate its labeled classes?"),
        ],
        fill=LIGHT_RED,
        accent=RED,
        note="Valid diagnostic; not the PPG12 Fig. 13 object.",
    )
    diagnostic_definition_card(
        draw,
        (1345, 955, 2465, 1308),
        "Physics/template-style inclusive diagnostic",
        [
            ("Signal =", "truth prompt photons."),
            ("Inclusive =", "inclusive-jet MC candidates after selections."),
            ("Question =", "signal-like BDT shape vs generic inclusive candidates."),
            ("Result =", "green overlay closely matches PPG12 Fig. 13."),
        ],
        fill=LIGHT_GREEN,
        accent=GREEN,
    )
    rounded(draw, (120, 1328, 2465, 1388), fill=(250, 250, 250), outline=LINE, radius=20, width=2)
    key_label = "Key comparison:"
    key_body = " red panel = truth-labeled background; green/reference panels = inclusive jet MC candidates after selections."
    draw.text((155, 1342), key_label, font=F["body_bold"], fill=INK)
    draw.text((155 + int(draw.textlength(key_label, font=F["body_bold"])) + 8, 1342), key_body, font=F["body"], fill=INK)
    im.convert("RGB").save(out, quality=95)
    return out


def slide3_direct_model_confirmation() -> Path:
    out = OUT / "slide03_direct_model_confirmation.png"
    im = canvas().convert("RGBA")
    draw = ImageDraw.Draw(im)
    draw_title(
        draw,
        "Same truth-labeled diagnostic, now using the PPG12-trained BDT",
        "This is the previous slide's red diagnostic: same truth-tagged signal/background rows, same plotting pipeline, PPG12 model applied out-of-box.",
    )

    s41 = clean_slide_thumb("old_deck_slide41.png")
    # Crop out old title/subtitle so the visual acts as evidence inside the new slide.
    crop = s41.crop((20, 108, 1048, 862))
    shadow_paste(im, crop, (115, 265, 1575, 1245), radius=18)
    label_bar(draw, (140, 218, 1490, 276), "Same rows + same plotting pipeline: two independently trained models", LIGHT_BLUE)

    bullet_card(
        draw,
        (1655, 245, 2410, 500),
        "Connection to previous slide",
        "This is the truth-labeled training diagnostic: signal is truth prompt photons; background is truth-labeled non-signal candidates.",
        fill=LIGHT_BLUE,
        accent=BLUE,
        title_font=F["label"],
        body_font=F["body_small"],
    )
    bullet_card(
        draw,
        (1655, 545, 2410, 800),
        "What changes here",
        "The rows and plotting code stay fixed. The PPG12-trained BDT is applied out-of-box to test whether it gives comparable separation.",
        fill=LIGHT_GREEN,
        accent=GREEN,
        title_font=F["label"],
        body_font=F["body_small"],
    )
    bullet_card(
        draw,
        (1655, 845, 2410, 1100),
        "Why this matters",
        "The shapes agree closely, so the model-comparison path is consistent. The slide-42 issue was the Fig. 13 inclusive histogram fill contract.",
        fill=LIGHT_AMBER,
        accent=AMBER,
        title_font=F["label"],
        body_font=F["body_small"],
    )
    im.convert("RGB").save(out, quality=95)
    return out


def slide4_same_row_summary() -> Path:
    out = OUT / "slide04_same_row_summary_and_backburner.png"
    im = canvas().convert("RGBA")
    draw = ImageDraw.Draw(im)
    draw_title(
        draw,
        "Same-row audits quantify the model-level agreement",
        "With common candidate rows and plotting conditions, the independent models give similar discrimination.",
    )

    s39 = clean_slide_thumb("old_deck_slide39.png")
    s40 = clean_slide_thumb("old_deck_slide40.png")
    p1 = padded_image(s39.crop((30, 116, 1032, 835)), left=10, top=10, right=10, bottom=16)
    p2 = padded_image(s40.crop((64, 145, 970, 838)), left=10, top=10, right=10, bottom=16)
    label_bar(draw, (120, 218, 1160, 276), "Score correlation: same candidates", LIGHT_BLUE)
    label_bar(draw, (1400, 218, 2440, 276), "ROC closure: similar discrimination", LIGHT_GREEN)
    shadow_paste(im, p1, (95, 305, 1215, 1015), radius=18)
    shadow_paste(im, p2, (1345, 305, 2465, 1015), radius=18)

    bullet_card(
        draw,
        (95, 1082, 790, 1335),
        "Score-by-score check",
        "Same candidate scored by both models. Pearson r = 0.976, so the score ordering is strongly aligned.",
        fill=LIGHT_BLUE,
        accent=BLUE,
        title_font=F["label"],
        body_font=F["body_small"],
    )
    bullet_card(
        draw,
        (930, 1082, 1625, 1335),
        "ROC check",
        "Discrimination power is similar: this-analysis AUC = 0.904; PPG12 model AUC = 0.913.",
        fill=LIGHT_GREEN,
        accent=GREEN,
        title_font=F["label"],
        body_font=F["body_small"],
    )
    bullet_card(
        draw,
        (1765, 1082, 2465, 1335),
        "What remains",
        "Perfect pp closure can be audited later: score calibration, sample composition, training details, and histogram normalization.",
        fill=LIGHT_AMBER,
        accent=AMBER,
        title_font=F["label"],
        body_font=F["body_small"],
    )
    im.convert("RGB").save(out, quality=95)
    return out


SCRIPTS = {
    "slide01_define_ppg12_reference_object.md": """# WP GammaJets Slide 1 Script - Define the PPG12 reference object

To start this section, I want to define exactly what comparison object I am using from PPG12.

The figure on the left is the BDT-score panel from the current PPG12 IAN Fig. 13 in the 22 to 28 GeV cluster ET bin, with no NPB cut. It is not just a model output shown in isolation. It is a unit-normalized histogram comparison between data, signal MC, inclusive MC, and the NPB-tagged data template, with the lower panel showing data minus inclusive MC.

The key point is that the signal curve is truth prompt photons, while the inclusive curve is whatever reconstructed candidates appear in inclusive jet MC after analysis selections. So this asks what the signal-like BDT distribution looks like compared to the generic candidate population expected from jet and inclusive samples.

To compare my output to this PPG12 figure, I need to reproduce the same fill contract: the same pT and eta bin, the no-NPB-cut convention, the same normalization logic, the same inclusive sample weighting, and the same stitch-window rules.
""",
    "slide02_overlay_exercise_three_way_comparison.md": """# WP GammaJets Slide 2 Script - What the PPG12 overlay exercise taught us

This slide shows the core lesson from the PPG12 overlay exercise.

On the left is the old slide-42-style comparison. The problem was that I was comparing to the PPG12 output histogram, but my inclusive curve had not yet been built with the same PPG12 fill contract. That made the comparison look more discrepant than a true object-to-object comparison should.

The middle panel is the regenerated comparison. Here the inclusive side uses raw-inclusive rows, PPG12 shower-shape sample weights, and the explicit PPG12 truth-jet stitch-window selection before the final unit normalization. That is the same mechanism we needed to emulate to compare to the IAN object.

The right panel is the original PPG12 Fig. 13 reference. The important point is that the fixed version now compares the right kind of object to the PPG12 object. There are still small residual differences, but this no longer looks like evidence that the model implementation or the Au+Au training is fundamentally broken.
""",
    "slide03_direct_model_confirmation.md": """# WP GammaJets Slide 3 Script - Direct model confirmation

This slide is the same-row confirmation check.

Here I am not comparing my output to the pre-filled PPG12 Fig. 13 histogram. I am taking the same reconstructed candidate rows and scoring them two ways: once with my current analysis BDT model and once with the PPG12-trained BDT model.

That isolates the model-consistency question. These are two independently trained photon-ID BDT models with the same intended discrimination target, run through the same plotting mechanism. Under that common-row mechanism, the score shapes are close. So the old slide 42 mismatch should not be interpreted as evidence that the BDT implementation, or the Au+Au training, was fundamentally broken.

The conclusion is that we needed two separate diagnostics. Slide 41 confirms model/scorer consistency. The corrected Fig. 13-style overlay tests whether we reproduced the PPG12 production histogram contract.
""",
    "slide04_same_row_summary_and_backburner.md": """# WP GammaJets Slide 4 Script - Same-row summary and back-burner follow-up

This final slide closes the loop on whether the pp baseline is consistent enough to move on.

The left panel is the same-row score correlation. Every point is the same candidate scored by both model paths. The right panel is the same-row ROC comparison. Together they show that the two model implementations are very close: the Pearson correlation is about 0.976, and the AUCs are 0.904 and 0.913.

That does not mean every PPG12 overlay is pixel-perfect. The remaining differences can still come from score calibration, sample composition, the exact TMVA versus XGBoost implementation details, or final histogram normalization choices. But those are now a back-burner closure audit, not the blocker for this working point.

The conclusion I would present is that the discrepancy from last week is understood at the level that matters for the meeting: the model comparison is consistent, the PPG12 histogram comparison needed the right fill contract, and the corrected overlay is consistent enough to move forward.
""",
}


def write_scripts() -> list[Path]:
    paths = []
    for name, text in SCRIPTS.items():
        path = OUT / name
        path.write_text(text.strip() + "\n")
        paths.append(path)
    return paths


def make_contact_sheet(slides: list[Path]) -> Path:
    out = OUT / "ppg12_overlay_story_contact_sheet.png"
    thumb_w, thumb_h = 640, 360
    sheet = Image.new("RGB", (thumb_w * 2 + 60, thumb_h * 2 + 110), "white")
    draw = ImageDraw.Draw(sheet)
    for i, path in enumerate(slides):
        img = Image.open(path).convert("RGB").resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        x = 20 + (i % 2) * (thumb_w + 20)
        y = 55 + (i // 2) * (thumb_h + 35)
        sheet.paste(img, (x, y))
        draw.text((x, y - 31), f"Slide {i + 1}: {path.stem}", font=font(20, bold=True), fill=INK)
        draw.rectangle((x, y, x + thumb_w, y + thumb_h), outline=LINE, width=2)
    sheet.save(out, quality=95)
    return out


def write_manifest(slides: list[Path], scripts: list[Path], contact: Path, fixed: Path) -> Path:
    manifest = {
        "slide_family": "WP GammaJets 6/3/26 PPG12 overlay story",
        "output_dir": str(OUT),
        "slides": [str(p) for p in slides],
        "scripts": [str(p) for p in scripts],
        "contact_sheet": str(contact),
        "generated_fixed_overlay_with_ratio": str(fixed),
        "source_assets": {
            "old_deck": {
                "presentation_id": "1GAoEcN9UGOkUVxTgT0m1vRhghl3O1jd28LwSHg9klBs",
                "title": "WP_GammaJets_6_1_26",
                "slides": {
                    "39": "g3e4e4b07144_0_84",
                    "40": "g3e4e4b07144_0_29",
                    "41": "g3e4e4b07144_0_235",
                    "42": "g3e8216e3697_0_0",
                },
            },
            "fixed_overlay_summary": str(SUMMARY),
            "ppg12_ian": "usefulDocs/PPG12_analysis_note_2026-05-21_v4_current_IAN.pdf, Fig. 13 BDT-score crop",
        },
        "key_metrics": {
            "same_row_audit": {"pearson_r": 0.976, "this_analysis_auc": 0.904, "ppg12_split_tmva_auc": 0.913},
            "truth_window_overlay": {
                "signal_rows_after_cuts": 4042616,
                "inclusive_rows_after_cuts": 926228,
                "this_inclusive_frac_score_lt_0p1": 0.2695,
                "ppg12_inclusive_frac_score_lt_0p1": 0.2334,
            },
        },
        "interpretation": [
            "Old slide-42-style comparison is stale for PPG12 histogram-equivalence claims.",
            "Corrected overlay uses raw-inclusive rows, PPG12 shower-shape weights, and truth-jet stitch-window selection.",
            "Same-row model audits remain evidence that the pp BDT implementation is broadly aligned.",
            "Residual perfect-overlay closure is a lower-priority follow-up, not a blocker for this working-point section.",
        ],
    }
    out = OUT / "ppg12_overlay_story_manifest.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n")
    return out


def run_visual_checks(slides: Iterable[Path]) -> None:
    for path in slides:
        img = Image.open(path)
        if img.size != (W, H):
            raise RuntimeError(f"{path} has size {img.size}, expected {(W, H)}")
        # Ensure slide is not blank by checking luminance variance.
        gray = np.asarray(img.convert("L"), dtype=np.float32)
        if float(gray.std()) < 8.0:
            raise RuntimeError(f"{path} appears visually blank")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fixed = make_fixed_overlay_with_ratio()
    slides = [
        slide1_reference_object(),
        slide2_overlay_exercise(fixed),
        slide3_direct_model_confirmation(),
        slide4_same_row_summary(),
    ]
    scripts = write_scripts()
    contact = make_contact_sheet(slides)
    manifest = write_manifest(slides, scripts, contact, fixed)
    run_visual_checks(slides)
    print("Generated:")
    for p in slides:
        print(p)
    print(contact)
    print(manifest)


if __name__ == "__main__":
    main()

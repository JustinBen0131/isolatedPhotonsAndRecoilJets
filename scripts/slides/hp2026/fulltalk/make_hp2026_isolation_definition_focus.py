#!/usr/bin/env python3
"""Local PNG-only candidate for the HP2026 isolation slide.

Redesign goal (Justin, 2026-06-10):
- the reconstructed isolation-energy cut definition was small and hard to find;
  make it the hero object so the audience sees the definition immediately;
- explain the physics, but do not carry analysis-reproduction detail;
- drop the ABCD 2x2 purity grid from this slide. The ABCD/purity region logic
  is already the subject of the very next slide (slide07 in the generator,
  "Purity is measured from data, not assumed", PPG12 Fig. 4), so the grid here
  is redundant and only competes with the isolation definition for space.

This renders a review candidate only. It does not mutate Google Slides.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw

import make_hp2026_fulltalk_candidates as full


ROOT = full.ROOT
OUT_DIR = ROOT / "outputs/manual-20260610-hp2026-isolation-definition-focus"
PNG_PATH = OUT_DIR / "hp2026_isolation_definition_focus_candidate.png"
SCRIPT_PATH = OUT_DIR / "hp2026_isolation_definition_focus_script.md"
MANIFEST_PATH = OUT_DIR / "manifest.json"


def _iso_formula_parts(size: int) -> list[tuple[str, int, int, Path]]:
    """E_T^{iso,reco} < 0.49 GeV + 0.037 E_T with true sub/superscripts."""
    sub = max(18, round(size * 0.52))
    sub_dy = round(size * 0.46)
    sup_dy = round(-size * 0.06)
    return [
        ("E", size, 0, full.TIMES_BOLD),
        ("T", sub, sub_dy, full.TIMES_BOLD),
        ("iso,reco", sub, sup_dy, full.TIMES_BOLD),
        ("  <  0.49 GeV + 0.037 ", size, 0, full.TIMES_BOLD),
        ("E", size, 0, full.TIMES_BOLD),
        ("T", sub, sub_dy, full.TIMES_BOLD),
    ]


def _formula_width(draw: ImageDraw.ImageDraw, parts: list[tuple[str, int, int, Path]]) -> int:
    # Subscript ("T") and superscript ("iso,reco") share the same x column;
    # account for that so centering matches what draw_formula_run renders.
    width = 0
    i = 0
    while i < len(parts):
        text, size, _dy, font_path = parts[i]
        w = full.text_box(draw, text, full.font(font_path, size))[0]
        if text == "T" and i + 1 < len(parts) and parts[i + 1][0] == "iso,reco":
            sup = parts[i + 1]
            sup_w = full.text_box(draw, sup[0], full.font(sup[3], sup[1]))[0]
            width += max(w, sup_w)
            i += 2
            continue
        width += w
        i += 1
    return width


def _draw_iso_formula(draw: ImageDraw.ImageDraw, parts: list[tuple[str, int, int, Path]], xy: tuple[int, int], fill) -> None:
    x, y = xy
    i = 0
    while i < len(parts):
        text, size, dy, font_path = parts[i]
        fnt = full.font(font_path, size)
        draw.text((x, y + dy), text, font=fnt, fill=fill)
        w = full.text_box(draw, text, fnt)[0]
        if text == "T" and i + 1 < len(parts) and parts[i + 1][0] == "iso,reco":
            sup_text, sup_size, sup_dy, sup_font = parts[i + 1]
            sup_fnt = full.font(sup_font, sup_size)
            draw.text((x, y + sup_dy), sup_text, font=sup_fnt, fill=fill)
            sup_w = full.text_box(draw, sup_text, sup_fnt)[0]
            x += max(w, sup_w)
            i += 2
            continue
        x += w
        i += 1


def _draw_centered(draw: ImageDraw.ImageDraw, text: str, cx: float, y: float, fnt, fill) -> None:
    tw, _ = full.text_box(draw, text, fnt)
    draw.text((cx - tw / 2, y), text, font=fnt, fill=fill)


def _draw_et_caption(draw: ImageDraw.ImageDraw, cx: float, y: float) -> None:
    """Centered 'E_T-dependent cutoff chosen for 80% signal efficiency' with a true E_T subscript."""
    main = full.font(full.TIMES, 42)
    sub = full.font(full.TIMES, 27)
    pre, rest = "E", "-dependent cutoff chosen for 80% signal efficiency"
    pre_w = full.text_box(draw, pre, main)[0]
    sub_w = full.text_box(draw, "T", sub)[0]
    rest_w = full.text_box(draw, rest, main)[0]
    x = cx - (pre_w + sub_w + rest_w) / 2
    draw.text((x, y), pre, font=main, fill=full.MUTED)
    draw.text((x + pre_w, y + 15), "T", font=sub, fill=full.MUTED)
    draw.text((x + pre_w + sub_w, y), rest, font=main, fill=full.MUTED)


def draw_isolation_reading_card(base: Image.Image, box: tuple[int, int, int, int]) -> None:
    """Single full-height RHS card: clean legend (what's plotted) over the isolation-cut definition.

    No 'How to read...' header, no redundant subtitle, no per-row sentences. Two labelled
    sections give organisation; the plot legend swatches speak for themselves.
    """
    draw = ImageDraw.Draw(base, "RGBA")
    full.shadow(base, box)
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=(255, 255, 255, 255), outline=(*full.PANEL_EDGE, 255), width=2)
    cx = (x0 + x1) / 2
    ix0 = x0 + 40
    ix1 = x1 - 36

    # ----- Section A: plot legend as colour-coded rows (clean left accent bar per data class) -----
    rows = [
        {"title": "Data, photon-like BDT score", "body": "passes NCB, preselection, and BDT cut", "swatch": (24, 24, 24), "bar": (96, 104, 116), "fill": (255, 255, 255), "marker": "points"},
        {"title": "Data, background-like BDT score", "body": "background-enriched sample below BDT threshold", "swatch": (207, 76, 62), "bar": (207, 76, 62), "fill": (255, 255, 255), "marker": "shade"},
        {"title": "Signal MC", "body": "isolated prompt-photon template", "swatch": (84, 104, 224), "bar": (84, 104, 224), "fill": (255, 255, 255), "marker": "shade"},
    ]
    ry = y0 + 34
    rh, gap = 166, 22
    for row in rows:
        ry0, ry1 = ry, ry + rh
        cy = (ry0 + ry1) // 2
        bar, sw = row["bar"], row["swatch"]
        draw.rounded_rectangle((ix0, ry0, ix1, ry1), radius=12, fill=(*row["fill"], 255), outline=(*bar, 95), width=2)
        draw.rounded_rectangle((ix0, ry0, ix0 + 11, ry1), radius=5, fill=(*bar, 255))
        if row["marker"] == "points":
            for px, py in ((ix0 + 56, cy - 39), (ix0 + 80, cy - 14), (ix0 + 56, cy + 31)):
                draw.ellipse((px - 8, py - 8, px + 8, py + 8), fill=(*sw, 255))
                draw.line((px, py - 21, px, py + 21), fill=(*sw, 160), width=3)
        else:
            draw.rounded_rectangle((ix0 + 38, cy - 35, ix0 + 112, cy + 35), radius=8, fill=(*sw, 70), outline=(*sw, 225), width=3)
        tx = ix0 + 136
        title_size = 40 if "BDT score" in row["title"] else 44
        draw.text((tx, cy - 52), row["title"], font=full.font(full.TIMES_BOLD, title_size), fill=full.INK)
        draw.text((tx, cy + 12), row["body"], font=full.font(full.TIMES, 38), fill=full.MUTED)
        ry += rh + gap

    # ----- divider -----
    dy = ry + 32
    draw.line((ix0, dy, ix1, dy), fill=(224, 230, 237, 255), width=2)

    # ----- Section B: the isolation cut (already defined earlier; here the result only) -----
    head_font = full.font(full.TIMES_BOLD, 44)
    sub_font = full.font(full.TIMES, 44)
    head_txt, sub_txt = "Isolation cut", "derived with signal MC"
    hw = full.text_box(draw, head_txt, head_font)[0]
    sw2 = full.text_box(draw, sub_txt, sub_font)[0]
    arrow_w = 94
    lx = cx - (hw + arrow_w + sw2) / 2
    ly = dy + 56
    draw.text((lx, ly), head_txt, font=head_font, fill=full.PHOTON_DARK)
    ax0, ax1, ay = lx + hw + 26, lx + hw + arrow_w - 26, ly + 30
    draw.line((ax0, ay, ax1, ay), fill=(*full.MUTED, 255), width=5)
    draw.polygon([(ax1 + 13, ay), (ax1 - 8, ay - 12), (ax1 - 8, ay + 12)], fill=(*full.MUTED, 255))
    draw.text((lx + hw + arrow_w, ly), sub_txt, font=sub_font, fill=full.INK)

    # Hero formula on the soft pill.
    size = 62
    parts = _iso_formula_parts(size)
    fw = _formula_width(draw, parts)
    fy = dy + 184
    px_pad, py_pad = 50, 30
    pill = (cx - fw / 2 - px_pad, fy - py_pad, cx + fw / 2 + px_pad, fy + size + py_pad)
    draw.rounded_rectangle(pill, radius=18, fill=(255, 251, 240, 255), outline=(*full.PHOTON, 150), width=2)
    _draw_iso_formula(draw, parts, (cx - fw / 2, fy), fill=full.INK)

    # Caption: ET-dependent cutoff, 80% signal efficiency.
    _draw_et_caption(draw, cx, dy + 324)


def build_slide() -> Image.Image:
    title = "Isolation defines the photon sample"
    subtitle = "Quiet prompt-like candidates separate from nearby jet activity."
    img = full.base_slide_hp2026_main(title, subtitle)
    full.add_top_right_sphenix_logo_like_slide2(img)
    # Full-height plot card; tight crop fills it. Narrower than before so the (height-limited)
    # square plot keeps its size while the right card gains room for conference-scale text.
    left_box = (132, 258, 1150, 1308)
    full.place_figure_tight(img, "fig3_isolation", left_box, inset=24, pad=24)
    draw = ImageDraw.Draw(img, "RGBA")
    full.outer_card_sidebar(draw, left_box, full.TEAL)
    # Wider full-height RHS card: colour-coded legend + the isolation-cut result.
    draw_isolation_reading_card(img, (1208, 258, 2390, 1308))
    full.draw_hp2026_identity_footer(img)
    return img


def write_script() -> None:
    SCRIPT_PATH.write_text(
        """# HP2026 Slide Script - Isolation defines the photon sample

Now that the photon-ID BDT has defined a prompt-like cluster, the next question is whether the region around that cluster is quiet. That is what isolation measures. A real prompt photon comes straight out of the hard scattering, so it should not drag a lot of nearby hadronic activity along with it, while fragmentation photons and neutral-meson backgrounds tend to sit inside a busier jet environment.

The plot on the left is the reconstructed isolation-energy distribution. The black points are the tight-ID data candidates after preselection. The red shaded shape is non-tight-ID data, so it is a background-enriched sideband. The blue shaded shape is the tight-ID prompt-photon signal simulation. The point to take away is the contrast: the prompt-photon template piles up at low isolation energy, while the background-enriched shape carries a long high-isolation tail.

The definition on the right is the one number to remember from this slide. A candidate counts as isolated when its reconstructed isolation energy falls below 0.49 GeV plus 0.037 times its transverse energy. I tune that threshold to keep about eighty percent prompt-photon isolation efficiency, so it is a deliberately loose, high-efficiency cut rather than an aggressive one.

So isolation does two things at once. It cleans the selected sample, and it gives me a second handle on the data. On the next slide I use that handle together with the photon-ID axis to build the sidebands and measure the purity directly from data.
""",
        encoding="utf-8",
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    img = build_slide()
    img.convert("RGB").save(PNG_PATH, "PNG")
    full.write_hp2026_main_header_spec(PNG_PATH)
    write_script()
    MANIFEST_PATH.write_text(
        json.dumps(
            {
                "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
                "google_slides_mutation": False,
                "deck": "HPslides_v1 (live deck slide 'Isolation defines the photon sample')",
                "source_baseline": "scripts/slides/hp2026/fulltalk/make_hp2026_fulltalk_candidates.py::slide06",
                "figure_source": "PPG12 paper Fig. 3 (reconstructed isolation energy), page 9 -> assets/paper_figures/fig3_isolation.png",
                "output_png": str(PNG_PATH.relative_to(ROOT)),
                "companion_script": str(SCRIPT_PATH.relative_to(ROOT)),
                "slide_intent": "Make the reconstructed isolation-energy cut definition the hero object; explain isolation physics without reproduction-level detail.",
                "changed": [
                    "Two full-height cards use the same vertical band as the adjacent purity slide: plot card (132,258,1150,1308) + wider RHS reading card (1208,258,2390,1308).",
                    "Plot card narrowed to reclaim wasted side margin; the height-limited square plot keeps its size while the RHS gains room for larger text.",
                    "Conference-readability pass: legend titles 38, bodies 30, lead 38, formula 60, caption 30.",
                    "Legend rows are now white with a clean left accent bar per data class (gray / red / indigo) and a matching outline, removing the extra pastel fills for projector readability.",
                    "Isolation cut section (cut already defined earlier in the deck): lead 'Isolation cut -> derived with signal MC', hero formula on a soft pill, caption 'E_T-dependent cutoff chosen for 80% signal efficiency'.",
                ],
                "removed": [
                    "'How to read the isolation distribution' header and its italic subtitle sentence (redundant).",
                    "Two-card RHS split; now one unified full-height card.",
                    "ABCD 2x2 purity-region grid (redundant: the ABCD/purity story is the very next slide, PPG12 Fig. 4).",
                    "non-isolated sideband offset detail (cut + 0.8 GeV) and the spoken bridge sentence (latter lives in the speaker script).",
                ],
                "kept": [
                    "Left PPG12 Fig. 3 isolation-energy figure unchanged.",
                    "Top-right 'How to read the isolation distribution' legend panel unchanged.",
                    "~80% isolation-efficiency motivation, as a single compact line.",
                ],
                "labels_note": "Figure carries 'sPHENIX Internal'. This is an internal/working candidate; the public HP version still needs the released/Preliminary-labeled isolation figure before it is public-ready.",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(PNG_PATH)
    print(SCRIPT_PATH)
    print(MANIFEST_PATH)


if __name__ == "__main__":
    main()

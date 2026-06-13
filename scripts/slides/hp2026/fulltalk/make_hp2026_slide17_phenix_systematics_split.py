#!/usr/bin/env python3
"""Local PNG candidates splitting Slide 17 into main and backup slides.

This generator intentionally does not mutate Google Slides.  It reuses the
existing HP2026 fulltalk style and paper-figure assets, then writes local PNG
previews and companion practice scripts.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from PIL import Image, ImageDraw

THIS = Path(__file__).resolve()
if str(THIS.parent) not in sys.path:
    sys.path.insert(0, str(THIS.parent))

import make_hp2026_fulltalk_candidates as ft  # noqa: E402


ROOT = ft.ROOT
OUTPUT = ROOT / "outputs/manual-20260611-slides15-17-header-normalized/phenix_comparison"
SCRIPT_DIR = OUTPUT / "speaker_scripts"


def _save_script(stem: str, title: str, body: str) -> Path:
    SCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    path = SCRIPT_DIR / f"{stem}_script.md"
    path.write_text(f"# {title}\n\n{body.strip()}\n", encoding="utf-8")
    return path


def _draw_side_bar_panel(
    img: Image.Image,
    box: tuple[int, int, int, int],
    accent: tuple[int, int, int],
    title: str,
    body_blocks: list[tuple[str, str]],
    *,
    title_size: int = 42,
    label_size: int = 31,
    body_size: int = 29,
    y_gap: int = 28,
) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    ft.shadow(img, box)
    draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 255), outline=(*ft.CARD_EDGE, 255), width=2)
    draw.rounded_rectangle((x0, y0, x0 + 13, y1), radius=6, fill=(*accent, 255))
    draw.text((x0 + 38, y0 + 28), title, font=ft.font(ft.TIMES_BOLD, title_size), fill=ft.INK)
    y = y0 + 96
    for lead, text in body_blocks:
        lead_font = ft.font(ft.TIMES_BOLD, label_size)
        body_font = ft.font(ft.TIMES, body_size)
        draw.text((x0 + 42, y), lead, font=lead_font, fill=ft.INK)
        lead_w, lead_h = ft.text_box(draw, lead, lead_font)
        y_after = ft.draw_wrapped(
            draw,
            text,
            (x0 + 56 + lead_w, y + 2),
            x1 - x0 - lead_w - 90,
            body_font,
            fill=ft.MUTED,
            line_gap=5,
        )
        y = max(y + lead_h, y_after) + y_gap


def _draw_comparison_readout_panel(
    img: Image.Image,
    box: tuple[int, int, int, int],
) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    ft.shadow(img, box)
    draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 255), outline=(*ft.CARD_EDGE, 255), width=2)
    draw.rounded_rectangle((x0, y0, x0 + 13, y1), radius=6, fill=(*ft.SPHENIX_BLUE, 255))
    draw.text((x0 + 38, y0 + 26), "Comparison readout", font=ft.font(ft.TIMES_BOLD, 44), fill=ft.INK)
    draw.line((x0 + 38, y0 + 86, x1 - 38, y0 + 86), fill=(222, 229, 237, 255), width=2)

    rows = [
        (
            "sPHENIX",
            ft.SPHENIX_BLUE,
            (241, 248, 254),
            "isolated prompt photons",
            "|η| < 0.7",
        ),
        (
            "PHENIX 2012",
            (177, 82, 201),
            (250, 244, 252),
            "direct photons, √s = 200 GeV",
            "PRD 86, 072008; no isolation; |η| < 0.25",
        ),
        (
            "Correction",
            (104, 99, 214),
            (246, 245, 255),
            "bin-width + rapidity-density rescaling",
            "puts PHENIX on the |η| < 0.7 basis",
        ),
    ]
    row_h = 104
    gap = 8
    y = y0 + 108
    label_w = 328
    for label, color, fill, main, sub in rows:
        row = (x0 + 38, y, x1 - 38, y + row_h)
        lf = ft.font(ft.TIMES_BOLD, 46)
        mf = ft.font(ft.TIMES_BOLD, 41)
        sf = ft.font(ft.TIMES, 35)
        draw.text((row[0] + 8, row[1] + row_h // 2), label, font=lf, fill=color, anchor="lm")
        draw.text((row[0] + label_w, row[1] + 16), main, font=mf, fill=ft.INK)
        draw.text((row[0] + label_w, row[1] + 61), sub, font=sf, fill=ft.MUTED)
        y += row_h + gap


def _draw_establishes_panel(
    img: Image.Image,
    box: tuple[int, int, int, int],
) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    ft.shadow(img, box)
    draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 255), outline=(*ft.CARD_EDGE, 255), width=2)
    draw.rounded_rectangle((x0, y0, x0 + 13, y1), radius=6, fill=(*ft.PHOTON, 255))
    draw.text((x0 + 38, y0 + 26), "What this comparison establishes", font=ft.font(ft.TIMES_BOLD, 42), fill=ft.INK)
    draw.line((x0 + 38, y0 + 86, x1 - 38, y0 + 86), fill=(229, 232, 236, 255), width=2)

    rows = [
        (
            "Consistency",
            "corrected PHENIX / sPHENIX is compatible with unity.",
        ),
        (
            "Reach",
            "sPHENIX extends acceptance and photon-energy reach.",
        ),
        (
            "Systematics",
            "low_ET_high_ET",
        ),
    ]
    y = y0 + 116
    row_h = 112
    gap = 14
    body_fill = (43, 51, 65)
    label_font = ft.font(ft.TIMES_BOLD, 46)
    # Body is a fixed column after the WIDEST label (like the top box), so the gap from
    # label to body is uniform and Consistency/Systematics no longer crowd their text.
    body_x_off = 94 + max(ft.text_box(draw, head, label_font)[0] for head, _ in rows) + 34
    for idx, (head, body) in enumerate(rows, start=1):
        row = (x0 + 38, y, x1 - 38, y + row_h)
        fill = (255, 252, 244) if idx == 3 else (250, 252, 254)
        draw.rounded_rectangle(row, radius=8, fill=(*fill, 255), outline=(224, 231, 238, 255), width=1)
        badge_r = 24
        cx = row[0] + 42
        cy = row[1] + row_h // 2
        draw.ellipse((cx - badge_r, cy - badge_r, cx + badge_r, cy + badge_r), fill=(*ft.PHOTON, 255))
        num = str(idx)
        nf = ft.font(ft.TIMES_BOLD, 34)
        draw.text((cx, cy), num, font=nf, fill=(255, 255, 255), anchor="mm")
        label_x = row[0] + 94
        draw.text((label_x, cy), head, font=label_font, fill=ft.INK, anchor="lm")
        body_x = row[0] + body_x_off
        body_font = ft.font(ft.TIMES, 41)
        body_max_w = row[2] - body_x - 28
        if body == "low_ET_high_ET":
            _draw_systematics_et_phrase(
                draw,
                row,
                body_x,
                body_max_w,
                body_font,
                fill=body_fill,
            )
        else:
            _draw_wrapped_centered_in_row(
                draw,
                body,
                body_x,
                row,
                body_max_w,
                body_font,
                fill=body_fill,
                line_gap=5,
            )
        y += row_h + gap


def _draw_wrapped_centered_in_row(
    draw: ImageDraw.ImageDraw,
    text: str,
    x: int,
    row: tuple[int, int, int, int],
    max_width: int,
    fnt,
    *,
    fill: tuple[int, int, int],
    line_gap: int = 4,
) -> None:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        trial = word if not current else f"{current} {word}"
        if ft.text_box(draw, trial, fnt)[0] <= max_width:
            current = trial
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    line_heights = [ft.text_box(draw, line, fnt)[1] for line in lines]
    total_h = sum(line_heights) + line_gap * max(0, len(lines) - 1)
    y = row[1] + ((row[3] - row[1]) - total_h) // 2 - 1
    for line, h in zip(lines, line_heights):
        draw.text((x, y), line, font=fnt, fill=fill)
        y += h + line_gap


def _draw_systematics_et_phrase(
    draw: ImageDraw.ImageDraw,
    row: tuple[int, int, int, int],
    x: int,
    max_width: int,
    main_font,
    *,
    fill: tuple[int, int, int],
) -> None:
    sub_font = ft.font(ft.TIMES, 29)
    lines = [
        [
            ("low E", main_font, 0),
            ("T", sub_font, 15),
            (": purity dominates", main_font, 0),
        ],
        [
            ("high E", main_font, 0),
            ("T", sub_font, 15),
            (": energy scale/resolution dominate", main_font, 0),
        ],
    ]

    def line_width(parts: list[tuple[str, object, int]]) -> int:
        return sum(ft.text_box(draw, text, fnt)[0] for text, fnt, _ in parts)

    main_h = ft.text_box(draw, "low E", main_font)[1]
    line_gap = 16
    total_h = main_h * len(lines) + line_gap * (len(lines) - 1)
    y = row[1] + ((row[3] - row[1]) - total_h) // 2 - 1
    for parts in lines:
        cursor = x
        width = line_width(parts)
        if width > max_width:
            # Keep the phrase bounded even if a font substitution changes metrics.
            cursor = x
        for text, fnt, dy in parts:
            draw.text((cursor, y + dy), text, font=fnt, fill=fill)
            cursor += ft.text_box(draw, text, fnt)[0]
        y += main_h + line_gap


def _draw_takeaway_panel(
    img: Image.Image,
    box: tuple[int, int, int, int],
    title: str,
    lines: list[str],
    *,
    accent: tuple[int, int, int] = ft.PHOTON,
) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    ft.shadow(img, box)
    draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 255), outline=(*ft.CARD_EDGE, 255), width=2)
    draw.rounded_rectangle((x0, y0, x0 + 13, y1), radius=6, fill=(*accent, 255))
    draw.text((x0 + 38, y0 + 26), title, font=ft.font(ft.TIMES_BOLD, 43), fill=ft.INK)

    y = y0 + 106
    bullet_font = ft.font(ft.TIMES_BOLD, 34)
    line_font = ft.font(ft.TIMES, 33)
    for line in lines:
        draw.ellipse((x0 + 44, y + 13, x0 + 62, y + 31), fill=(*accent, 255))
        y = ft.draw_wrapped(
            draw,
            line,
            (x0 + 78, y),
            x1 - x0 - 118,
            line_font if not line.startswith("Dominant") else bullet_font,
            fill=ft.INK if line.startswith("Dominant") else ft.MUTED,
            line_gap=7,
        ) + 30


def _draw_systematics_card(
    img: Image.Image,
    box: tuple[int, int, int, int],
    title: str,
    rows: list[tuple[tuple[int, int, int], str, str]],
) -> None:
    draw = ImageDraw.Draw(img, "RGBA")
    x0, y0, x1, y1 = box
    ft.shadow(img, box)
    draw.rounded_rectangle(box, radius=10, fill=(255, 255, 255, 255), outline=(*ft.CARD_EDGE, 255), width=2)
    draw.rounded_rectangle((x0, y0, x0 + 13, y1), radius=6, fill=(*ft.TEAL, 255))
    draw.text((x0 + 38, y0 + 26), title, font=ft.font(ft.TIMES_BOLD, 42), fill=ft.INK)
    y = y0 + 102
    for color, head, body in rows:
        draw.rounded_rectangle((x0 + 42, y + 6, x0 + 70, y + 34), radius=6, fill=(*color, 255))
        draw.text((x0 + 88, y), head, font=ft.font(ft.TIMES_BOLD, 31), fill=ft.INK)
        y = ft.draw_wrapped(
            draw,
            body,
            (x0 + 88, y + 40),
            x1 - x0 - 132,
            ft.font(ft.TIMES, 29),
            fill=ft.MUTED,
            line_gap=5,
        ) + 36


def render_main_phenix_slide() -> tuple[Path, Path]:
    title = "RHIC comparison: isolated prompt photons"
    img = ft.base_slide_hp2026_main(title)
    draw = ImageDraw.Draw(img, "RGBA")
    ft.add_top_right_sphenix_logo_like_slide2(img)

    fig_box = (86, 270, 1286, 1308)
    fig_placed, fig_panel = ft.place_figure_snug_panel(
        img,
        "fig9_phenix",
        fig_box,
        pad=18,
        crop_pad=8,
        trim_bottom=0,
    )

    _draw_comparison_readout_panel(
        img,
        (1248, 270, 2390, 724),
    )

    _draw_establishes_panel(
        img,
        (1248, 764, 2390, 1308),
    )

    ft.draw_hp2026_identity_footer(img)

    png = OUTPUT / "slide15_rhic_comparison_header_normalized.png"
    OUTPUT.mkdir(parents=True, exist_ok=True)
    img.convert("RGB").save(png, "PNG")
    ft.write_hp2026_main_header_spec(png)
    script = _save_script(
        "slide17_main_phenix_comparison_no_systematics",
        title,
        """
        This is the short RHIC-context check.  The blue points are the corrected
        sPHENIX isolated prompt-photon cross section in |eta gamma| less than
        0.7.  The PHENIX reference is the 2012 direct-photon measurement at
        midrapidity, published in Physical Review D 86, 072008.  That
        measurement had no isolation requirement, used |eta gamma| less than
        0.25, and reported bin-center cross sections.

        The purple points show the PHENIX result after the comparison correction
        used in the paper: a bin-width correction plus a rapidity-density
        rescaling to |eta gamma| less than 0.7.

        The main point is the ratio panel: corrected PHENIX divided by sPHENIX
        is consistent with unity within uncertainties.  So the corrected
        sPHENIX p+p baseline agrees with RHIC legacy photons, while extending
        the rapidity acceptance and transverse-energy reach.

        For the uncertainty pattern, I only want to quote the headline here:
        purity is the important low-transverse-energy contribution, while
        energy scale and resolution matter most at high transverse energy.  The
        full budget is moved to the backup slide.
        """,
    )
    return png, script


def render_backup_systematics_slide() -> tuple[Path, Path]:
    img = ft.base_slide(
        "Backup: systematic uncertainty budget",
        "The full source-by-source budget supports the compact main-slide statement.",
    )
    ft.add_top_right_sphenix_logo_like_slide2(img)

    fig_placed, fig_panel = ft.place_figure_snug_panel(
        img,
        "fig7_systematics",
        (88, 318, 1748, 1292),
        pad=20,
        crop_pad=8,
        trim_bottom=0,
    )

    _draw_systematics_card(
        img,
        (1800, 318, 2392, 808),
        "How to read it",
        [
            (ft.PHOTON, "Colored lines", "one systematic source varied through the corrected cross section."),
            (ft.INK, "Black envelope", "total systematic uncertainty, including luminosity normalization."),
            (ft.TEAL, "Y-axis", "fractional shift from the nominal corrected result."),
        ],
    )

    _draw_systematics_card(
        img,
        (1800, 850, 2392, 1292),
        "Dominant pattern",
        [
            (ft.PHOTON_DARK, "Low transverse energy", "purity and closure terms are the largest residual-background controls."),
            (ft.SPHENIX_BLUE, "High transverse energy", "energy scale and resolution dominate because the spectrum is steeply falling."),
        ],
    )

    ft.draw_hp2026_identity_footer(img)

    png = OUTPUT / "slide17_backup_systematics_budget.png"
    OUTPUT.mkdir(parents=True, exist_ok=True)
    img.convert("RGB").save(png, "PNG")
    script = _save_script(
        "slide17_backup_systematics_budget",
        "Backup: systematic uncertainty budget",
        """
        This backup slide is here in case someone asks what controls the
        uncertainty budget.  Each colored line corresponds to varying one
        systematic source and propagating that variation through the corrected
        cross section.  The black envelope is the total systematic uncertainty,
        including the luminosity normalization.

        The physics-readable structure is that the low-transverse-energy bins
        are most sensitive to purity and closure, while the high-transverse
        energy region is most sensitive to energy scale and resolution.  That
        is why the main-talk comparison slide only quotes the compact
        low-energy versus high-energy pattern.
        """,
    )
    return png, script


def write_manifest(outputs: list[tuple[str, Path, Path]]) -> Path:
    manifest = {
        "purpose": "Local split of Slide 17 into main PHENIX comparison and backup systematics slides.",
        "google_slides_mutated": False,
        "source_generator_reused": str(
            (ROOT / "scripts/slides/hp2026/fulltalk/make_hp2026_fulltalk_candidates.py").relative_to(ROOT)
        ),
        "source_assets": {
            "phenix_comparison": str(ft.figure_path("fig9_phenix").relative_to(ROOT)),
            "systematics_budget": str(ft.figure_path("fig7_systematics").relative_to(ROOT)),
        },
        "outputs": [
            {"role": role, "png": str(png.relative_to(ROOT)), "script": str(script.relative_to(ROOT))}
            for role, png, script in outputs
        ],
    }
    path = OUTPUT / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return path


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    outputs = [
        ("main_slide17", *render_main_phenix_slide()),
        ("backup_systematics", *render_backup_systematics_slide()),
    ]
    manifest = write_manifest(outputs)
    for role, png, script in outputs:
        print(f"{role}: {png}")
        print(f"{role}_script: {script}")
    print(f"manifest: {manifest}")


if __name__ == "__main__":
    main()

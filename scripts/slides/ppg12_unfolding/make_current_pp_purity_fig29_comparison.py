#!/usr/bin/env python3
"""Build a PPG12 Fig. 29 style pp purity comparison slide.

Left: current-IAN Fig. 29 left-panel crop.
Right: current recovered pp data raw and leakage-corrected ABCD purity points.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont


ROOT_DIR = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
OUT_DIR = ROOT_DIR / "dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_basev3e_20260611/purity_current_pp"
POINTS_CSV = OUT_DIR / "current_pp_abcd_purity_raw_vs_leakage_corrected_points.csv"
MISMATCH_DIAGNOSTIC = OUT_DIR / "ppg12_fig29_purity_input_mismatch_diagnostic.md"
IAN_PDF = ROOT_DIR / "usefulDocs/PPG12_analysis_note_2026-05-21_v4_current_IAN.pdf"
IAN_PAGE_RENDER = OUT_DIR / "ppg12_current_ian_purity_page-039.png"
IAN_FIG29_LEFT = OUT_DIR / "ppg12_current_ian_fig29_left_purity.png"
PDFTOPPM = Path("/Users/patsfan753/.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/pdftoppm")


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
    ]
    for candidate in candidates:
        if Path(candidate).exists():
            return ImageFont.truetype(candidate, size=size)
    return ImageFont.load_default()


def ensure_ian_fig29_crop() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if IAN_FIG29_LEFT.exists() and IAN_FIG29_LEFT.stat().st_size > 50_000:
        return
    if not IAN_PAGE_RENDER.exists() or IAN_PAGE_RENDER.stat().st_size < 100_000:
        prefix = OUT_DIR / "ppg12_current_ian_purity_page"
        subprocess.run(
            [str(PDFTOPPM), "-png", "-f", "39", "-l", "39", "-r", "220", str(IAN_PDF), str(prefix)],
            check=True,
        )
    page = Image.open(IAN_PAGE_RENDER).convert("RGB")
    # Current-IAN page 39, Figure 29 left panel only.
    page.crop((285, 245, 920, 885)).save(IAN_FIG29_LEFT)


def trim_white_margins(img: Image.Image, tolerance: int = 248, pad: int = 8) -> Image.Image:
    rgb = img.convert("RGB")
    arr = np.asarray(rgb)
    mask = np.any(arr < tolerance, axis=2)
    if not mask.any():
        return rgb
    ys, xs = np.where(mask)
    left = max(int(xs.min()) - pad, 0)
    upper = max(int(ys.min()) - pad, 0)
    right = min(int(xs.max()) + pad + 1, rgb.width)
    lower = min(int(ys.max()) + pad + 1, rgb.height)
    return rgb.crop((left, upper, right, lower))


def read_points() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with POINTS_CSV.open() as handle:
        for row in csv.DictReader(handle):
            parsed = {key: float(value) if key != "correction_ok" else float(value == "True") for key, value in row.items()}
            rows.append(parsed)
    return rows


def render_current_panel() -> tuple[Image.Image, dict]:
    rows_all = read_points()
    rows = [r for r in rows_all if r["pt_lo"] >= 10 and r["pt_hi"] <= 26]
    hidden = [r for r in rows_all if r["pt_hi"] > 26]

    x = np.array([(r["pt_lo"] + r["pt_hi"]) / 2.0 for r in rows])
    ex = np.array([(r["pt_hi"] - r["pt_lo"]) / 2.0 for r in rows])
    y_raw = np.array([r["raw"] for r in rows])
    y_corr = np.array([r["corrected"] for r in rows])
    ey_raw = np.array([r["raw_err"] for r in rows])
    ey_corr = np.array([r["corrected_err"] for r in rows])

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig, ax = plt.subplots(figsize=(5.75, 5.75), dpi=230)
    fig.patch.set_facecolor("white")
    ax.errorbar(
        x,
        y_corr,
        xerr=ex,
        yerr=ey_corr,
        linestyle="None",
        marker="o",
        markersize=7.5,
        markerfacecolor="white",
        markeredgecolor="#144dff",
        markeredgewidth=1.5,
        ecolor="#144dff",
        elinewidth=1.1,
        capsize=2.5,
        label="w/ sig. leak. corr.",
    )
    ax.errorbar(
        x,
        y_raw,
        xerr=ex,
        yerr=ey_raw,
        linestyle="None",
        marker="o",
        markersize=7.0,
        markerfacecolor="black",
        markeredgecolor="black",
        ecolor="black",
        elinewidth=1.1,
        capsize=2.5,
        label="w/o sig. leak. corr.",
    )
    ax.set_xlim(10, 35)
    ax.set_ylim(0, 1.2)
    ax.set_xlabel(r"$E_T^{\gamma,\mathrm{rec}}$ [GeV]", fontsize=16)
    ax.set_ylabel("Purity", fontsize=17)
    ax.tick_params(which="both", labelsize=13, length=5)
    ax.minorticks_on()
    ax.text(0.07, 0.95, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="top", fontsize=13.5)
    ax.text(0.07, 0.885, r"$p{+}p$ $\sqrt{s}=200$ GeV", transform=ax.transAxes, ha="left", va="top", fontsize=12.8)
    ax.text(0.07, 0.820, "RecoilJets ABCD diagnostic", transform=ax.transAxes, ha="left", va="top", fontsize=12.4)
    ax.legend(loc="lower left", bbox_to_anchor=(0.07, 0.08), frameon=False, fontsize=12.0, handlelength=1.4)
    fig.subplots_adjust(left=0.17, right=0.98, top=0.98, bottom=0.15)
    panel_path = OUT_DIR / "current_pp_fig29_style_purity_overlay_panel.png"
    fig.savefig(panel_path, dpi=230)
    plt.close(fig)
    return Image.open(panel_path).convert("RGB"), {
        "points_csv": str(POINTS_CSV),
        "drawn_bins": [[r["pt_lo"], r["pt_hi"]] for r in rows],
        "hidden_bins": [[r["pt_lo"], r["pt_hi"], r["a"]] for r in hidden],
        "style": "PPG12 Figure 29 left-panel convention, no fit band on current output",
    }


def fit_on_white(img: Image.Image, size: tuple[int, int]) -> Image.Image:
    card = Image.new("RGBA", size, (255, 255, 255, 255))
    scale = min(size[0] / img.width, size[1] / img.height)
    fitted = img.resize((max(1, int(img.width * scale)), max(1, int(img.height * scale))), Image.Resampling.LANCZOS)
    card.alpha_composite(fitted.convert("RGBA"), ((size[0] - fitted.width) // 2, (size[1] - fitted.height) // 2))
    return card


def draw_arrow_bullets(
    draw: ImageDraw.ImageDraw,
    *,
    xy: tuple[int, int],
    bullets: list[str],
    max_width: int,
    body_font: ImageFont.FreeTypeFont,
    font_px: int,
) -> tuple[int, list[dict]]:
    x, y = xy
    nodes: list[dict] = []
    arrow_w = 28
    line_gap = 12
    bullet_gap = 20
    for idx, text in enumerate(bullets):
        words = text.split()
        lines: list[str] = []
        current = ""
        for word in words:
            trial = word if not current else f"{current} {word}"
            if draw.textlength(trial, font=body_font) <= max_width - 70:
                current = trial
            else:
                if current:
                    lines.append(current)
                current = word
        if current:
            lines.append(current)
        top_y = y
        arrow_y = y + 15
        draw.polygon([(x, arrow_y), (x + arrow_w, arrow_y + 13), (x, arrow_y + 26)], fill=(36, 104, 168))
        text_x = x + 50
        for line in lines:
            draw.text((text_x, y), line, font=body_font, fill=(23, 32, 51))
            y += font_px + line_gap
        bbox = (x, top_y, x + max_width, y - line_gap)
        nodes.append(
            {
                "kind": "text",
                "name": f"arrow bullet {idx + 1}",
                "role": "audience",
                "text": text,
                "bbox": list(bbox),
                "font_px": font_px,
                "title_axis_align": "left",
            }
        )
        y += bullet_gap
    return y, nodes


def compose_slide(current_panel: Image.Image, current_meta: dict) -> Path:
    ensure_ian_fig29_crop()
    reference = trim_white_margins(Image.open(IAN_FIG29_LEFT).convert("RGB"))
    current = trim_white_margins(current_panel)

    slide = Image.new("RGB", (2560, 1440), "white")
    draw = ImageDraw.Draw(slide)
    nodes: list[dict] = []

    title_font = font(76, bold=True)
    lead_font = font(52)
    lead_bold = font(52, bold=True)
    panel_label_font = font(42, bold=True)
    bullet_font = font(44)

    title = "Photon Purity — Signal-Leakage Correction"
    title_x = 70
    draw.text((title_x, 42), title, font=title_font, fill=(0, 0, 0))
    nodes.append({"kind": "text", "name": "slide title", "role": "title", "text": title, "bbox": list(draw.textbbox((title_x, 42), title, font=title_font)), "font_px": 76})

    lead = "Data ABCD purity"
    rest = "with and without signal-leakage correction"
    lead_y = 168
    draw.text((title_x, lead_y), lead, font=lead_bold, fill=(0, 0, 0))
    lead_w = draw.textlength(lead, font=lead_bold)
    draw.text((title_x + int(lead_w) + 22, lead_y), rest, font=lead_font, fill=(0, 0, 0))
    nodes.append({"kind": "text", "name": "lead sentence", "role": "audience", "text": f"{lead} {rest}", "bbox": [title_x, lead_y, 1680, lead_y + 62], "font_px": 52})

    card_w, card_h = 940, 840
    left_xy = (240, 292)
    right_xy = (1380, 292)
    left_card = fit_on_white(reference, (card_w, card_h))
    right_card = fit_on_white(current, (card_w, card_h))
    slide.paste(left_card.convert("RGB"), left_xy, left_card)
    slide.paste(right_card.convert("RGB"), right_xy, right_card)
    nodes.extend(
        [
            {"kind": "image", "name": "left purity reference", "bbox": [left_xy[0], left_xy[1], left_xy[0] + card_w, left_xy[1] + card_h], "repeated_group": "purity_plots"},
            {"kind": "image", "name": "right purity current", "bbox": [right_xy[0], right_xy[1], right_xy[0] + card_w, right_xy[1] + card_h], "repeated_group": "purity_plots"},
        ]
    )

    left_label = "PPG12 IAN Fig. 29 reference"
    right_label = "Current RecoilJets diagnostic"
    label_y = left_xy[1] + card_h + 8
    for xy, label in [(left_xy, left_label), (right_xy, right_label)]:
        draw.text((xy[0] + (card_w - int(draw.textlength(label, font=panel_label_font))) // 2, label_y), label, font=panel_label_font, fill=(0, 0, 0))
    nodes.extend(
        [
            {"kind": "text", "name": "left panel label", "role": "audience", "text": left_label, "bbox": [left_xy[0], label_y, left_xy[0] + card_w, label_y + 54], "font_px": 42},
            {"kind": "text", "name": "right panel label", "role": "audience", "text": right_label, "bbox": [right_xy[0], label_y, right_xy[0] + card_w, label_y + 54], "font_px": 42},
        ]
    )

    _, bullet_nodes = draw_arrow_bullets(
        draw,
        xy=(70, 1230),
        bullets=[
            "Same visual convention as the current IAN, blue open markers are leakage-corrected and black filled markers are raw ABCD.",
            "Current pp points use the recovered complete table-QA data ROOT; the 26-35 GeV bin is retained in CSV but omitted here because A=39.",
        ],
        max_width=2300,
        body_font=bullet_font,
        font_px=44,
    )
    nodes.extend(bullet_nodes)

    out = OUT_DIR / "slide_current_pp_purity_fig29_comparison.png"
    slide.save(out)
    manifest = {
        "output_png": str(out),
        "generator": str(Path(__file__).resolve()),
        "ian_source_pdf": str(IAN_PDF),
        "ian_source_date": "2026-05-21 v4 current IAN",
        "ian_reference_page": 39,
        "ian_reference_crop": str(IAN_FIG29_LEFT),
        "current_output": current_meta,
        "truthfulness_note": "The current RHS uses RecoilJets ABCD diagnostic counts and signal-MC leakage fractions. It is not an apples-to-apples PPG12 Fig. 29 purity reproduction because the current ROOT lacks the PPG12 Photon_final/CalculatePhotonYield contract and high-ET purity bins.",
    }
    out.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    out.with_suffix(".layout_nodes.json").write_text(
        json.dumps(
            {
                "minimum_audience_font_px": 34,
                "minimum_title_font_px": 54,
                "minimum_plot_annotation_font_px": 24,
                "nodes": nodes,
            },
            indent=2,
        )
        + "\n"
    )
    out.with_suffix(".speaker_script.md").write_text(
        "\n".join(
            [
                "# Speaker Script",
                "",
                "This slide compares the current-IAN photon purity figure to the same raw-versus-leakage-corrected ABCD purity diagnostic from the recovered pp output.",
                "The left panel is the current PPG12 IAN Figure 29 reference. The right panel uses the same marker convention, with blue open points for the signal-leakage-corrected purity and black filled points for the raw ABCD purity.",
                "The current output is a point-level diagnostic, not yet a fit-replacement slide. The high-pT 26-35 GeV bin is kept in the sidecar CSV but omitted here because the recovered table-QA root has only 39 A-region entries there.",
                "",
            ]
        )
    )
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--force-recoiljets-diagnostic",
        action="store_true",
        help="Render the non-apples-to-apples RecoilJets ABCD diagnostic anyway.",
    )
    args = parser.parse_args()
    if not args.force_recoiljets_diagnostic:
        raise SystemExit(
            "Refusing to render this as a PPG12 Fig. 29 comparison: current pp output "
            "uses RecoilJets ABCD diagnostic histograms and lacks the PPG12 Photon_final/"
            "CalculatePhotonYield purity contract. See "
            f"{MISMATCH_DIAGNOSTIC}. Re-run with --force-recoiljets-diagnostic only "
            "for a clearly labeled diagnostic slide."
        )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    panel, meta = render_current_panel()
    out = compose_slide(panel, meta)
    print(out)


if __name__ == "__main__":
    main()

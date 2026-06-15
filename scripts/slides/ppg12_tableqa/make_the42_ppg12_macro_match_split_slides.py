#!/usr/bin/env python3
"""Build split PPG12 reference vs THE-42 macro-match QA slide candidates."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt
from PIL import Image, ImageChops, ImageDraw, ImageFont, ImageOps

import sys

SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = next(p for p in SCRIPT_PATH.parents if (p / "AGENTS.md").exists())
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.append(str(SCRIPTS_DIR))

from slides.common.slide_defaults import SLIDE_HEIGHT_PX, SLIDE_WIDTH_PX  # noqa: E402


CAMPAIGN = REPO_ROOT / "dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_basev3e_20260611"
OUT_DIR = CAMPAIGN / "comparison_slide"
IAN_DIR = OUT_DIR / "ian_extract"
REF_DIR = IAN_DIR / "clean_ref_panels"
ROOT_DIR = CAMPAIGN / "merged_roots"
INCLUSIVE_CACHE = CAMPAIGN / "inclusive_sample_hist_cache/the42_ppg12_tableqa_v1_inclusive_sample_projectx_hists.json"
PLOTTER_PATH = REPO_ROOT / "scripts/plotting/pp_currentian/make_the42_ppg12_tableqa_v1_tables.py"

PPG12_IAN_PAGE19 = IAN_DIR / "ppg12_ian_page19_300dpi.png"
PPG12_E11_REF = REF_DIR / "ppg12_fig13_e11_to_e33_clean.png"
PPG12_BDT_REF = REF_DIR / "ppg12_fig13_bdt_clean.png"
PPG12_IAN_LABEL = "PPG12 IAN v4 (May 21, 2026)"


def load_plotter():
    spec = importlib.util.spec_from_file_location("tableqa_plotter", PLOTTER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import plotter: {PLOTTER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def font(size: int, *, bold: bool = False, italic: bool = False) -> ImageFont.FreeTypeFont:
    candidates: list[str] = []
    if bold and italic:
        candidates.extend([
            "/System/Library/Fonts/Supplemental/Times New Roman Bold Italic.ttf",
            "/Library/Fonts/Times New Roman Bold Italic.ttf",
        ])
    elif bold:
        candidates.extend([
            "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf",
            "/Library/Fonts/Times New Roman Bold.ttf",
        ])
    elif italic:
        candidates.extend([
            "/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf",
            "/Library/Fonts/Times New Roman Italic.ttf",
        ])
    else:
        candidates.extend([
            "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
            "/Library/Fonts/Times New Roman.ttf",
        ])
    for candidate in candidates:
        if Path(candidate).exists():
            return ImageFont.truetype(candidate, size=size)
    return ImageFont.truetype("/System/Library/Fonts/Times.ttc", size=size)


def text_bbox(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, fnt: ImageFont.ImageFont) -> tuple[int, int, int, int]:
    return tuple(int(v) for v in draw.textbbox(xy, text, font=fnt))


def draw_rich_sphenix(draw: ImageDraw.ImageDraw, xy: tuple[int, int], size: int, fill: tuple[int, int, int]) -> None:
    x, y = xy
    f_sphenix = font(size, bold=True, italic=True)
    f_rest = font(size)
    draw.text((x, y), "sPHENIX", font=f_sphenix, fill=fill)
    draw.text((x + int(draw.textlength("sPHENIX", font=f_sphenix)) + 16, y), "Internal QA", font=f_rest, fill=fill)


def crop_whitespace(img: Image.Image, *, bg=(255, 255, 255), pad: int = 8) -> Image.Image:
    rgb = img.convert("RGB")
    gray = ImageOps.grayscale(ImageChops.difference(rgb, Image.new("RGB", rgb.size, bg)))
    bbox = gray.point(lambda p: 255 if p > 12 else 0).getbbox()
    if not bbox:
        return rgb
    x0, y0, x1, y1 = bbox
    return rgb.crop((max(0, x0 - pad), max(0, y0 - pad), min(rgb.size[0], x1 + pad), min(rgb.size[1], y1 + pad)))


def prepare_reference_panels() -> None:
    REF_DIR.mkdir(parents=True, exist_ok=True)
    if not PPG12_IAN_PAGE19.exists():
        raise FileNotFoundError(f"Render PPG12 IAN page 19 first: {PPG12_IAN_PAGE19}")
    page = Image.open(PPG12_IAN_PAGE19).convert("RGB")
    # These crops keep the actual IAN panels intact while excluding neighboring
    # table-panel labels that appear in the dense multi-panel figure layout.
    crops = {
        PPG12_E11_REF: (250, 560, 970, 1200),
        PPG12_BDT_REF: (1575, 560, 2280, 1200),
    }
    for out, box in crops.items():
        page.crop(box).save(out)


def clean_the42_panel(var: str, out: Path) -> dict:
    plotter = load_plotter()
    table = [t for t in plotter.TABLES if t["slug"] == "fig13_style_pt22_28_cut0"][0]
    xlim, rebin = plotter.ppg12_axis_settings(var)
    files = {
        "data": plotter.open_root(ROOT_DIR / "RecoilJets_pp_ALL_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"),
        "signal_mc": plotter.open_root(ROOT_DIR / "RecoilJets_photonjet5plus10plus20_MERGED.root"),
    }
    cache = plotter.load_inclusive_cache(INCLUSIVE_CACHE, use_stitched_inclusive=False)
    data = plotter.norm_arrays(
        plotter.get_hist(files["data"], "Photon_4_GeV_plus_MBD_NS_geq_1", var, table["pt_token"], table["cut"]),
        rebin,
        xlim,
    )
    sig = plotter.norm_arrays(
        plotter.get_hist(files["signal_mc"], "SIM", var, table["pt_token"], table["cut"]),
        rebin,
        xlim,
    )
    inc = plotter.norm_payload(plotter.cache_hist(cache, "current_ian_jet8to40", var, table["pt_token"], table["cut"]), rebin, xlim)
    npb_scale = plotter.npb_tail_scale(files, table)
    npb = plotter.norm_arrays(
        plotter.get_hist(files["data"], "Photon_4_GeV_plus_MBD_NS_geq_1", var, table["pt_token"], "cut4"),
        rebin,
        xlim,
    )
    npb = plotter.scale_arrays(npb, npb_scale)
    chi2 = plotter.chi2_summary(data, inc)

    fig, ax = plt.subplots(figsize=(7.7, 6.4))
    plotter.draw_shape(ax, data, label="Data", color="black", marker="o", markersize=5.2)
    plotter.draw_shape(ax, sig, label="Signal MC", color="red", linewidth=2.0)
    plotter.draw_shape(ax, inc, label="Inclusive MC", color="blue", linewidth=2.0)
    plotter.draw_shape(ax, npb, label="NPB-tagged data", color="#238b1e", linewidth=2.0)
    plotter.add_panel_annotation(ax, table, chi2, fontsize=14.5, linespacing=1.08, include_chi2=False)
    ax.legend(loc="upper right", frameon=False, fontsize=15, handlelength=1.9)
    ax.set_xlim(*xlim)
    ax.set_ylim(bottom=0)
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax * 1.10)
    ax.tick_params(axis="both", labelsize=15, direction="in", top=True, right=True)
    ax.set_xlabel(var, fontsize=18)
    ax.set_ylabel("normalized counts", fontsize=18)
    fig.subplots_adjust(left=0.14, right=0.96, bottom=0.14, top=0.96)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return {"png": str(out), "var": var, "chi2": chi2, "npb_scale": npb_scale}


def fit_plot(path: Path, size: tuple[int, int], *, crop: bool = True, pad: int = 8) -> Image.Image:
    img = Image.open(path).convert("RGB")
    if crop:
        img = crop_whitespace(img, pad=pad)
    return ImageOps.contain(img, size, method=Image.Resampling.LANCZOS)


def draw_text_box(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    *,
    fill: str,
    outline: str,
    title: str,
    body: str,
    title_color: str,
    body_color: str,
    title_font_px: int = 38,
    body_font_px: int = 37,
) -> list[dict]:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=18, fill=fill, outline=outline, width=2)
    title_f = font(title_font_px, bold=True)
    body_f = font(body_font_px)
    draw.text((x0 + 28, y0 + 22), title, font=title_f, fill=title_color)
    words = body.split()
    lines: list[str] = []
    cur = ""
    max_w = x1 - x0 - 56
    for word in words:
        trial = word if not cur else f"{cur} {word}"
        if draw.textlength(trial, font=body_f) <= max_w:
            cur = trial
        else:
            if cur:
                lines.append(cur)
            cur = word
    if cur:
        lines.append(cur)
    yy = y0 + 74
    for line in lines[:3]:
        draw.text((x0 + 28, yy), line, font=body_f, fill=body_color)
        yy += 42
    return [
        {"kind": "text", "name": title, "role": "audience", "text": title, "bbox": list(text_bbox(draw, (x0 + 28, y0 + 22), title, title_f)), "font_px": title_font_px},
        {"kind": "text", "name": f"{title} body", "role": "audience", "text": body, "bbox": [x0 + 28, y0 + 74, x1 - 28, yy], "font_px": body_font_px},
    ]


def draw_arrow_bullets(
    draw: ImageDraw.ImageDraw,
    *,
    xy: tuple[int, int],
    bullets: list[str],
    max_width: int,
    font_px: int = 40,
) -> list[dict]:
    x, y = xy
    body_f = font(font_px)
    nodes: list[dict] = []
    line_gap = 12
    bullet_gap = 20
    arrow_w = 24
    for idx, text in enumerate(bullets):
        words = text.split()
        lines: list[str] = []
        current = ""
        for word in words:
            trial = word if not current else f"{current} {word}"
            if draw.textlength(trial, font=body_f) <= max_width - 58:
                current = trial
            else:
                if current:
                    lines.append(current)
                current = word
        if current:
            lines.append(current)

        arrow_y = y + 14
        draw.polygon(
            [(x, arrow_y), (x + arrow_w, arrow_y + 12), (x, arrow_y + 24)],
            fill="#2468a8",
        )
        text_x = x + 42
        top_y = y
        for line in lines:
            draw.text((text_x, y), line, font=body_f, fill="#172033")
            y += font_px + line_gap
        bbox = (x, top_y, x + max_width, y - line_gap)
        nodes.append({
            "kind": "text",
            "name": f"arrow bullet {idx + 1}",
            "role": "audience",
            "text": text,
            "bbox": list(bbox),
            "font_px": font_px,
            "title_axis_align": "left",
        })
        y += bullet_gap
    return nodes


def compose_slide(
    *,
    slug: str,
    title: str,
    reference_plot: Path,
    output_plot: Path,
    bullets: list[str],
    speaker_readout: str,
    manifest_extra: dict,
) -> dict:
    slide = Image.new("RGB", (SLIDE_WIDTH_PX, SLIDE_HEIGHT_PX), "white")
    draw = ImageDraw.Draw(slide)
    nodes: list[dict] = []
    title_f = font(68, bold=True)
    title_xy = (82, 52)
    draw.text(title_xy, title, font=title_f, fill="#111827")
    title_box = text_bbox(draw, title_xy, title, title_f)
    nodes.append({"kind": "text", "name": "title", "role": "title", "text": title, "bbox": list(title_box), "font_px": 68, "title_anchor": True})

    nodes.extend(draw_arrow_bullets(draw, xy=(82, 164), bullets=bullets, max_width=2150, font_px=40))

    plot_top, plot_bottom = 330, 1360
    left_frame = (82, plot_top, 1222, plot_bottom)
    right_frame = (1338, plot_top, 2480, plot_bottom)
    for box, label, accent, path, frame_fill, frame_outline, crop_plot in [
        (left_frame, "IAN reference, May 21, 2026", "#525866", reference_plot, "#ffffff", "#c3cad5", True),
        (right_frame, "This analysis base pp output", "#2468a8", output_plot, "#ffffff", "#a8c8e8", False),
    ]:
        x0, y0, x1, y1 = box
        draw.rounded_rectangle(box, radius=16, fill=frame_fill, outline=frame_outline, width=3)
        label_f = font(45, bold=True)
        label_w = min(x1 - x0 - 44, int(draw.textlength(label, font=label_f)) + 64)
        draw.rounded_rectangle((x0 + 22, y0 + 16, x0 + 22 + label_w, y0 + 74), radius=10, fill=accent)
        draw.text((x0 + 50, y0 + 20), label, font=label_f, fill="white")
        fitted = fit_plot(path, (x1 - x0 - 80, y1 - y0 - 116), crop=crop_plot, pad=18)
        px = x0 + (x1 - x0 - fitted.size[0]) // 2
        py = y0 + 94 + (y1 - y0 - 116 - fitted.size[1]) // 2
        slide.paste(fitted, (px, py))
        nodes.append({"kind": "box", "name": f"{label} plot frame", "bbox": [x0, y0, x1, y1], **({"title_axis_align": "left"} if x0 == 82 else {})})

    out_png = OUT_DIR / f"{slug}.png"
    slide.save(out_png)
    layout_path = OUT_DIR / f"{slug}.layout_nodes.json"
    layout = {"title_axis_x": 82, "title_axis_tolerance_px": 12, "nodes": nodes}
    layout_path.write_text(json.dumps(layout, indent=2) + "\n")
    script_path = OUT_DIR / f"{slug}_speaker_script.md"
    script_path.write_text(
        "\n".join([
            "# Speaker Script",
            "",
            f"This slide compares the {PPG12_IAN_LABEL} reference panel to the repaired RecoilJets {manifest_extra['variable_label']} output.",
            bullets[0],
            bullets[1],
            speaker_readout,
            "",
        ])
    )
    manifest_path = OUT_DIR / f"{slug}_manifest.json"
    manifest = {
        "png": str(out_png),
        "layout_nodes": str(layout_path),
        "speaker_script": str(script_path),
        "reference_plot": str(reference_plot),
        "output_plot": str(output_plot),
        "ppg12_ian_page": str(PPG12_IAN_PAGE19),
        "ppg12_ian_label": PPG12_IAN_LABEL,
        "inclusive_cache": str(INCLUSIVE_CACHE),
        "slide_bullets": bullets,
        "speaker_readout": speaker_readout,
        **manifest_extra,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def make_slides() -> list[dict]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    prepare_reference_panels()
    e11_out = OUT_DIR / "the42_ppg12_macro_matched_fig13_e11_to_e33_clean.png"
    bdt_out = OUT_DIR / "the42_ppg12_macro_matched_fig13_bdt_clean.png"
    e11_meta = clean_the42_panel("e11_to_e33", e11_out)
    bdt_meta = clean_the42_panel("bdt", bdt_out)
    manifests = [
        compose_slide(
            slug="the42_ppg12_macro_match_e11_slide",
            title="E11/E33 pp QA against PPG12 IAN v4",
            reference_plot=PPG12_E11_REF,
            output_plot=e11_out,
            bullets=[
                f"Shows the {PPG12_IAN_LABEL} reference and repaired pp output for 22 < E_T < 28 GeV, no NPB cut.",
                "The repaired output recovers the same class hierarchy with NPB control low, inclusive broad, signal high.",
            ],
            speaker_readout="The comparison is meant as a qualitative contract check: the repaired inclusive route now gives the expected PPG12 ordering.",
            manifest_extra={"variable": "e11_to_e33", "variable_label": "E11/E33", "the42_panel": e11_meta},
        ),
        compose_slide(
            slug="the42_ppg12_macro_match_bdt_slide",
            title="BDT-score pp QA against PPG12 IAN v4",
            reference_plot=PPG12_BDT_REF,
            output_plot=bdt_out,
            bullets=[
                f"Shows the {PPG12_IAN_LABEL} reference and repaired pp output for the same Figure 13 BDT-score panel.",
                "The score hierarchy matches with NPB-tagged data and inclusive MC low, while signal MC concentrates high.",
            ],
            speaker_readout="The readout is the same ordering and shape logic as PPG12, not a claim of identical bin-by-bin yields.",
            manifest_extra={"variable": "bdt", "variable_label": "BDT-score", "the42_panel": bdt_meta},
        ),
    ]
    return manifests


if __name__ == "__main__":
    for manifest in make_slides():
        print(manifest["png"])

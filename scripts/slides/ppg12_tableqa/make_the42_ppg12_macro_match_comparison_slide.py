#!/usr/bin/env python3
"""Build a THE-42 slide comparing PPG12 reference panels to macro-matched output."""

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
REPAIR_DIR = CAMPAIGN / "repair_diagnostics"
ROOT_DIR = CAMPAIGN / "merged_roots"
INCLUSIVE_CACHE = CAMPAIGN / "inclusive_sample_hist_cache/the42_ppg12_tableqa_v1_inclusive_sample_projectx_hists.json"
PLOTTER_PATH = REPO_ROOT / "scripts/plotting/pp_currentian/make_the42_ppg12_tableqa_v1_tables.py"

PPG12_E11_REF = Path(
    "/var/folders/l3/f02nw86n5cn0tpf_zstf0ypr0000gn/T/TemporaryItems/"
    "NSIRD_screencaptureui_XLagwx/Screenshot 2026-06-12 at 6.58.11\u202fPM.png"
)
PPG12_BDT_REF = Path(
    "/var/folders/l3/f02nw86n5cn0tpf_zstf0ypr0000gn/T/TemporaryItems/"
    "NSIRD_screencaptureui_fhQLdU/Screenshot 2026-06-11 at 12.09.40\u202fPM.png"
)
THE42_E11 = REPAIR_DIR / "the42_ppg12_macro_matched_fig13_e11_to_e33.png"


def load_plotter():
    spec = importlib.util.spec_from_file_location("tableqa_plotter", PLOTTER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import plotter: {PLOTTER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def font(size: int, *, bold: bool = False, italic: bool = False) -> ImageFont.FreeTypeFont:
    candidates = []
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
    box = draw.textbbox(xy, text, font=fnt)
    return tuple(int(v) for v in box)


def draw_rich_sphenix(draw: ImageDraw.ImageDraw, xy: tuple[int, int], size: int, fill: tuple[int, int, int]) -> tuple[int, int, int, int]:
    x, y = xy
    f_sphenix = font(size, bold=True, italic=True)
    f_rest = font(size)
    draw.text((x, y), "sPHENIX", font=f_sphenix, fill=fill)
    w = int(draw.textlength("sPHENIX", font=f_sphenix))
    draw.text((x + w + 16, y), "Internal QA", font=f_rest, fill=fill)
    return text_bbox(draw, xy, "sPHENIX Internal QA", f_rest)


def crop_whitespace(img: Image.Image, *, bg=(255, 255, 255), pad: int = 12) -> Image.Image:
    rgb = img.convert("RGB")
    bg_img = Image.new("RGB", rgb.size, bg)
    # Use a fast approximate crop based on non-white pixels.
    gray = ImageOps.grayscale(ImageChops.difference(rgb, bg_img))
    bbox = gray.point(lambda p: 255 if p > 12 else 0).getbbox()
    if not bbox:
        return rgb
    x0, y0, x1, y1 = bbox
    x0 = max(0, x0 - pad)
    y0 = max(0, y0 - pad)
    x1 = min(rgb.size[0], x1 + pad)
    y1 = min(rgb.size[1], y1 + pad)
    return rgb.crop((x0, y0, x1, y1))


def fit_image(path: Path, size: tuple[int, int], *, crop: bool = True, crop_top_px: int = 0) -> Image.Image:
    img = Image.open(path).convert("RGB")
    if crop:
        img = crop_whitespace(img)
    if crop_top_px > 0 and crop_top_px < img.size[1] - 20:
        img = img.crop((0, crop_top_px, img.size[0], img.size[1]))
    fitted = ImageOps.contain(img, size, method=Image.Resampling.LANCZOS)
    canvas = Image.new("RGB", size, "white")
    canvas.paste(fitted, ((size[0] - fitted.size[0]) // 2, (size[1] - fitted.size[1]) // 2))
    return canvas


def contain_image(path: Path, size: tuple[int, int], *, crop: bool = True, crop_top_px: int = 0) -> Image.Image:
    img = Image.open(path).convert("RGB")
    if crop:
        img = crop_whitespace(img)
    if crop_top_px > 0 and crop_top_px < img.size[1] - 20:
        img = img.crop((0, crop_top_px, img.size[0], img.size[1]))
    return ImageOps.contain(img, size, method=Image.Resampling.LANCZOS)


def generate_bdt_panel() -> Path:
    plotter = load_plotter()
    table = [t for t in plotter.TABLES if t["slug"] == "fig13_style_pt22_28_cut0"][0]
    var = "bdt"
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
    sig = plotter.norm_arrays(plotter.get_hist(files["signal_mc"], "SIM", var, table["pt_token"], table["cut"]), rebin, xlim)
    inc_payload = plotter.cache_hist(cache, "current_ian_jet8to40", var, table["pt_token"], table["cut"])
    inc = plotter.norm_payload(inc_payload, rebin, xlim)
    npb_scale = plotter.npb_tail_scale(files, table)
    npb = plotter.norm_arrays(
        plotter.get_hist(files["data"], "Photon_4_GeV_plus_MBD_NS_geq_1", var, table["pt_token"], "cut4"),
        rebin,
        xlim,
    )
    npb = plotter.scale_arrays(npb, npb_scale)
    chi2 = plotter.chi2_summary(data, inc)

    fig = plt.figure(figsize=(5.0, 6.25))
    sub = fig.add_gridspec(2, 1, height_ratios=[4.0, 1.15], hspace=0.03, left=0.16, right=0.96, bottom=0.12, top=0.95)
    ax = fig.add_subplot(sub[0, 0])
    rax = fig.add_subplot(sub[1, 0], sharex=ax)
    plotter.draw_shape(ax, data, label="Data", color="black", marker="o", markersize=4.0)
    plotter.draw_shape(ax, sig, label="Signal MC", color="red", linewidth=1.6)
    plotter.draw_shape(ax, inc, label="Inclusive MC", color="blue", linewidth=1.6)
    plotter.draw_shape(ax, npb, label="NPB-tagged data", color="#238b1e", linewidth=1.6)
    plotter.add_panel_annotation(ax, table, chi2)
    ax.legend(loc="upper right", frameon=False, fontsize=10, handlelength=1.8)
    ax.set_xlim(*xlim)
    ax.set_ylim(bottom=0)
    ax.tick_params(axis="both", labelsize=10, direction="in", top=True, right=True, labelbottom=False)
    ax.set_ylabel("normalized counts", fontsize=12)
    plotter.draw_residual(rax, data, inc, xlim)
    rax.set_xlabel("bdt", fontsize=11)
    rax.tick_params(axis="both", labelsize=10, direction="in", top=True, right=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "the42_ppg12_macro_matched_fig13_bdt.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)

    manifest = {
        "png": str(out),
        "source": "macro-matched THE-42 baseV3E table-QA using projectx inclusive sample cache",
        "inclusive_cache": str(INCLUSIVE_CACHE),
        "npb_scale": npb_scale,
        "chi2": chi2,
    }
    (OUT_DIR / "the42_ppg12_macro_matched_fig13_bdt_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return out


def draw_round_rect(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], fill: str, outline: str, width: int = 2, radius: int = 18) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def draw_panel(
    slide: Image.Image,
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    image_path: Path,
    label: str,
    tag: str,
    *,
    accent: str,
    crop: bool = True,
    crop_top_px: int = 0,
) -> dict:
    x0, y0, x1, y1 = box
    draw_round_rect(draw, box, fill="#ffffff", outline="#d8dde6", width=2, radius=12)
    inner = (x0 + 10, y0 + 10, x1 - 10, y1 - 10)
    fitted = fit_image(image_path, (inner[2] - inner[0], inner[3] - inner[1]), crop=crop, crop_top_px=crop_top_px)
    slide.paste(fitted, (inner[0], inner[1]))
    tag_font = font(29, bold=True)
    label_font = font(29, bold=True)
    chip_y0 = y0 + 12
    draw.rounded_rectangle((x0 + 16, chip_y0, x0 + 122, chip_y0 + 42), radius=9, fill=accent)
    draw.text((x0 + 34, chip_y0 + 4), tag, font=tag_font, fill="white")
    draw.text((x0 + 138, chip_y0 + 4), label, font=label_font, fill="#172033")
    return {
        "kind": "box",
        "name": f"panel {tag} {label}",
        "bbox": [x0, y0, x1, y1],
        **({"title_axis_align": "left"} if x0 == 82 else {}),
    }


def draw_plot_column(
    slide: Image.Image,
    draw: ImageDraw.ImageDraw,
    *,
    x0: int,
    y0: int,
    col_w: int,
    max_plot_h: int,
    image_path: Path,
    label: str,
    tag: str,
    accent: str,
) -> dict:
    chip_font = font(27, bold=True)
    label_font = font(29, bold=True)
    chip = (x0, y0, x0 + 92, y0 + 40)
    draw.rounded_rectangle(chip, radius=9, fill=accent)
    draw.text((x0 + 17, y0 + 3), tag, font=chip_font, fill="white")
    draw.text((x0 + 108, y0 + 4), label, font=label_font, fill="#172033")

    plot_top = y0 + 52
    fitted = contain_image(image_path, (col_w, max_plot_h - 52), crop=True)
    px = x0 + (col_w - fitted.size[0]) // 2
    py = plot_top
    frame = (px - 8, py - 8, px + fitted.size[0] + 8, py + fitted.size[1] + 8)
    draw.rounded_rectangle(frame, radius=10, fill="#ffffff", outline="#d8dde6", width=2)
    slide.paste(fitted, (px, py))
    return {
        "kind": "box",
        "name": f"plot column {tag} {label}",
        "bbox": [frame[0], frame[1], frame[2], frame[3]],
    }


def draw_bullet(draw: ImageDraw.ImageDraw, x: int, y: int, text: str, width: int, *, fill="#172033") -> tuple[int, int, int, int]:
    f = font(38)
    bullet_f = font(42, bold=True)
    words = text.split()
    lines: list[str] = []
    cur = ""
    for word in words:
        trial = word if not cur else f"{cur} {word}"
        if draw.textlength(trial, font=f) <= width:
            cur = trial
        else:
            if cur:
                lines.append(cur)
            cur = word
    if cur:
        lines.append(cur)
    draw.text((x, y + 1), "•", font=bullet_f, fill="#2468a8")
    yy = y
    for line in lines:
        draw.text((x + 34, yy), line, font=f, fill=fill)
        yy += 47
    return (x, y, x + width + 34, yy)


def make_slide() -> dict:
    bdt_panel = generate_bdt_panel()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    slide = Image.new("RGB", (SLIDE_WIDTH_PX, SLIDE_HEIGHT_PX), "#f7f9fc")
    draw = ImageDraw.Draw(slide)

    title_font = font(68, bold=True)
    title = "THE-42 pp QA recovers PPG12 class ordering in E11/E33 and BDT"
    title_x, title_y = 82, 58
    draw.text((title_x, title_y), title, font=title_font, fill="#111827")
    title_box = text_bbox(draw, (title_x, title_y), title, title_font)
    draw.rounded_rectangle((82, 132, 2480, 143), radius=5, fill="#2468a8")
    draw_rich_sphenix(draw, (2125, 54), 30, (35, 43, 58))

    left_x0, col_w, gutter = 82, 830, 24
    row_h = 512
    row1_y, row2_y = 195, 745
    right_x0, right_x1 = 1810, 2480

    nodes = [
        {
            "kind": "text",
            "name": "title",
            "role": "title",
            "text": title,
            "bbox": list(title_box),
            "font_px": 68,
            "title_anchor": True,
        },
    ]

    grid_x0, grid_x1 = 28, 2532
    grid_y0, grid_y1 = 154, 1410
    col_gap = 14
    col_w = (grid_x1 - grid_x0 - 3 * col_gap) // 4
    col_x = [grid_x0 + i * (col_w + col_gap) for i in range(4)]
    max_plot_h = grid_y1 - grid_y0
    panel_nodes = [
        draw_plot_column(slide, draw, x0=col_x[0], y0=grid_y0, col_w=col_w, max_plot_h=max_plot_h, image_path=PPG12_E11_REF, label="PPG12 E11/E33", tag="REF", accent="#525866"),
        draw_plot_column(slide, draw, x0=col_x[1], y0=grid_y0, col_w=col_w, max_plot_h=max_plot_h, image_path=THE42_E11, label="THE-42 E11/E33", tag="OUT", accent="#2468a8"),
        draw_plot_column(slide, draw, x0=col_x[2], y0=grid_y0, col_w=col_w, max_plot_h=max_plot_h, image_path=PPG12_BDT_REF, label="PPG12 BDT", tag="REF", accent="#525866"),
        draw_plot_column(slide, draw, x0=col_x[3], y0=grid_y0, col_w=col_w, max_plot_h=max_plot_h, image_path=bdt_panel, label="THE-42 BDT", tag="OUT", accent="#2468a8"),
    ]
    nodes.extend(panel_nodes)

    readout = (82, 1138, 2480, 1244)
    draw.rounded_rectangle(readout, radius=12, fill="#eef6ff", outline="#b9d2ee", width=2)
    readout_font = font(40, bold=True)
    readout_text = "Sample-level jet8-40 cache passes ordering check; NPB is a control overlay; do not compare bin-by-bin."
    readout_xy = (readout[0] + 26, readout[1] + 28)
    draw.text(readout_xy, readout_text, font=readout_font, fill="#123c69")
    nodes.append({
        "kind": "text",
        "name": "bottom interpretation readout",
        "role": "audience",
        "text": readout_text,
        "bbox": list(text_bbox(draw, readout_xy, readout_text, readout_font)),
        "font_px": 40,
        "title_band_exception": True,
    })

    out_png = OUT_DIR / "the42_ppg12_macro_match_comparison_slide.png"
    slide.save(out_png)

    layout = {
        "title_axis_x": title_x,
        "title_axis_tolerance_px": 12,
        "nodes": nodes,
    }
    layout_path = OUT_DIR / "the42_ppg12_macro_match_comparison_slide.layout_nodes.json"
    layout_path.write_text(json.dumps(layout, indent=2) + "\n")

    script_path = OUT_DIR / "the42_ppg12_macro_match_comparison_slide_speaker_script.md"
    script_path.write_text(
        "\n".join([
            "# Speaker Script",
            "",
            "This slide shows the specific check that resolved the confusing pp table-QA behavior.",
            "Across E11 over E33 and BDT score, the THE-42 output recovers the same class ordering seen in the PPG12 reference once inclusive MC is plotted from the sample-level cache instead of the double-weighted final stitched ROOT.",
            "The NPB-tagged curve is shown as a control overlay here, not as a rejection applied to the plotted candidates.",
            "The point is not bin-by-bin equality. PPG12 used the base_v1E configuration, while this checkpoint is our baseV3E RecoilJets table-QA output. The important result is that the plotting and weighting contract is coherent and the expected ordering is visible.",
            "",
        ])
    )

    manifest_path = OUT_DIR / "the42_ppg12_macro_match_comparison_slide_manifest.json"
    manifest = {
        "png": str(out_png),
        "layout_nodes": str(layout_path),
        "speaker_script": str(script_path),
        "inputs": {
            "ppg12_e11_reference": str(PPG12_E11_REF),
            "ppg12_bdt_reference": str(PPG12_BDT_REF),
            "the42_e11_macro_matched": str(THE42_E11),
            "the42_bdt_macro_matched": str(bdt_panel),
            "inclusive_cache": str(INCLUSIVE_CACHE),
        },
        "claim": "The PPG12-style hierarchy is qualitatively recovered after avoiding inclusive finalStitch double weighting and matching the PPG12 plotting contract.",
        "caveat": "Qualitative hierarchy check, not bin-by-bin equality; PPG12 reference uses base_v1E while THE-42 uses baseV3E/RecoilJets table-QA.",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    manifest = make_slide()
    print(manifest["png"])

#!/usr/bin/env python3
"""Build pp PPG12/baseV3E stage-flow slide candidates.

These slides show one variable across the table-QA stages in a single 1x3 row:
all candidates, after preselection/NPB, and after tight ID.  The plot canvases
are regenerated for the slide geometry instead of pasted from the full table.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from PIL import Image, ImageChops, ImageDraw, ImageFont, ImageOps


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = next(p for p in SCRIPT_PATH.parents if (p / "AGENTS.md").exists())
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.append(str(SCRIPTS_DIR))

from slides.common.slide_defaults import SLIDE_HEIGHT_PX, SLIDE_WIDTH_PX  # noqa: E402


CAMPAIGN = REPO_ROOT / "dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_basev3e_20260611"
OUT_DIR = CAMPAIGN / "stage_flow_slide"
ROOT_DIR = CAMPAIGN / "merged_roots"
INCLUSIVE_CACHE = CAMPAIGN / "inclusive_sample_hist_cache/the42_ppg12_tableqa_v1_inclusive_sample_projectx_hists.json"
PLOTTER_PATH = REPO_ROOT / "scripts/plotting/pp_currentian/make_the42_ppg12_tableqa_v1_tables.py"

DATA_ROOT = ROOT_DIR / "RecoilJets_pp_ALL_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
SIGNAL_ROOT = ROOT_DIR / "RecoilJets_photonjet5plus10plus20_MERGED.root"
INCLUSIVE_STITCHED_ROOT = ROOT_DIR / "RecoilJets_jet5plus8plus12plus20plus30plus40_MERGED.root"
TOPDIR_DATA = "Photon_4_GeV_plus_MBD_NS_geq_1"
PT_TOKEN = "1535"
PT_LABEL = r"$15<E_T<35$ GeV"
PT_LABEL_TEXT = "15 < E_T < 35 GeV"
MODEL_LABEL = "PPG12/baseV3E pp output"
IAN_LABEL = "PPG12 IAN v4 (May 21, 2026)"

STAGES = [
    ("cut0", "Before preselection", "all candidates"),
    ("cut1", "After preselection", "preselection + NPB cut"),
    ("cut2", "After tight ID", "tight BDT ID"),
]

VARS = {
    "bdt": {
        "axis": "BDT score",
        "title": "BDT score stage flow is coherent in the pp checkpoint",
        "slug": "bdt",
        "takeaways": [
            "Read left to right - all candidates, preselection plus NPB cut, then tight photon ID.",
            "Data moves toward the signal-MC score shape after tight ID; this is PPG12/baseV3E, not the separately trained pp BDT.",
        ],
    },
    "e11_to_e33": {
        "axis": r"$E_{11}/E_{33}$",
        "title": "E11/E33 sharpens through the same pp selection chain",
        "slug": "e11_to_e33",
        "takeaways": [
            "Read left to right - all candidates, preselection plus NPB cut, then tight photon ID.",
            "The tight-ID selection sharpens E11/E33 and leaves data closer to signal MC while inclusive MC stays broader.",
        ],
    },
}


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


def wrap_text(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.ImageFont, max_width: int) -> list[str]:
    lines: list[str] = []
    current = ""
    for word in text.split():
        trial = word if not current else f"{current} {word}"
        if draw.textlength(trial, font=fnt) <= max_width:
            current = trial
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    return lines


def draw_arrow_bullets(
    draw: ImageDraw.ImageDraw,
    *,
    xy: tuple[int, int],
    bullets: list[str],
    max_width: int,
    font_px: int,
) -> tuple[int, list[dict]]:
    x, y = xy
    fnt = font(font_px)
    nodes: list[dict] = []
    for idx, text in enumerate(bullets):
        top = y
        arrow_y = y + 11
        draw.polygon([(x, arrow_y), (x + 24, arrow_y + 12), (x, arrow_y + 24)], fill="#2468a8")
        lines = wrap_text(draw, text, fnt, max_width - 54)
        for line in lines:
            draw.text((x + 44, y), line, font=fnt, fill="#172033")
            y += font_px + 8
        nodes.append({
            "kind": "text",
            "name": f"arrow bullet {idx + 1}",
            "role": "audience",
            "text": text,
            "bbox": [x, top, x + max_width, y - 8],
            "font_px": font_px,
            "title_axis_align": "left",
        })
        y += 12
    return y, nodes


def stage_table(stage: str) -> dict:
    labels = {cut: (title, detail) for cut, title, detail in STAGES}
    title, detail = labels[stage]
    return {
        "pt_token": PT_TOKEN,
        "pt_label": PT_LABEL,
        "cut": stage,
        "cut_label": detail,
        "stage_title": title,
        "include_npb_template": stage == "cut0",
    }


def panel_arrays(plotter, files: dict, inclusive_cache: dict, var: str, stage: str) -> dict:
    xlim, rebin = plotter.ppg12_axis_settings(var)
    data = plotter.norm_arrays(plotter.get_hist(files["data"], TOPDIR_DATA, var, PT_TOKEN, stage), rebin, xlim)
    sig = plotter.norm_arrays(plotter.get_hist(files["signal_mc"], "SIM", var, PT_TOKEN, stage), rebin, xlim)
    inc_payload = plotter.cache_hist(inclusive_cache, "current_ian_jet8to40", var, PT_TOKEN, stage)
    inc = plotter.norm_payload(inc_payload, rebin, xlim)
    npb = None
    npb_scale = None
    if stage == "cut0":
        npb_scale = plotter.npb_tail_scale(files, stage_table(stage))
        npb = plotter.norm_arrays(plotter.get_hist(files["data"], TOPDIR_DATA, var, PT_TOKEN, "cut4"), rebin, xlim)
        npb = plotter.scale_arrays(npb, npb_scale)
    return {
        "xlim": xlim,
        "rebin": rebin,
        "data": data,
        "signal_mc": sig,
        "inclusive_mc": inc,
        "npb_tagged_data": npb,
        "npb_scale": npb_scale,
        "inclusive_stats": plotter.payload_stats_after_transform(inc_payload, rebin, xlim),
    }


def draw_stage_plot(plotter, files: dict, inclusive_cache: dict, var: str, stage: str, out: Path) -> dict:
    info = panel_arrays(plotter, files, inclusive_cache, var, stage)
    fig, ax = plt.subplots(figsize=(8.0, 8.5))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    plotter.draw_shape(ax, info["data"], label="Data", color="black", marker="o", markersize=5.8)
    plotter.draw_shape(ax, info["signal_mc"], label="Signal MC", color="red", linewidth=2.3)
    plotter.draw_shape(ax, info["inclusive_mc"], label="Inclusive MC", color="blue", linewidth=2.3)
    if info["npb_tagged_data"] is not None:
        plotter.draw_shape(ax, info["npb_tagged_data"], label="NPB-tagged data", color="#238b1e", linewidth=2.2)

    table = stage_table(stage)
    plotter.add_panel_annotation(ax, table, None, fontsize=17.5, linespacing=1.08, include_chi2=False)
    ax.legend(loc="upper right", frameon=False, fontsize=15.5, handlelength=1.8)
    ax.set_xlim(*info["xlim"])
    ax.set_ylim(bottom=0)
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax * 1.12)
    ax.set_xlabel(VARS[var]["axis"], fontsize=21)
    ax.set_ylabel("normalized counts", fontsize=21)
    ax.tick_params(axis="both", labelsize=17, direction="in", top=True, right=True)
    fig.subplots_adjust(left=0.15, right=0.965, bottom=0.13, top=0.965)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)

    entries = {}
    for key in ("data", "signal_mc", "inclusive_mc", "npb_tagged_data"):
        arrays = info[key]
        entries[key] = float(arrays[1].sum()) if arrays is not None else 0.0
    return {
        "png": str(out),
        "stage": stage,
        "stage_label": table["stage_title"],
        "xlim": list(info["xlim"]),
        "rebin": info["rebin"],
        "normalized_integrals_visible": entries,
        "npb_scale": info["npb_scale"],
        "inclusive_stats": info["inclusive_stats"],
    }


def crop_whitespace(img: Image.Image, *, bg=(255, 255, 255), pad: int = 4) -> Image.Image:
    rgb = img.convert("RGB")
    diff = ImageChops.difference(rgb, Image.new("RGB", rgb.size, bg))
    gray = ImageOps.grayscale(diff)
    bbox = gray.point(lambda p: 255 if p > 10 else 0).getbbox()
    if not bbox:
        return rgb
    x0, y0, x1, y1 = bbox
    return rgb.crop((max(0, x0 - pad), max(0, y0 - pad), min(rgb.size[0], x1 + pad), min(rgb.size[1], y1 + pad)))


def fit_plot(path: Path, size: tuple[int, int]) -> Image.Image:
    img = Image.open(path).convert("RGB")
    img = crop_whitespace(img, pad=10)
    return ImageOps.contain(img, size, method=Image.Resampling.LANCZOS)


def compose_slide(var: str, panel_meta: list[dict]) -> dict:
    spec = VARS[var]
    slide = Image.new("RGB", (SLIDE_WIDTH_PX, SLIDE_HEIGHT_PX), "white")
    draw = ImageDraw.Draw(slide)
    nodes: list[dict] = []

    title_f = font(74, bold=True)
    title_xy = (82, 48)
    draw.text(title_xy, spec["title"], font=title_f, fill="#111827")
    title_box = text_bbox(draw, title_xy, spec["title"], title_f)
    nodes.append({"kind": "text", "name": "title", "role": "title", "text": spec["title"], "bbox": list(title_box), "font_px": 74, "title_anchor": True})

    bullet_bottom, bullet_nodes = draw_arrow_bullets(
        draw,
        xy=(82, 162),
        bullets=spec["takeaways"],
        max_width=2280,
        font_px=43,
    )
    nodes.extend(bullet_nodes)

    panel_top = max(328, bullet_bottom + 20)
    panel_bottom = 1384
    left = 82
    gutter = 34
    panel_w = (SLIDE_WIDTH_PX - 2 * left - 2 * gutter) // 3
    label_f = font(38, bold=True)
    small_f = font(26)
    for idx, meta in enumerate(panel_meta):
        x0 = left + idx * (panel_w + gutter)
        y0 = panel_top
        x1 = x0 + panel_w
        y1 = panel_bottom
        draw.rounded_rectangle((x0, y0, x1, y1), radius=14, fill="#ffffff", outline="#d0d6df", width=2)
        stage_label = meta["stage_label"]
        draw.text((x0 + 18, y0 + 15), stage_label, font=label_f, fill="#173b63")
        draw.text((x0 + 18, y0 + 55), f"{MODEL_LABEL}, 15 < ET < 35 GeV", font=small_f, fill="#4b5563")
        fitted = fit_plot(Path(meta["png"]), (panel_w - 20, y1 - y0 - 82))
        px = x0 + (panel_w - fitted.size[0]) // 2
        py = y0 + 76 + (y1 - y0 - 82 - fitted.size[1]) // 2
        py = min(py + 44, y1 - fitted.size[1] - 18)
        slide.paste(fitted, (px, py))
        nodes.append({"kind": "box", "name": f"{stage_label} plot frame", "bbox": [x0, y0, x1, y1], **({"title_axis_align": "left"} if idx == 0 else {})})
        nodes.append({"kind": "text", "name": f"{stage_label} panel label", "role": "audience", "text": stage_label, "bbox": list(text_bbox(draw, (x0 + 18, y0 + 15), stage_label, label_f)), "font_px": 38})

    slug = f"ppg12_basev3e_pp_stage_flow_{spec['slug']}_slide"
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
            f"This slide shows {spec['axis']} through the pp table-QA selection flow.",
            f"The inputs are the repaired {MODEL_LABEL} ROOTs, with inclusive MC drawn from the sample-level no-double-weight cache.",
            "The point is not to claim this is the separately trained pp BDT. It is the coherent PPG12/baseV3E infrastructure checkpoint that the AuAu table-QA output follows.",
            "",
        ])
    )
    manifest_path = OUT_DIR / f"{slug}_manifest.json"
    manifest = {
        "png": str(out_png),
        "layout_nodes": str(layout_path),
        "speaker_script": str(script_path),
        "variable": var,
        "axis_label": spec["axis"],
        "model_label": MODEL_LABEL,
        "ian_reference": IAN_LABEL,
        "pt_token": PT_TOKEN,
        "pt_label": PT_LABEL_TEXT,
        "stages": panel_meta,
        "inputs": {
            "data": str(DATA_ROOT),
            "signal_mc": str(SIGNAL_ROOT),
            "inclusive_cache": str(INCLUSIVE_CACHE),
            "stitched_inclusive_not_used_for_plot": str(INCLUSIVE_STITCHED_ROOT),
        },
        "contract_note": (
            "Current slide uses stock PPG12/baseV3E score path. "
            "A separately trained pp BDT view requires a model export/rescore and is not represented by these ROOT histograms."
        ),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def make_slides() -> list[dict]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plotter = load_plotter()
    files = {
        "data": plotter.open_root(DATA_ROOT),
        "signal_mc": plotter.open_root(SIGNAL_ROOT),
    }
    inclusive_cache = plotter.load_inclusive_cache(INCLUSIVE_CACHE, use_stitched_inclusive=False)
    manifests: list[dict] = []
    for var in ("bdt", "e11_to_e33"):
        panel_meta = []
        for stage, _, _ in STAGES:
            panel_png = OUT_DIR / "panels" / f"ppg12_basev3e_pp_{var}_{PT_TOKEN}_{stage}.png"
            panel_meta.append(draw_stage_plot(plotter, files, inclusive_cache, var, stage, panel_png))
        manifests.append(compose_slide(var, panel_meta))
    combined = OUT_DIR / "ppg12_basev3e_pp_stage_flow_manifest.json"
    combined.write_text(json.dumps({"slides": manifests}, indent=2) + "\n")
    return manifests


if __name__ == "__main__":
    for manifest in make_slides():
        print(manifest["png"])

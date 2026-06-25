#!/usr/bin/env python3
"""Build an interim AuAu/pp E11/E33 stage-flow slide candidate."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageChops, ImageDraw, ImageFont, ImageOps


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = next(p for p in SCRIPT_PATH.parents if (p / "AGENTS.md").exists())
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.append(str(SCRIPTS_DIR))

from slides.common.slide_defaults import SLIDE_HEIGHT_PX, SLIDE_WIDTH_PX  # noqa: E402


AU_AU_CAMPAIGN = REPO_ROOT / "dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_auauDataNewV008_b002_20260612"
AU_AU_CACHE = AU_AU_CAMPAIGN / "interim_complete_runs/the42_b002_interim_complete_run_shower_shape_hists.json"
PP_CAMPAIGN = REPO_ROOT / "dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_basev3e_20260611"
PP_ROOT_DIR = PP_CAMPAIGN / "merged_roots"
PP_INCLUSIVE_CACHE = PP_CAMPAIGN / "inclusive_sample_hist_cache/the42_ppg12_tableqa_v1_inclusive_sample_projectx_hists.json"
PLOTTER_PATH = REPO_ROOT / "scripts/plotting/pp_currentian/make_the42_ppg12_tableqa_v1_tables.py"
OUT_DIR = AU_AU_CAMPAIGN / "interim_complete_runs/e11_pp_reference_slide"

PP_DATA_ROOT = PP_ROOT_DIR / "RecoilJets_pp_ALL_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
PP_SIGNAL_ROOT = PP_ROOT_DIR / "RecoilJets_photonjet5plus10plus20_MERGED.root"
PP_STITCHED_INCLUSIVE_ROOT = PP_ROOT_DIR / "RecoilJets_jet5plus8plus12plus20plus30plus40_MERGED.root"
PP_TOPDIR_DATA = "Photon_4_GeV_plus_MBD_NS_geq_1"

VAR = "e11_to_e33"
PT_TOKEN = "1535"
PT_LABEL = "15 < E_T < 35 GeV"
STAGES = [
    ("cut0", "Before preselection"),
    ("cut1", "After preselection"),
    ("cut2", "After tight ID"),
]
CENTRALITIES = [
    ("cent0_20", "AuAu 0-20%", "#111111"),
    ("cent20_50", "AuAu 20-50%", "#1f77b4"),
    ("cent50_80", "AuAu 50-80%", "#d95f02"),
]


def font(size: int, *, bold: bool = False, italic: bool = False) -> ImageFont.FreeTypeFont:
    names: list[str]
    if bold and italic:
        names = ["Times New Roman Bold Italic.ttf", "Times New Roman Bold Italic.ttf"]
    elif bold:
        names = ["Times New Roman Bold.ttf", "Times New Roman Bold.ttf"]
    elif italic:
        names = ["Times New Roman Italic.ttf", "Times New Roman Italic.ttf"]
    else:
        names = ["Times New Roman.ttf", "Times New Roman.ttf"]
    roots = [Path("/System/Library/Fonts/Supplemental"), Path("/Library/Fonts")]
    for root in roots:
        for name in names:
            path = root / name
            if path.exists():
                return ImageFont.truetype(str(path), size=size)
    return ImageFont.truetype("/System/Library/Fonts/Times.ttc", size=size)


def text_bbox(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, fnt: ImageFont.ImageFont) -> list[int]:
    return [int(v) for v in draw.textbbox(xy, text, font=fnt)]


def load_plotter():
    spec = importlib.util.spec_from_file_location("tableqa_plotter", PLOTTER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import pp plotter from {PLOTTER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_auau_cache() -> dict:
    payload = json.loads(AU_AU_CACHE.read_text())
    if payload.get("schema") != "THE42_B002_INTERIM_COMPLETE_RUN_SHOWER_SHAPES_V1":
        raise RuntimeError(f"Unexpected AuAu cache schema in {AU_AU_CACHE}")
    return payload


def rebin(x: np.ndarray, y: np.ndarray, e: np.ndarray, factor: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if factor <= 1 or len(x) < factor:
        return x, y, e
    n = len(x) // factor
    trim = n * factor
    xr = x[:trim].reshape(n, factor)
    yr = y[:trim].reshape(n, factor)
    er = e[:trim].reshape(n, factor)
    return xr.mean(axis=1), yr.sum(axis=1), np.sqrt(np.sum(er * er, axis=1))


def normalize_visible(
    x: np.ndarray,
    y: np.ndarray,
    e: np.ndarray,
    xlim: tuple[float, float] = (0.0, 1.0),
    factor: int = 4,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x, y, e = rebin(x, y, e, factor)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (x >= xlim[0]) & (x <= xlim[1])
    x, y, e = x[mask], y[mask], e[mask]
    total = float(np.sum(y))
    if total > 0:
        y = y / total
        e = e / total
    return x, y, e


def auau_arrays(cache: dict, cent: str, stage: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    item = cache["hists"][VAR][cent][stage]
    x = np.asarray(item["x"], dtype=float)
    y = np.asarray(item["y"], dtype=float)
    e = np.asarray(item["e"], dtype=float)
    return normalize_visible(x, y, e, factor=4)


def add_sphenix_label(ax, *, collision: str, include_pt: bool = False) -> None:
    ax.text(
        0.05,
        0.92,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.8,
    )
    detail = collision
    if include_pt:
        detail += f"\n{PT_LABEL}, |eta| < 0.7"
    ax.text(
        0.05,
        0.78,
        detail,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.0,
        linespacing=1.05,
    )


def draw_auau_panel(cache: dict, cent: str, color: str, stage: str, out: Path) -> dict:
    x, y, e = auau_arrays(cache, cent, stage)
    fig, ax = plt.subplots(figsize=(6.9, 2.25), dpi=220)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.step(x, y, where="mid", color=color, linewidth=1.8)
    ax.errorbar(x, y, yerr=e, fmt="o", color=color, markersize=2.6, elinewidth=0.65, capsize=0)
    ax.set_xlim(0, 1)
    ymax = float(np.nanmax(y)) if y.size else 1.0
    ax.set_ylim(0, max(0.02, ymax * 1.28))
    ax.grid(axis="y", color="#e6e9ef", linewidth=0.55)
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=True, labelsize=8.4, pad=1)
    ax.set_xlabel(r"$E_{11}/E_{33}$", fontsize=10.0, labelpad=1)
    ax.set_ylabel("norm.", fontsize=9.2, labelpad=1)
    if cent == "cent0_20" and stage == "cut0":
        add_sphenix_label(ax, collision="Au+Au sqrt(sNN)=200 GeV", include_pt=True)
    fig.subplots_adjust(left=0.075, right=0.992, top=0.985, bottom=0.205)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return {"png": str(out), "entries": float(np.sum(y)), "max_norm_bin": ymax}


def pp_stage_table(stage: str) -> dict:
    labels = {
        "cut0": "all candidates",
        "cut1": "preselection + NPB cut",
        "cut2": "tight BDT ID",
    }
    return {
        "pt_token": PT_TOKEN,
        "pt_label": r"$15<E_T<35$ GeV",
        "cut": stage,
        "cut_label": labels[stage],
        "stage_title": dict(STAGES)[stage],
        "include_npb_template": stage == "cut0",
    }


def draw_pp_panel(plotter, files: dict, inclusive_cache: dict, stage: str, out: Path) -> dict:
    xlim, rebin_factor = plotter.ppg12_axis_settings(VAR)
    data = plotter.norm_arrays(plotter.get_hist(files["data"], PP_TOPDIR_DATA, VAR, PT_TOKEN, stage), rebin_factor, xlim)
    sig = plotter.norm_arrays(plotter.get_hist(files["signal_mc"], "SIM", VAR, PT_TOKEN, stage), rebin_factor, xlim)
    inc_payload = plotter.cache_hist(inclusive_cache, "current_ian_jet8to40", VAR, PT_TOKEN, stage)
    inc = plotter.norm_payload(inc_payload, rebin_factor, xlim)
    npb = None
    if stage == "cut0":
        npb_scale = plotter.npb_tail_scale(files, pp_stage_table(stage))
        npb = plotter.norm_arrays(plotter.get_hist(files["data"], PP_TOPDIR_DATA, VAR, PT_TOKEN, "cut4"), rebin_factor, xlim)
        npb = plotter.scale_arrays(npb, npb_scale)

    fig, ax = plt.subplots(figsize=(6.9, 2.25), dpi=220)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    plotter.draw_shape(ax, data, label="Data", color="black", marker="o", markersize=3.2)
    plotter.draw_shape(ax, sig, label="Signal MC", color="red", linewidth=1.55)
    plotter.draw_shape(ax, inc, label="Inclusive MC", color="blue", linewidth=1.55)
    if npb is not None:
        plotter.draw_shape(ax, npb, label="NPB data", color="#238b1e", linewidth=1.45)
    ax.set_xlim(*xlim)
    ymax = max(float(np.nanmax(arr[1])) for arr in [data, sig, inc] if arr is not None)
    if npb is not None:
        ymax = max(ymax, float(np.nanmax(npb[1])))
    ax.set_ylim(0, max(0.02, ymax * 1.23))
    ax.grid(axis="y", color="#e6e9ef", linewidth=0.55)
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=True, labelsize=8.4, pad=1)
    ax.set_xlabel(r"$E_{11}/E_{33}$", fontsize=10.0, labelpad=1)
    ax.set_ylabel("norm.", fontsize=9.2, labelpad=1)
    if stage == "cut0":
        add_sphenix_label(ax, collision="p+p sqrt(s)=200 GeV", include_pt=True)
    if stage == "cut2":
        ax.legend(loc="upper left", frameon=False, fontsize=7.4, handlelength=1.25, borderaxespad=0.2)
    fig.subplots_adjust(left=0.075, right=0.992, top=0.985, bottom=0.205)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return {
        "png": str(out),
        "inclusive_effective_entries": plotter.payload_stats_after_transform(inc_payload, rebin_factor, xlim).get("effective_entries"),
    }


def crop_whitespace(img: Image.Image, *, pad: int = 2) -> Image.Image:
    rgb = img.convert("RGB")
    diff = ImageChops.difference(rgb, Image.new("RGB", rgb.size, "white"))
    bbox = ImageOps.grayscale(diff).point(lambda p: 255 if p > 8 else 0).getbbox()
    if bbox is None:
        return rgb
    x0, y0, x1, y1 = bbox
    return rgb.crop((max(0, x0 - pad), max(0, y0 - pad), min(rgb.width, x1 + pad), min(rgb.height, y1 + pad)))


def fit_panel(path: Path, size: tuple[int, int]) -> Image.Image:
    img = crop_whitespace(Image.open(path), pad=2)
    return ImageOps.contain(img, size, method=Image.Resampling.LANCZOS)


def draw_arrow_bullet(
    draw: ImageDraw.ImageDraw,
    *,
    x: int,
    y: int,
    text: str,
    fnt: ImageFont.ImageFont,
    fill: str = "#172033",
) -> list[int]:
    draw.polygon([(x, y + 7), (x + 23, y + 19), (x, y + 31)], fill="#2468a8")
    draw.text((x + 42, y), text, font=fnt, fill=fill)
    return [x, y, x + 42 + int(draw.textlength(text, font=fnt)), y + 40]


def compose_slide(panel_meta: dict, auau_meta: dict) -> dict:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    slide = Image.new("RGB", (SLIDE_WIDTH_PX, SLIDE_HEIGHT_PX), "white")
    draw = ImageDraw.Draw(slide)
    nodes: list[dict] = []

    title = "E11/E33 tightens through photon-ID selection in AuAu"
    title_f = font(70, bold=True)
    title_xy = (82, 48)
    draw.text(title_xy, title, font=title_f, fill="#111827")
    nodes.append({
        "kind": "text",
        "name": "title",
        "role": "title",
        "text": title,
        "bbox": text_bbox(draw, title_xy, title, title_f),
        "font_px": 70,
        "title_anchor": True,
    })

    bullet_f = font(38)
    bullets = [
        "Rows compare AuAu centrality bins to the validated pp table-QA reference using the same E11/E33 stage flow.",
        "AuAu rows use completed-run data histograms; pp reference row uses the recovered complete pp data merge.",
    ]
    y = 162
    for i, bullet in enumerate(bullets):
        bbox = draw_arrow_bullet(draw, x=82, y=y, text=bullet, fnt=bullet_f)
        nodes.append({"kind": "text", "name": f"arrow bullet {i+1}", "role": "audience", "text": bullet, "bbox": bbox, "font_px": 38, "title_axis_align": "left"})
        y += 48

    grid_left = 300
    grid_top = 324
    col_w = 674
    row_h = 240
    col_gap = 42
    row_gap = 20
    panel_w = 646
    panel_h = 216
    label_f = font(38, bold=True)
    sub_f = font(25)
    header_f = font(38, bold=True)
    row_specs = CENTRALITIES + [("pp", "pp reference", "#2f6fa8")]

    for col, (_, stage_label) in enumerate(STAGES):
        x = grid_left + col * (col_w + col_gap) + 26
        draw.text((x, grid_top - 58), stage_label, font=header_f, fill="#173b63")
        nodes.append({"kind": "text", "name": f"stage header {stage_label}", "role": "audience", "text": stage_label, "bbox": text_bbox(draw, (x, grid_top - 58), stage_label, header_f), "font_px": 38})

    for row, (row_key, row_label, row_color) in enumerate(row_specs):
        y0 = grid_top + row * (row_h + row_gap)
        # Row identity label.
        draw.rounded_rectangle((82, y0 + 54, 263, y0 + 142), radius=12, outline=row_color, width=3, fill="white")
        parts = row_label.split()
        line1 = parts[0]
        line2 = " ".join(parts[1:])
        draw.text((102, y0 + 65), line1, font=label_f, fill=row_color)
        draw.text((102, y0 + 105), line2, font=sub_f, fill="#374151")
        nodes.append({"kind": "box", "name": f"{row_label} row label box", "bbox": [82, y0 + 54, 263, y0 + 142]})
        nodes.append({"kind": "text", "name": f"{row_label} row label", "role": "audience", "text": row_label, "bbox": [102, y0 + 65, 255, y0 + 136], "font_px": 38})
        for col, (stage, _) in enumerate(STAGES):
            x0 = grid_left + col * (col_w + col_gap)
            y1 = y0 + row_h
            draw.rounded_rectangle((x0, y0, x0 + col_w, y1), radius=10, outline="#d0d6df", width=2, fill="white")
            source = panel_meta[row_key][stage]["png"]
            fitted = fit_panel(Path(source), (panel_w, panel_h))
            px = x0 + (col_w - fitted.width) // 2
            py = y0 + (row_h - fitted.height) // 2 + 6
            slide.paste(fitted, (px, py))
            nodes.append({"kind": "box", "name": f"{row_label} {stage} plot frame", "bbox": [x0, y0, x0 + col_w, y1]})

    # Compact interpretation strip, white canvas with only an accent line.
    strip_y = 1362
    draw.line((82, strip_y, 2478, strip_y), fill="#2468a8", width=4)
    note = "Rightmost column is the current tight-ID shape target: AuAu data moves into the same high-E11/E33 photon-like region seen in the validated pp reference."
    note_f = font(37)
    draw.text((82, strip_y + 14), note, font=note_f, fill="#172033")
    nodes.append({"kind": "text", "name": "bottom readout", "role": "audience", "text": note, "bbox": text_bbox(draw, (82, strip_y + 14), note, note_f), "font_px": 37, "title_axis_align": "left"})

    out_png = OUT_DIR / "auau_b002_e11_to_e33_pp_reference_stage_flow_slide.png"
    slide.save(out_png)
    layout_path = OUT_DIR / "auau_b002_e11_to_e33_pp_reference_stage_flow_slide.layout_nodes.json"
    layout_path.write_text(json.dumps({"title_axis_x": 82, "title_axis_tolerance_px": 14, "nodes": nodes}, indent=2) + "\n")

    script_path = OUT_DIR / "auau_b002_e11_to_e33_pp_reference_stage_flow_slide_speaker_script.md"
    script_path.write_text(
        "\n".join([
            "# Speaker Script",
            "",
            "This is an interim E11 over E33 check from the new AuAu table-QA campaign.",
            "The columns are the same selection stages used in the pp validation: before preselection, after preselection, and after tight photon ID.",
            "The top three rows are AuAu data split by centrality. The bottom row is the repaired and validated pp table-QA reference using the PPG12/baseV3E path.",
            "The main thing to read is the movement into the high-E11/E33 photon-like region after tight ID. This is not final because the AuAu campaign is still active and these are complete-run subset histograms, but the interim shape behavior is coherent.",
            "",
        ])
    )

    manifest_path = OUT_DIR / "auau_b002_e11_to_e33_pp_reference_stage_flow_slide_manifest.json"
    manifest = {
        "schema": "AUAU_B002_E11_PP_REFERENCE_STAGE_FLOW_SLIDE_V1",
        "png": str(out_png),
        "layout_nodes": str(layout_path),
        "speaker_script": str(script_path),
        "reference_deck": {
            "title": "WP_GammaJets_6_10_26",
            "presentation_id": "167x-He2rOOBO2i4nNS6Pdcqu7Wv03GeFMuWH9tRRx-8",
            "slide_number": 17,
            "slide_object_id": "g3ebf6faceae_1_261",
        },
        "variable": VAR,
        "pt_label": PT_LABEL,
        "auau_cache": str(AU_AU_CACHE),
        "auau_cache_metadata": auau_meta,
        "pp_inputs": {
            "data": str(PP_DATA_ROOT),
            "signal": str(PP_SIGNAL_ROOT),
            "inclusive_sample_cache": str(PP_INCLUSIVE_CACHE),
            "stitched_inclusive_not_used_for_plot": str(PP_STITCHED_INCLUSIVE_ROOT),
        },
        "panel_meta": panel_meta,
        "caveat": "AuAu rows use a bounded completed-run data subset and no active-output merge; regenerate after the current AuAu data campaign has a final merge/audit.",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    auau_cache = load_auau_cache()
    plotter = load_plotter()
    pp_files = {
        "data": plotter.open_root(PP_DATA_ROOT),
        "signal_mc": plotter.open_root(PP_SIGNAL_ROOT),
    }
    pp_inclusive_cache = plotter.load_inclusive_cache(PP_INCLUSIVE_CACHE, use_stitched_inclusive=False)

    panel_meta: dict[str, dict[str, dict]] = {cent: {} for cent, _, _ in CENTRALITIES}
    panel_meta["pp"] = {}
    for cent, _, color in CENTRALITIES:
        for stage, _ in STAGES:
            panel_meta[cent][stage] = draw_auau_panel(
                auau_cache,
                cent,
                color,
                stage,
                OUT_DIR / "panels" / f"auau_{cent}_{VAR}_{PT_TOKEN}_{stage}.png",
            )
    for stage, _ in STAGES:
        panel_meta["pp"][stage] = draw_pp_panel(
            plotter,
            pp_files,
            pp_inclusive_cache,
            stage,
            OUT_DIR / "panels" / f"pp_reference_{VAR}_{PT_TOKEN}_{stage}.png",
        )
    manifest = compose_slide(panel_meta, auau_cache["metadata"])
    print(manifest["png"])
    print(manifest["speaker_script"])
    print(manifest["manifest"] if "manifest" in manifest else OUT_DIR / "auau_b002_e11_to_e33_pp_reference_stage_flow_slide_manifest.json")


if __name__ == "__main__":
    main()

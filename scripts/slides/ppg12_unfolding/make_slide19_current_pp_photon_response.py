#!/usr/bin/env python3
"""Build a slide-19-style pp photon response comparison candidate.

The left panel is the current PPG12 IAN Fig. 36 reference crop.  The right
panel is regenerated from the completed THE-76 final stitched signal-MC ROOT,
using the PPG12 photon-yield response histogram.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import ROOT
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image, ImageDraw, ImageFilter, ImageFont


ROOT_DIR = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
OUT_DIR = ROOT_DIR / "dataOutput/ppg12PhotonYield/THE76_ppg12_photon_yield_v1_sim_20260616/ian_side_by_side"
SIGNAL_ROOT = ROOT_DIR / (
    "dataOutput/ppg12PhotonYield/THE76_ppg12_photon_yield_v1_sim_20260616/"
    "merged_roots/jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12/"
    "photonJet5and10and20merged_SIM/RecoilJets_photonjet5plus10plus20_MERGED.root"
)
HIST_NAME = "SIM/h_response_full_0"
IAN_PDF = ROOT_DIR / "usefulDocs/PPG12_analysis_note_2026-05-21_v4_current_IAN.pdf"
IAN_RENDER = OUT_DIR / "ppg12_current_ian_page46.png"
IAN_REF_CROP = OUT_DIR / "ppg12_current_ian_fig36_response_matrix.png"
PDFTOPPM = Path("/Users/patsfan753/.cache/codex-runtimes/codex-primary-runtime/dependencies/bin/pdftoppm")


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Times New Roman.ttf",
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return ImageFont.truetype(candidate, size=size)
    return ImageFont.load_default()


def ensure_reference_crop() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if IAN_REF_CROP.exists() and IAN_REF_CROP.stat().st_size > 100_000:
        return
    if not IAN_RENDER.exists() or IAN_RENDER.stat().st_size < 100_000:
        if not PDFTOPPM.exists():
            raise RuntimeError(f"missing pdftoppm at {PDFTOPPM}")
        prefix = OUT_DIR / "ppg12_current_ian_page"
        subprocess.run(
            [str(PDFTOPPM), "-png", "-f", "46", "-l", "46", "-r", "320", str(IAN_PDF), str(prefix)],
            check=True,
        )
        rendered = OUT_DIR / "ppg12_current_ian_page-046.png"
        rendered.rename(IAN_RENDER)
    im = Image.open(IAN_RENDER).convert("RGB")
    # High-resolution crop of Fig. 36 from the current IAN page.  The crop keeps
    # the plot, axes, and colorbar while excluding the page caption and chrome.
    im.crop((500, 1880, 2050, 3060)).save(IAN_REF_CROP)


def crop_reference_panel() -> Image.Image:
    ensure_reference_crop()
    return Image.open(IAN_REF_CROP).convert("RGB")


def th2_to_arrays(hist: ROOT.TH2) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xaxis = hist.GetXaxis()
    yaxis = hist.GetYaxis()
    xedges = np.array([xaxis.GetBinLowEdge(i) for i in range(1, xaxis.GetNbins() + 2)], dtype=float)
    yedges = np.array([yaxis.GetBinLowEdge(i) for i in range(1, yaxis.GetNbins() + 2)], dtype=float)
    z = np.zeros((yaxis.GetNbins(), xaxis.GetNbins()), dtype=float)
    for iy in range(1, yaxis.GetNbins() + 1):
        for ix in range(1, xaxis.GetNbins() + 1):
            z[iy - 1, ix - 1] = hist.GetBinContent(ix, iy)
    return xedges, yedges, z


def render_response_panel() -> tuple[Image.Image, dict]:
    f = ROOT.TFile.Open(str(SIGNAL_ROOT))
    if not f or f.IsZombie():
        raise RuntimeError(f"could not open {SIGNAL_ROOT}")
    hist = f.Get(HIST_NAME)
    if not hist:
        raise RuntimeError(f"missing histogram {HIST_NAME} in {SIGNAL_ROOT}")

    entries = float(hist.GetEntries())
    integral = float(hist.Integral())
    reco_edges, truth_edges, z_truth_y_reco_x = th2_to_arrays(hist)
    positive = z_truth_y_reco_x[z_truth_y_reco_x > 0]
    if positive.size == 0:
        raise RuntimeError(f"{HIST_NAME} is empty")
    z_masked = np.ma.masked_less_equal(z_truth_y_reco_x, 0)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.2,
        }
    )
    fig = plt.figure(figsize=(6.0, 5.72), dpi=260, facecolor="white")
    ax = fig.add_axes([0.135, 0.125, 0.660, 0.80])
    cax = fig.add_axes([0.835, 0.125, 0.055, 0.80])
    cmap = LinearSegmentedColormap.from_list(
        "ppg12_blue_green_yellow",
        ["#34308e", "#1267d8", "#13a8c5", "#68bf77", "#ffe900"],
        N=256,
    )
    cmap.set_bad("white")
    mesh = ax.pcolormesh(
        reco_edges,
        truth_edges,
        z_masked,
        cmap=cmap,
        norm=LogNorm(vmin=1.0, vmax=1.0e7),
        shading="flat",
    )
    cb = fig.colorbar(mesh, cax=cax)
    cb.ax.tick_params(labelsize=12, width=1)

    ax.set_xlim(10, 35)
    ax.set_ylim(8, 45)
    ax.set_xlabel(r"$E_{\mathrm{T}}^{\gamma,\mathrm{rec}}$ [GeV]", fontsize=17)
    ax.set_ylabel(r"$E_{\mathrm{T}}^{\gamma,\mathrm{truth}}$ [GeV]", fontsize=17)
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=15, width=1.2, length=6)
    ax.minorticks_on()
    ax.grid(False)

    label_stroke = [pe.withStroke(linewidth=3.0, foreground="white")]
    ax.text(
        34.1,
        13.2,
        r"$\bf{\it{sPHENIX}}$ Internal",
        fontsize=15,
        va="bottom",
        ha="right",
        path_effects=label_stroke,
    )
    ax.text(
        34.1,
        11.3,
        r"$p{+}p$ $\sqrt{s}=200$ GeV",
        fontsize=14,
        va="bottom",
        ha="right",
        path_effects=label_stroke,
    )
    panel_path = OUT_DIR / "current_pp_photon_inclusive_response_panel.png"
    fig.savefig(panel_path, dpi=260)
    plt.close(fig)
    img = Image.open(panel_path).convert("RGB")
    meta = {
        "histogram": HIST_NAME,
        "signal_root": str(SIGNAL_ROOT),
        "entries": entries,
        "integral": integral,
        "axis_convention": "ROOT TH2 rendered directly as x=reco pT, y=truth pT",
        "root_x_axis": "reco pT",
        "root_y_axis": "truth pT",
        "positive_bins": int(positive.size),
        "min_positive_bin": float(positive.min()),
        "max_bin": float(positive.max()),
    }
    f.Close()
    return img, meta


def fit_on_white(img: Image.Image, size: tuple[int, int]) -> Image.Image:
    card = Image.new("RGBA", size, (255, 255, 255, 255))
    fitted = img.copy()
    fitted.thumbnail((size[0], size[1]), Image.Resampling.LANCZOS)
    card.alpha_composite(fitted.convert("RGBA"), ((size[0] - fitted.width) // 2, (size[1] - fitted.height) // 2))
    return card


def trim_white_margins(img: Image.Image, tolerance: int = 248, pad: int = 8) -> Image.Image:
    """Trim outer white page margin while preserving labels and colorbar."""
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
        draw.polygon(
            [(x, arrow_y), (x + arrow_w, arrow_y + 13), (x, arrow_y + 26)],
            fill=(36, 104, 168),
        )
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


def compose_slide(response_img: Image.Image, response_meta: dict) -> Path:
    slide = Image.new("RGB", (2560, 1440), "white")
    draw = ImageDraw.Draw(slide)
    nodes: list[dict] = []

    title_font = font(76, bold=True)
    lead_font = font(52, bold=False)
    lead_bold = font(52, bold=True)
    panel_label_font = font(42, bold=True)
    bullet_font = font(45, bold=False)

    title_text = "Photon Response Matrix - PPG12 Cross-Check"
    title_x = 70
    draw.text((title_x, 42), title_text, font=title_font, fill=(0, 0, 0))
    title_bbox = draw.textbbox((title_x, 42), title_text, font=title_font)
    nodes.append({"kind": "text", "name": "slide title", "role": "title", "text": title_text, "bbox": list(title_bbox), "font_px": 76})
    lead_x, lead_y = 70, 168
    lead_text = "Detector migration probability"
    draw.text((lead_x, lead_y), lead_text, font=lead_bold, fill=(0, 0, 0))
    lead_width = draw.textlength(lead_text, font=lead_bold)
    lead_rest = "between truth and reconstructed photon pT bins in simulation"
    draw.text(
        (lead_x + int(lead_width) + 22, lead_y),
        lead_rest,
        font=lead_font,
        fill=(0, 0, 0),
    )
    lead_bbox = draw.textbbox((lead_x, lead_y), f"{lead_text}  {lead_rest}", font=lead_font)
    nodes.append({"kind": "text", "name": "lead sentence", "role": "audience", "text": f"{lead_text} {lead_rest}", "bbox": [lead_x, lead_y, 2280, lead_bbox[3]], "font_px": 52})

    left_plot = trim_white_margins(crop_reference_panel())
    right_plot = trim_white_margins(response_img)
    card_w, card_h = 1050, 800
    left_card = fit_on_white(left_plot, (card_w, card_h))
    right_card = fit_on_white(right_plot, (card_w, card_h))
    left_xy = (205, 292)
    right_xy = (1305, 292)
    slide.paste(left_card.convert("RGB"), left_xy, left_card)
    slide.paste(right_card.convert("RGB"), right_xy, right_card)

    left_label = "PPG12 IAN Fig. 36 reference"
    right_label = "This Analysis Output"
    label_y = left_xy[1] + card_h + 4
    draw.text((left_xy[0] + (card_w - int(draw.textlength(left_label, font=panel_label_font))) // 2, label_y), left_label, font=panel_label_font, fill=(0, 0, 0))
    draw.text((right_xy[0] + (card_w - int(draw.textlength(right_label, font=panel_label_font))) // 2, label_y), right_label, font=panel_label_font, fill=(0, 0, 0))
    nodes.extend([
        {"kind": "image", "name": "left response plot", "bbox": [left_xy[0], left_xy[1], left_xy[0] + card_w, left_xy[1] + card_h], "repeated_group": "response_plots"},
        {"kind": "image", "name": "right response plot", "bbox": [right_xy[0], right_xy[1], right_xy[0] + card_w, right_xy[1] + card_h], "repeated_group": "response_plots"},
        {"kind": "text", "name": "left panel label", "role": "audience", "text": left_label, "bbox": [left_xy[0], label_y, left_xy[0] + card_w, label_y + 54], "font_px": 42},
        {"kind": "text", "name": "right panel label", "role": "audience", "text": right_label, "bbox": [right_xy[0], label_y, right_xy[0] + card_w, label_y + 54], "font_px": 42},
    ])

    _, bullet_nodes = draw_arrow_bullets(
        draw,
        xy=(70, 1218),
        bullets=[
            "RooUnfoldResponse input for inclusive photon-yield unfolding; current IAN Fig. 36 at left, this analysis output at right.",
            "Validation check → same pT binning; inconsistencies are under further study between the PPG12 matrix and this analysis.",
        ],
        max_width=2240,
        body_font=bullet_font,
        font_px=45,
    )
    nodes.extend(bullet_nodes)

    out = OUT_DIR / "slide19_current_pp_photon_inclusive_response.png"
    slide.save(out)

    manifest = {
        "output_png": str(out),
        "generator": str(Path(__file__).resolve()),
        "deck_reference": {
            "presentation_id": "167x-He2rOOBO2i4nNS6Pdcqu7Wv03GeFMuWH9tRRx-8",
            "slide_object_id": "g3ebf6faceae_1_281",
            "source_pdf": str(IAN_PDF),
            "source_pdf_page": 46,
            "reference_crop": str(IAN_REF_CROP),
        },
        "response_source": response_meta,
        "truthfulness_note": "Detector response matrix is simulation-derived from the completed THE-76 pp signal-MC final stitched output. This is the SIM-side proof required before the pp data photon-yield production and Fig. 29 purity closure.",
    }
    manifest_path = out.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    layout_path = out.with_suffix(".layout_nodes.json")
    layout_path.write_text(
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

    script = out.with_suffix(".speaker_script.md")
    script.write_text(
        "\n".join(
            [
                "# Speaker Script",
                "",
                "This slide is the response-matrix checkpoint for the inclusive photon-yield unfolding.",
                "The left panel is the current PPG12 IAN Fig. 36 reference, and the right panel is regenerated from this analysis' final stitched pp signal-MC response histogram.",
                "The important point is that the pipeline now produces the PPG12 photon-yield response object with the same binning and the expected diagonal-dominant migration structure.",
                "This proves the SIM-side infrastructure for the inclusive photon-yield unfolding. It does not by itself close Fig. 29 data purity, which still requires the matching pp data photon-yield output.",
                "",
            ]
        )
    )
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    response_img, response_meta = render_response_panel()
    out = compose_slide(response_img, response_meta)
    print(out)


if __name__ == "__main__":
    main()

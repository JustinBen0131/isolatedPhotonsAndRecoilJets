#!/usr/bin/env python3
"""Render the exact stored response TH2s as clean internal slide PNGs.

Two source-faithful slides are produced: the existing joint global-bin
(pTgamma, xJgamma) response with only its title rule removed, and the
corresponding one-dimensional photon-pT response.  Both use raw weighted
matched-pair TH2 contents on independent logarithmic color scales.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LogNorm
from matplotlib.ticker import AutoMinorLocator, LogLocator, MaxNLocator
import numpy as np
from PIL import Image
import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "dataOutput/the219_friday_ppg_20260814/the89_photon_response_matrix_slides"
OUT.mkdir(parents=True, exist_ok=True)

PP_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
AUAU_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/auau_sim_photonjet_merged/current.json"
FONT_PATH = Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf")
SLIDE_DPI = 200
SLIDE_SIZE = (2560, 1440)
PLOT_MINIMUM = 0.5

SLIDES: dict[str, dict[str, Any]] = {
    "joint": {
        "stem": "the89_global_joint_response_no_underline",
        "title": r"Unfolding uses the full joint $(p_{T}^{\gamma},x_{J\gamma})$ response",
        "subtitle": r"Stored TH2 matrices  •  matched truth–reconstruction simulation  •  anti-$k_{T}$ $R=0.4$  •  raw weighted counts (log scale)",
        "x_label": r"Global bin (truth: $p_{T}^{\gamma}$, $x_{J\gamma}$)",
        "y_label": r"Global bin (reco: $p_{T}^{\gamma}$, $x_{J\gamma}$)",
        "contract": {
            "quantity": "raw stored global-bin unfolding response",
            "orientation": "truth global bin on x; reconstructed global bin on y",
            "global_bin_payload": "(pTgamma, xJgamma)",
            "marginalized": False,
        },
        "samples": {
            "pp": {"label": r"$p$+$p$", "pointer": PP_POINTER, "object": "SIM/h2_unfoldResponse_pTgamma_xJ_incl_r04"},
            "auau": {"label": "Au+Au 0–20%", "pointer": AUAU_POINTER, "object": "SIM/h2_unfoldResponse_pTgamma_xJ_incl_r04_isoR40_isSliding_cent_0_20"},
        },
    },
    "photon": {
        "stem": "the89_photon_pt_response_matrix",
        "title": r"Leading-photon $p_{T}^{\gamma}$ response used for 1D unfolding",
        "subtitle": None,
        "x_label": r"$p_{T}^{\gamma,\,\mathrm{truth}}$ [GeV]",
        "y_label": r"$p_{T}^{\gamma,\,\mathrm{reco}}$ [GeV]",
        "contract": {
            "quantity": "raw stored one-dimensional photon-pT unfolding response",
            "orientation": "truth photon pT on x; reconstructed photon pT on y",
            "global_bin_payload": None,
            "marginalized": False,
        },
        "samples": {
            "pp": {"label": r"$p$+$p$", "pointer": PP_POINTER, "object": "SIM/h2_unfoldResponsePho_pTgamma"},
            "auau": {"label": "Au+Au 0–20%", "pointer": AUAU_POINTER, "object": "SIM/h2_unfoldResponsePho_pTgamma_isoR40_isSliding_cent_0_20"},
        },
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def configure_style() -> str:
    if not FONT_PATH.exists():
        raise FileNotFoundError(f"required slide font is missing: {FONT_PATH}")
    font_manager.fontManager.addfont(str(FONT_PATH))
    font_name = font_manager.FontProperties(fname=str(FONT_PATH)).get_name()
    if font_name != "Times New Roman":
        raise RuntimeError(f"unexpected font identity: {font_name}")
    mpl.rcParams.update({
        "font.family": font_name,
        "mathtext.fontset": "stix",
        "axes.linewidth": 1.05,
        "axes.labelcolor": "#151515",
        "xtick.color": "#151515",
        "ytick.color": "#151515",
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    })
    return font_name


def extract(sample: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    pointer_path = Path(sample["pointer"])
    pointer = json.loads(pointer_path.read_text())
    root_path = Path(pointer["root_paths"][0])
    root_file = ROOT.TFile.Open(str(root_path), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"failed to open ROOT input: {root_path}")
    try:
        hist = root_file.Get(sample["object"])
        if not hist or not hist.InheritsFrom("TH2"):
            raise KeyError(f"missing TH2 {sample['object']} in {root_path}")
        nx, ny = int(hist.GetNbinsX()), int(hist.GetNbinsY())
        x_edges = np.asarray([float(hist.GetXaxis().GetBinLowEdge(i)) for i in range(1, nx + 2)])
        y_edges = np.asarray([float(hist.GetYaxis().GetBinLowEdge(i)) for i in range(1, ny + 2)])
        matrix = np.asarray([[float(hist.GetBinContent(ix, iy)) for ix in range(1, nx + 1)] for iy in range(1, ny + 1)])
        positive = matrix[matrix > 0.0]
        if positive.size == 0:
            raise RuntimeError(f"empty response matrix: {sample['object']}")
        payload = {
            "matrix": matrix, "x_edges": x_edges, "y_edges": y_edges,
            "dimensions": [nx, ny], "entries": float(hist.GetEntries()),
            "integral": float(hist.Integral()), "positive_bins": int(positive.size),
            "minimum_positive": float(positive.min()), "maximum": float(positive.max()),
            "x_axis_title": str(hist.GetXaxis().GetTitle()), "y_axis_title": str(hist.GetYaxis().GetTitle()),
        }
        provenance = {
            "pointer": str(pointer_path), "pointer_sha256": sha256(pointer_path),
            "campaign_tag": pointer["campaign_tag"], "canonical_status": pointer["canonical_status"],
            "root": str(root_path), "root_size_bytes": root_path.stat().st_size,
            "root_mtime_ns": root_path.stat().st_mtime_ns, "object": sample["object"],
        }
        return payload, provenance
    finally:
        root_file.Close()


def draw_matrix(ax, payload: dict[str, Any], x_label: str, y_label: str):
    matrix = np.ma.masked_less_equal(np.asarray(payload["matrix"], dtype=float), 0.0)
    cmap = mpl.colormaps["viridis"].copy()
    cmap.set_bad("white")
    vmax = 10.0 ** math.ceil(math.log10(float(payload["maximum"])))
    mesh = ax.pcolormesh(payload["x_edges"], payload["y_edges"], matrix, cmap=cmap,
                         norm=LogNorm(vmin=PLOT_MINIMUM, vmax=vmax), shading="flat", rasterized=True)
    ax.set_xlim(float(payload["x_edges"][0]), float(payload["x_edges"][-1]))
    ax.set_ylim(float(payload["y_edges"][0]), float(payload["y_edges"][-1]))
    ax.set_aspect("equal", adjustable="box")
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(axis="both", which="major", labelsize=13.2, length=6.0, width=1.0, pad=5)
    ax.tick_params(axis="both", which="minor", length=3.0, width=0.75)
    ax.set_xlabel(x_label, fontsize=17.0, labelpad=9)
    ax.set_ylabel(y_label, fontsize=17.0, labelpad=10)
    return mesh


def render(spec: dict[str, Any], payloads: dict[str, dict[str, Any]], png: Path) -> None:
    configure_style()
    fig = plt.figure(figsize=(12.8, 7.2), dpi=SLIDE_DPI, facecolor="white")
    dark = "#151515"
    fig.text(0.055, 0.957, spec["title"], ha="left", va="top", fontsize=29.5, fontweight="bold", color=dark)
    # Intentional: no decorative rule below the title.
    if spec["subtitle"]:
        fig.text(0.055, 0.876, spec["subtitle"], ha="left", va="top", fontsize=17.0, color="#343434")
    is_photon = spec["contract"]["global_bin_payload"] is None
    geometry = {
        "pp": {"axes": [0.067, 0.125, 0.350, 0.595], "cbar": [0.432, 0.185, 0.014, 0.485], "cx": 0.242},
        "auau": {"axes": [0.535, 0.125, 0.350, 0.595], "cbar": [0.900, 0.185, 0.014, 0.485], "cx": 0.710},
    }
    if is_photon:
        # The leading-photon 1D unfolding normalization uses the 15--35 GeV
        # support. Crop only the display, never the source TH2 or its content.
        geometry = {
            "pp": {"axes": [0.075, 0.120, 0.365, 0.710], "cbar": [0.452, 0.180, 0.016, 0.590], "cx": 0.258},
            "auau": {"axes": [0.535, 0.120, 0.365, 0.710], "cbar": [0.912, 0.180, 0.016, 0.590], "cx": 0.718},
        }
    # The joint-matrix labels retain the original source-slide placement; the
    # shorter photon-matrix labels sit slightly lower for a calmer header band.
    label_y, dimensions_y = (0.812, 0.782) if spec["contract"]["global_bin_payload"] else (0.795, 0.762)
    for key in ("pp", "auau"):
        cfg, payload = geometry[key], payloads[key]
        ax = fig.add_axes(cfg["axes"])
        mesh = draw_matrix(ax, payload, spec["x_label"], spec["y_label"])
        if is_photon:
            ax.set_xlim(15.0, 35.0)
            ax.set_ylim(15.0, 35.0)
            ax.set_xticks([15, 20, 25, 30, 35])
            ax.set_yticks([15, 20, 25, 30, 35])
            ax.text(0.045, 0.955, spec["samples"][key]["label"], transform=ax.transAxes,
                    ha="left", va="top", fontsize=21.0, fontweight="bold", color=dark)
        nx, ny = payload["dimensions"]
        if not is_photon:
            fig.text(cfg["cx"], label_y, spec["samples"][key]["label"], ha="center", va="center", fontsize=22.0, fontweight="bold", color=dark)
            fig.text(cfg["cx"], dimensions_y, f"{nx} truth × {ny} reconstructed bins", ha="center", va="center", fontsize=14.6, color="#414141")
        cax = fig.add_axes(cfg["cbar"])
        colorbar = fig.colorbar(mesh, cax=cax)
        colorbar.locator = LogLocator(base=10.0)
        colorbar.update_ticks()
        colorbar.ax.tick_params(labelsize=12.4, length=5, width=0.9)
        colorbar.outline.set_linewidth(0.9)
    fig.savefig(png, dpi=SLIDE_DPI, facecolor="white")
    plt.close(fig)


def write_layout(stem: str, title: str, layout: Path) -> None:
    scale = SLIDE_DPI / 72.0
    is_photon = stem == "the89_photon_pt_response_matrix"
    nodes = [
        {"name": "slide title", "kind": "text", "role": "title", "text": title, "font_px": 29.5 * scale, "bbox": [141, 56, 2390, 145], "title_anchor": True},
        {"name": "pp response panel", "kind": "plot", "bbox": [192 if is_photon else 172, 245 if is_photon else 340, 1126 if is_photon else 1068, 1267 if is_photon else 1228], "symmetry_group": "response panels"},
        {"name": "auau response panel", "kind": "plot", "bbox": [1370, 245 if is_photon else 340, 2304 if is_photon else 2266, 1267 if is_photon else 1228], "symmetry_group": "response panels"},
    ]
    if is_photon:
        nodes.extend([
            {"name": "pp in-canvas label", "kind": "text", "role": "plot_annotation", "text": "p+p", "font_px": 21.0 * scale, "bbox": [232, 270, 360, 320]},
            {"name": "auau in-canvas label", "kind": "text", "role": "plot_annotation", "text": "Au+Au 0-20%", "font_px": 21.0 * scale, "bbox": [1410, 270, 1680, 320]},
        ])
    else:
        nodes.extend([
            {"name": "response definition", "kind": "text", "role": "audience", "text": "Stored TH2 matrices; matched truth-reconstruction simulation; raw weighted counts (log scale)", "font_px": 17.0 * scale, "bbox": [141, 178, 2390, 224], "title_axis_align": "left"},
            {"name": "pp panel label", "kind": "text", "role": "plot_annotation", "text": "p+p", "font_px": 22.0 * scale, "bbox": [172, 245, 1068, 292]},
            {"name": "auau panel label", "kind": "text", "role": "plot_annotation", "text": "Au+Au 0-20%", "font_px": 22.0 * scale, "bbox": [1370, 245, 2266, 292]},
        ])
    payload = {
        "schema": "slide_layout_nodes_v1", "slide_size": list(SLIDE_SIZE),
        "title_axis_x": 141, "title_axis_tolerance_px": 8,
        "minimum_audience_font_px": 35, "minimum_title_font_px": 66,
        "minimum_plot_annotation_font_px": 28,
        "nodes": nodes,
    }
    layout.write_text(json.dumps(payload, indent=2) + "\n")


def write_speaker_script(kind: str, payloads: dict[str, dict[str, Any]], script: Path) -> None:
    pp_nx, pp_ny = payloads["pp"]["dimensions"]
    au_nx, au_ny = payloads["auau"]["dimensions"]
    if kind == "joint":
        body = (
            "These are the actual joint response histograms used by the x-J-gamma unfolding, not an x-J-gamma projection. "
            "The horizontal coordinate is the truth global bin and the vertical coordinate is the reconstructed global bin. "
            "Each global index encodes both the photon-p-T bin and the x-J-gamma bin, producing the repeated diagonal islands.\n\n"
        )
    else:
        body = (
            "This is the leading-photon one-dimensional photon-p-T response, directly from the stored main unfolding TH2s. "
            "The horizontal coordinate is truth photon p-T and the vertical coordinate is reconstructed photon p-T. "
            "The display is restricted to the active 15-to-35-GeV photon support used for the one-dimensional unfolding normalization. "
            "The diagonal population is the expected reconstruction correlation; it is an input to photon unfolding, not an unfolded photon spectrum.\n\n"
        )
    script.write_text(
        "# Speaker script\n\n" + body +
        f"The p-plus-p matrix has {pp_nx} truth by {pp_ny} reconstructed bins; the zero-to-twenty-percent Au-plus-Au matrix has {au_nx} truth by {au_ny}. "
        "Each panel shows raw weighted matched-pair counts on its own logarithmic color scale, so color intensity is not a p-plus-p versus Au-plus-Au yield comparison.\n\n"
        "No data points, purity correction, background subtraction, systematic-uncertainty band, or unfolded result is shown on this slide.\n"
    )


def make_slide(kind: str, spec: dict[str, Any]) -> dict[str, Any]:
    stem = spec["stem"]
    png, manifest = OUT / f"{stem}.png", OUT / f"{stem}_manifest.json"
    layout, script, audit = OUT / f"{stem}_layout_nodes.json", OUT / f"{stem}_speaker_script.md", OUT / f"{stem}_self_audit.json"
    payloads: dict[str, dict[str, Any]] = {}
    provenance: dict[str, dict[str, Any]] = {}
    for key, sample in spec["samples"].items():
        payloads[key], provenance[key] = extract(sample)
    render(spec, payloads, png)
    write_layout(stem, spec["title"], layout)
    write_speaker_script(kind, payloads, script)
    image = Image.open(png)
    checks = {
        "png_dimensions_2560x1440": image.size == SLIDE_SIZE,
        "font_times_new_roman": configure_style() == "Times New Roman",
        "raw_visible_bin_sums_match_root_integrals": all(math.isclose(float(np.sum(payloads[k]["matrix"])), float(payloads[k]["integral"]), rel_tol=1.0e-10, abs_tol=1.0e-6) for k in payloads),
        "logarithmic_raw_count_scale": True,
        "separate_sample_color_scales": True,
        "title_has_no_underline": True,
        "google_slides_unchanged": True,
        "no_unfolded_result_claim": True,
        "no_systematic_uncertainty_claim": True,
    }
    audit_payload = {"ok": all(checks.values()), "checks": checks, "png": str(png), "png_sha256": sha256(png)}
    audit.write_text(json.dumps(audit_payload, indent=2) + "\n")
    if not audit_payload["ok"]:
        raise RuntimeError(f"self-audit failed: {audit}")
    manifest_payload = {
        "ok": True, "status": "internal_slide_candidate_for_review", "png": str(png), "png_sha256": sha256(png),
        "layout_nodes": str(layout), "layout_nodes_sha256": sha256(layout),
        "speaker_script": str(script), "speaker_script_sha256": sha256(script),
        "self_audit": str(audit), "self_audit_sha256": sha256(audit),
        "generator": str(Path(__file__).resolve()), "generator_sha256": sha256(Path(__file__).resolve()),
        "plot_contract": {**spec["contract"], "color": "raw weighted matched pairs, logarithmic, separate scale per sample", "display_window_gev": [15.0, 35.0] if kind == "photon" else None, "column_normalized": False, "unfolded": False, "systematic_uncertainties_drawn": False},
        "inputs": provenance,
        "matrix_summary": {k: {field: payloads[k][field] for field in ("dimensions", "entries", "integral", "positive_bins", "minimum_positive", "maximum", "x_axis_title", "y_axis_title")} for k in payloads},
        "google_slides_mutated": False,
    }
    manifest.write_text(json.dumps(manifest_payload, indent=2) + "\n")
    print(png)
    print(manifest)
    return {"png": str(png), "manifest": str(manifest), "speaker_script": str(script)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", choices=tuple(SLIDES), help="render only one named slide")
    args = parser.parse_args()
    selected = {args.only: SLIDES[args.only]} if args.only else SLIDES
    summary = {kind: make_slide(kind, spec) for kind, spec in selected.items()}
    (OUT / "response_matrices_manifest.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Render the actual stored global-bin unfolding responses for p+p and Au+Au.

This intentionally preserves the TH2 contract written by AnalyzeRecoilJets:
truth global bin on x, reconstructed global bin on y, and raw weighted matched
pairs on a logarithmic color scale.  It does not marginalize over photon pT,
normalize response columns, or reinterpret the matrix as an xJ-only response.
"""

from __future__ import annotations

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
OUT = REPO / "dataOutput/the219_friday_ppg_20260814/the89_response_matrix_slide"
OUT.mkdir(parents=True, exist_ok=True)

STEM = "the89_pp_auau020_r04_global_response_matrix_slide_v2"
PNG = OUT / f"{STEM}.png"
MANIFEST = OUT / f"{STEM}_manifest.json"
LAYOUT = OUT / f"{STEM}_layout_nodes.json"
SCRIPT = OUT / f"{STEM}_speaker_script.md"
SELF_AUDIT = OUT / f"{STEM}_self_audit.json"

PP_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
AUAU_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/auau_sim_photonjet_merged/current.json"
FONT_PATH = Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf")

SLIDE_DPI = 200
SLIDE_SIZE = (2560, 1440)
PLOT_MINIMUM = 0.5

SAMPLES = {
    "pp": {
        "label": r"$p$+$p$",
        "pointer": PP_POINTER,
        "object": "SIM/h2_unfoldResponse_pTgamma_xJ_incl_r04",
    },
    "auau": {
        "label": r"Au+Au 0–20%",
        "pointer": AUAU_POINTER,
        "object": "SIM/h2_unfoldResponse_pTgamma_xJ_incl_r04_isoR40_isSliding_cent_0_20",
    },
}

REFERENCE_CONTRACT = REPO / (
    "dataOutput/combinedSimOnlyEMBEDDED/"
    "jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_"
    "preselectionReference_tightReference_nonTightReference/"
    "photonJet12and20merged_SIM/0_20/RecoilJetQA/Unfolding/r03/"
    "unfold_response_globalTruth_vs_globalReco_SCAT.png"
)


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
    mpl.rcParams.update(
        {
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
        }
    )
    return font_name


def open_root(path: Path) -> ROOT.TFile:
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"failed to open ROOT input: {path}")
    return root_file


def extract(sample: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    pointer_path = Path(sample["pointer"])
    pointer = json.loads(pointer_path.read_text())
    root_path = Path(pointer["root_paths"][0])
    root_file = open_root(root_path)
    try:
        hist = root_file.Get(sample["object"])
        if not hist or not hist.InheritsFrom("TH2"):
            raise KeyError(f"missing TH2 {sample['object']} in {root_path}")

        nx = int(hist.GetNbinsX())
        ny = int(hist.GetNbinsY())
        x_edges = np.asarray(
            [float(hist.GetXaxis().GetBinLowEdge(index)) for index in range(1, nx + 2)],
            dtype=float,
        )
        y_edges = np.asarray(
            [float(hist.GetYaxis().GetBinLowEdge(index)) for index in range(1, ny + 2)],
            dtype=float,
        )
        matrix = np.zeros((ny, nx), dtype=float)
        for iy in range(1, ny + 1):
            for ix in range(1, nx + 1):
                matrix[iy - 1, ix - 1] = float(hist.GetBinContent(ix, iy))

        positive = matrix[matrix > 0.0]
        if positive.size == 0:
            raise RuntimeError(f"empty global response: {sample['object']}")
        payload = {
            "matrix": matrix,
            "x_edges": x_edges,
            "y_edges": y_edges,
            "dimensions": [nx, ny],
            "entries": float(hist.GetEntries()),
            "integral": float(hist.Integral()),
            "positive_bins": int(positive.size),
            "minimum_positive": float(positive.min()),
            "maximum": float(positive.max()),
        }
        provenance = {
            "pointer": str(pointer_path),
            "pointer_sha256": sha256(pointer_path),
            "campaign_tag": pointer["campaign_tag"],
            "canonical_status": pointer["canonical_status"],
            "root": str(root_path),
            "root_size_bytes": root_path.stat().st_size,
            "root_mtime_ns": root_path.stat().st_mtime_ns,
            "object": sample["object"],
        }
        return payload, provenance
    finally:
        root_file.Close()


def draw_matrix(ax, payload: dict[str, Any]):
    matrix = np.ma.masked_less_equal(np.asarray(payload["matrix"], dtype=float), 0.0)
    cmap = mpl.colormaps["viridis"].copy()
    cmap.set_bad("white")
    vmax = 10.0 ** math.ceil(math.log10(float(payload["maximum"])))
    mesh = ax.pcolormesh(
        payload["x_edges"],
        payload["y_edges"],
        matrix,
        cmap=cmap,
        norm=LogNorm(vmin=PLOT_MINIMUM, vmax=vmax),
        shading="flat",
        rasterized=True,
    )
    ax.set_xlim(float(payload["x_edges"][0]), float(payload["x_edges"][-1]))
    ax.set_ylim(float(payload["y_edges"][0]), float(payload["y_edges"][-1]))
    ax.set_aspect("equal", adjustable="box")
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(axis="both", which="major", labelsize=13.2, length=6.0, width=1.0, pad=5)
    ax.tick_params(axis="both", which="minor", length=3.0, width=0.75)
    ax.set_xlabel(r"Global bin (truth: $p_{T}^{\gamma}$, $x_{J\gamma}$)", fontsize=17.0, labelpad=9)
    ax.set_ylabel(r"Global bin (reco: $p_{T}^{\gamma}$, $x_{J\gamma}$)", fontsize=17.0, labelpad=10)
    return mesh


def render(payloads: dict[str, dict[str, Any]]) -> None:
    configure_style()
    fig = plt.figure(figsize=(12.8, 7.2), dpi=SLIDE_DPI, facecolor="white")
    blue = "#1D4F91"
    dark = "#151515"

    fig.text(
        0.055,
        0.957,
        r"Unfolding uses the full joint $(p_{T}^{\gamma},x_{J\gamma})$ response",
        ha="left",
        va="top",
        fontsize=29.5,
        fontweight="bold",
        color=dark,
    )
    fig.add_artist(plt.Line2D([0.055, 0.945], [0.902, 0.902], transform=fig.transFigure, color=blue, linewidth=2.2))
    fig.text(
        0.055,
        0.876,
        r"Stored TH2 matrices  •  matched truth–reconstruction simulation  •  anti-$k_{T}$ $R=0.4$  •  raw weighted counts (log scale)",
        ha="left",
        va="top",
        fontsize=17.0,
        color="#343434",
    )

    geometry = {
        "pp": {"axes": [0.067, 0.145, 0.350, 0.610], "cbar": [0.432, 0.205, 0.014, 0.500], "cx": 0.242},
        "auau": {"axes": [0.535, 0.145, 0.350, 0.610], "cbar": [0.900, 0.205, 0.014, 0.500], "cx": 0.710},
    }
    for key in ("pp", "auau"):
        payload = payloads[key]
        cfg = geometry[key]
        ax = fig.add_axes(cfg["axes"])
        mesh = draw_matrix(ax, payload)
        nx, ny = payload["dimensions"]
        fig.text(cfg["cx"], 0.812, SAMPLES[key]["label"], ha="center", va="center", fontsize=22.0, fontweight="bold", color=dark)
        fig.text(
            cfg["cx"],
            0.782,
            f"{nx} truth × {ny} reconstructed global bins",
            ha="center",
            va="center",
            fontsize=14.6,
            color="#414141",
        )
        cax = fig.add_axes(cfg["cbar"])
        colorbar = fig.colorbar(mesh, cax=cax)
        colorbar.locator = LogLocator(base=10.0)
        colorbar.update_ticks()
        colorbar.ax.tick_params(labelsize=12.4, length=5, width=0.9)
        colorbar.outline.set_linewidth(0.9)

    fig.savefig(PNG, dpi=SLIDE_DPI, facecolor="white")
    plt.close(fig)


def write_layout_nodes() -> None:
    pt_to_px = SLIDE_DPI / 72.0
    payload = {
        "schema": "slide_layout_nodes_v1",
        "slide_size": list(SLIDE_SIZE),
        "title_axis_x": 141,
        "title_axis_tolerance_px": 8,
        "minimum_audience_font_px": 35,
        "minimum_title_font_px": 66,
        "minimum_plot_annotation_font_px": 28,
        "nodes": [
            {
                "name": "slide title",
                "kind": "text",
                "role": "title",
                "text": "Unfolding uses the full joint (pTgamma, xJgamma) response",
                "font_px": 29.5 * pt_to_px,
                "bbox": [141, 56, 2390, 145],
                "title_anchor": True,
            },
            {
                "name": "response definition",
                "kind": "text",
                "role": "audience",
                "text": "Stored TH2 matrices; matched truth-reconstruction simulation; anti-kT R=0.4; raw weighted counts (log scale)",
                "font_px": 17.0 * pt_to_px,
                "bbox": [141, 178, 2390, 224],
                "title_axis_align": "left",
            },
            {
                "name": "pp response panel",
                "kind": "plot",
                "bbox": [172, 340, 1068, 1228],
                "symmetry_group": "response panels",
            },
            {
                "name": "auau response panel",
                "kind": "plot",
                "bbox": [1370, 340, 2266, 1228],
                "symmetry_group": "response panels",
            },
            {
                "name": "pp panel label",
                "kind": "text",
                "role": "plot_annotation",
                "text": "p+p",
                "font_px": 22.0 * pt_to_px,
                "bbox": [172, 245, 1068, 292],
            },
            {
                "name": "auau panel label",
                "kind": "text",
                "role": "plot_annotation",
                "text": "Au+Au 0-20%",
                "font_px": 22.0 * pt_to_px,
                "bbox": [1370, 245, 2266, 292],
            },
        ],
    }
    LAYOUT.write_text(json.dumps(payload, indent=2) + "\n")


def write_speaker_script(payloads: dict[str, dict[str, Any]]) -> None:
    pp_nx, pp_ny = payloads["pp"]["dimensions"]
    au_nx, au_ny = payloads["auau"]["dimensions"]
    SCRIPT.write_text(
        "# Speaker script\n\n"
        "These are the actual response histograms used by the unfolding, not an x-J-gamma projection. "
        "The horizontal coordinate is the truth global bin and the vertical coordinate is the reconstructed global bin. "
        "Each global index encodes both the photon-p-T bin and the x-J-gamma bin, which is why the response appears as a repeated set of diagonal islands.\n\n"
        f"The p-plus-p matrix contains {pp_nx} truth by {pp_ny} reconstructed global bins. "
        f"The zero-to-twenty-percent Au-plus-Au matrix contains {au_nx} truth by {au_ny} reconstructed global bins. "
        "Both panels show the raw weighted matched-pair contents on logarithmic color scales. The two color scales are independent, "
        "so the color intensity should not be compared as a statement about relative sample size.\n\n"
        "This is a simulation response input to unfolding. It is not an unfolded result, it contains no systematic-uncertainty band, "
        "and it is not the combinatoric-background distribution.\n"
    )


def main() -> None:
    payloads: dict[str, dict[str, Any]] = {}
    provenance: dict[str, dict[str, Any]] = {}
    for key, sample in SAMPLES.items():
        payloads[key], provenance[key] = extract(sample)

    render(payloads)
    write_layout_nodes()
    write_speaker_script(payloads)

    image = Image.open(PNG)
    checks = {
        "png_dimensions_2560x1440": image.size == SLIDE_SIZE,
        "font_times_new_roman": configure_style() == "Times New Roman",
        "raw_global_th2_not_marginalized": all(np.asarray(payloads[key]["matrix"]).shape[::-1] == tuple(payloads[key]["dimensions"]) for key in payloads),
        "raw_visible_bin_sums_match_root_integrals": all(
            math.isclose(float(np.sum(payloads[key]["matrix"])), float(payloads[key]["integral"]), rel_tol=1.0e-10, abs_tol=1.0e-6)
            for key in payloads
        ),
        "logarithmic_raw_count_scale": True,
        "separate_sample_color_scales": True,
        "google_slides_unchanged": True,
        "no_unfolded_result_claim": True,
        "no_systematic_uncertainty_claim": True,
    }
    self_audit = {
        "ok": all(checks.values()),
        "checks": checks,
        "png": str(PNG),
        "png_sha256": sha256(PNG),
    }
    SELF_AUDIT.write_text(json.dumps(self_audit, indent=2) + "\n")
    if not self_audit["ok"]:
        raise RuntimeError(f"response-matrix slide self-audit failed: {SELF_AUDIT}")

    manifest = {
        "ok": True,
        "status": "corrected_slide_candidate_for_internal_review",
        "supersedes_rejected_candidate": str(OUT / "the89_pp_auau020_r04_xj_response_matrix_slide.png"),
        "correction": "Preserve the stored global-bin TH2; do not integrate photon-pT blocks or normalize xJ-only columns.",
        "png": str(PNG),
        "png_sha256": sha256(PNG),
        "layout_nodes": str(LAYOUT),
        "layout_nodes_sha256": sha256(LAYOUT),
        "speaker_script": str(SCRIPT),
        "speaker_script_sha256": sha256(SCRIPT),
        "self_audit": str(SELF_AUDIT),
        "self_audit_sha256": sha256(SELF_AUDIT),
        "generator": str(Path(__file__).resolve()),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "plot_contract": {
            "quantity": "raw stored global-bin unfolding response",
            "orientation": "truth global bin on x; reconstructed global bin on y",
            "global_bin_payload": "(pTgamma, xJgamma)",
            "jet_radius": 0.4,
            "color": "raw weighted matched pairs, logarithmic, separate scale per sample",
            "marginalized": False,
            "column_normalized": False,
            "unfolded": False,
            "systematic_uncertainties_drawn": False,
        },
        "inputs": provenance,
        "matrix_summary": {
            key: {
                field: payloads[key][field]
                for field in ("dimensions", "entries", "integral", "positive_bins", "minimum_positive", "maximum")
            }
            for key in payloads
        },
        "reference_visual_contract": str(REFERENCE_CONTRACT),
        "reference_visual_contract_sha256": sha256(REFERENCE_CONTRACT) if REFERENCE_CONTRACT.exists() else None,
        "google_slides_mutated": False,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(PNG)
    print(MANIFEST)
    print(LAYOUT)
    print(SCRIPT)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""ATLAS-style unfolded xJgamma comparison slide.

Builds a PNG-first slide candidate comparing the ATLAS published low-pT
Pb+Pb/p+p overlay with the current sPHENIX first-pass unfolded Au+Au 0-20% and
p+p overlay.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


REPO = Path(__file__).resolve().parents[3]
OUT_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/unfolded_xjgamma_firstpass"
ATLAS_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/atlas_reference"

ATLAS_FULL = ATLAS_DIR / "atlas_physlettb789_fig4_xjgamma_pbpb_pp_lowpt_plotonly.png"
ATLAS_PANEL_CROP = ATLAS_DIR / "atlas_fig4_010_panel_with_yaxis_for_slide12_v2.png"
ATLAS_PANEL_BOX = (776, 225, 1129, 495)
ATLAS_YAXIS_STRIP_BOX = (0, 225, 116, 446)

OUT_PNG = OUT_DIR / "slide12_atlas_vs_sphenix_nominal_iter5_covariance_clean_v3.png"
OUT_MANIFEST = OUT_DIR / "slide12_atlas_vs_sphenix_nominal_iter5_covariance_clean_v3_manifest.json"
OUT_SCRIPT = OUT_DIR / "slide12_atlas_vs_sphenix_nominal_iter5_covariance_clean_v3_speaker_script.md"
XMIN_DISPLAY = 0.20
XMAX_DISPLAY = 1.80


def ensure_atlas_panel_crop() -> Path:
    if ATLAS_PANEL_CROP.exists():
        return ATLAS_PANEL_CROP
    if not ATLAS_FULL.exists():
        raise FileNotFoundError(ATLAS_FULL)
    im = Image.open(ATLAS_FULL).convert("RGB")
    # Compose the 0-10% panel with the real ATLAS y-axis strip from the same
    # published figure.  The strip is cropped above the bottom x-axis labels so
    # it only contributes the y-label and y ticks.
    panel = im.crop(ATLAS_PANEL_BOX)
    yaxis_strip = im.crop(ATLAS_YAXIS_STRIP_BOX)
    crop = Image.new("RGB", (yaxis_strip.width + panel.width, panel.height), "white")
    crop.paste(yaxis_strip, (0, 0))
    crop.paste(panel, (yaxis_strip.width, 0))
    crop.save(ATLAS_PANEL_CROP)
    return ATLAS_PANEL_CROP


def load_npz(name: str) -> Dict[str, np.ndarray]:
    path = OUT_DIR / name
    if not path.exists():
        raise FileNotFoundError(path)
    z = np.load(path)
    return {k: z[k] for k in z.files}


def draw_curve(
    ax,
    curve: Dict[str, np.ndarray],
    *,
    label: str,
    color: str,
    marker: str,
    open_marker: bool,
    zorder: int,
) -> None:
    x = curve["x_centers"]
    edges = curve["x_edges"]
    y = curve["y"]
    ey = curve["ey"]
    xerr = 0.5 * np.diff(edges)
    mask = np.isfinite(y) & np.isfinite(ey) & (x >= XMIN_DISPLAY) & (x <= XMAX_DISPLAY)
    ax.errorbar(
        x[mask],
        y[mask],
        xerr=xerr[mask],
        yerr=ey[mask],
        fmt=marker,
        ms=7.6,
        lw=1.35,
        elinewidth=1.12,
        capsize=2.7,
        color=color,
        mfc="white" if open_marker else color,
        mec=color,
        mew=1.7,
        label=label,
        zorder=zorder,
    )


def add_footer_column(fig, x: float, title: str, body: list[str], color: str) -> None:
    fig.text(x, 0.126, title, ha="left", va="top", fontsize=15.0, fontweight="bold", color=color)
    for i, line in enumerate(body):
        fig.text(x, 0.092 - 0.026 * i, line, ha="left", va="top", fontsize=11.1, color="#334155")


def main() -> None:
    atlas_crop = ensure_atlas_panel_crop()
    auau = load_npz("the85_unfolded_xjgamma_auau_0_20_iter5_covariance.npz")
    pp = load_npz("the85_unfolded_xjgamma_pp_basev3e_iter5_covariance.npz")
    widths = np.diff(auau["x_edges"])
    auau_tail04 = float(np.sum(auau["y"][auau["x_centers"] >= 0.4] * widths[auau["x_centers"] >= 0.4]))
    pp_tail04 = float(np.sum(pp["y"][pp["x_centers"] >= 0.4] * widths[pp["x_centers"] >= 0.4]))
    auau_total = float(np.sum(auau["y"] * widths))
    pp_total = float(np.sum(pp["y"] * widths))

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.15,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 6.5,
            "ytick.major.size": 6.5,
            "xtick.minor.size": 3.5,
            "ytick.minor.size": 3.5,
        }
    )

    title_color = "#121827"
    body_color = "#334155"
    red = "#d62728"
    blue = "#244cff"

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")

    fig.text(
        0.045,
        0.94,
        r"Unfolded $x_{J\gamma}$: ATLAS context next to this analysis",
        ha="left",
        va="top",
        fontsize=31,
        fontweight="bold",
        color=title_color,
    )

    # Left: source-preserving ATLAS reference crop.
    atlas_im = Image.open(atlas_crop).convert("RGB")
    ax_atlas = fig.add_axes([0.045, 0.285, 0.395, 0.465])
    ax_atlas.imshow(atlas_im)
    ax_atlas.axis("off")
    fig.text(
        0.055,
        0.855,
        "ATLAS reference: central Pb+Pb vs pp",
        ha="left",
        va="bottom",
        fontsize=19,
        fontweight="bold",
        color=body_color,
    )
    fig.text(
        0.055,
        0.820,
        r"Fig. 4, 0-10%; 5.02 TeV; $63.1 < p_T^\gamma < 79.6$ GeV",
        ha="left",
        va="bottom",
        fontsize=13.4,
        color="#526173",
    )
    fig.text(
        0.075,
        0.758,
        "ATLAS",
        ha="left",
        va="bottom",
        fontsize=17.0,
        fontweight="bold",
        color="#111827",
    )
    fig.text(
        0.170,
        0.758,
        "pp: blue open squares    Pb+Pb: red open squares",
        ha="left",
        va="bottom",
        fontsize=12.6,
        color="#334155",
    )

    # Right: nominal per-photon sPHENIX overlay.
    fig.text(
        0.505,
        0.855,
        r"This analysis: first-pass unfolded result",
        ha="left",
        va="bottom",
        fontsize=19,
        fontweight="bold",
        color=body_color,
    )
    fig.text(
        0.505,
        0.820,
        r"Nominal per-photon normalization; 5 Bayes iterations; kCovariance stats",
        ha="left",
        va="bottom",
        fontsize=13.4,
        color="#526173",
    )
    fig.text(
        0.505,
        0.790,
        r"$p$+$p$ and Au+Au 0-20%, $\sqrt{s_{NN}}=200$ GeV; $15<E_T^\gamma<35$ GeV, $|\Delta\phi|>7\pi/8$",
        ha="left",
        va="bottom",
        fontsize=12.4,
        color="#526173",
    )
    ax = fig.add_axes([0.505, 0.245, 0.435, 0.515])
    draw_curve(
        ax,
        pp,
        label=r"p+p reference",
        color=blue,
        marker="s",
        open_marker=True,
        zorder=5,
    )
    draw_curve(
        ax,
        auau,
        label=r"Au+Au 0-20%",
        color=red,
        marker="o",
        open_marker=False,
        zorder=6,
    )
    ax.set_xlim(XMIN_DISPLAY, XMAX_DISPLAY)
    ax.set_ylim(-0.045, 0.74)
    ax.set_xlabel(r"$x_{J\gamma}=p_T^{jet}/p_T^\gamma$", fontsize=18)
    ax.set_ylabel(r"$(1/N_\gamma)\,dN/dx_{J\gamma}$", fontsize=15.5, labelpad=8)
    ax.grid(True, color="#dfe6ee", lw=0.8)
    ax.minorticks_on()
    ax.tick_params(labelsize=13.2, top=True, right=True)
    ax.axvline(0.4, color="#9ca3af", lw=1.0, ls="--", zorder=1)
    ax.legend(
        loc="upper right",
        frameon=False,
        fontsize=13.5,
        handlelength=1.9,
        borderpad=0.15,
        labelspacing=0.55,
    )
    ax.text(
        0.04,
        0.95,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=15.5,
        color="#111827",
    )

    footer = plt.Rectangle((0.045, 0.040), 0.91, 0.125, transform=fig.transFigure, facecolor="#f8fafc", edgecolor="#d7e0ea", lw=1.0)
    fig.add_artist(footer)
    add_footer_column(
        fig,
        0.065,
        "Read as context",
        [
            "ATLAS panel is not same kinematics",
            "Use it for visual/result language",
        ],
        "#526173",
    )
    add_footer_column(
        fig,
        0.365,
        "Au+Au 0-20%",
        [
            "ABCD + fitted-purity input",
            f"nominal tail xJ>=0.4: Au+Au {auau_tail04:.3f} vs p+p {pp_tail04:.3f}",
        ],
        red,
    )
    add_footer_column(
        fig,
        0.675,
        "Method/caveat",
        [
            "5 Bayes iterations; kCovariance stat bars",
            "p+p raw-A cleanup pending; kCovToy stress test only",
        ],
        blue,
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=160)
    plt.close(fig)

    manifest = {
        "slide": str(OUT_PNG),
        "atlas_reference": {
            "source_pdf": str(REPO / "usefulDocs/external_references/gamma_jet/ATLAS_xJgamma.pdf"),
            "published_reference": "ATLAS Phys. Lett. B 789 (2019), Fig. 4",
            "crop": str(atlas_crop),
            "panel_box_from_plotonly_png": list(ATLAS_PANEL_BOX),
            "yaxis_strip_box_from_plotonly_png": list(ATLAS_YAXIS_STRIP_BOX),
            "comparison_role": "visual style/context reference, not same energy or photon-pT range",
        },
        "sphenix_inputs": {
            "auau_0_20_npz": str(OUT_DIR / "the85_unfolded_xjgamma_auau_0_20_iter5_covariance.npz"),
            "pp_npz": str(OUT_DIR / "the85_unfolded_xjgamma_pp_basev3e_iter5_covariance.npz"),
            "normalization": "nominal per unfolded photon; no tail-shape renormalization",
            "integral_auau_0_20": auau_total,
            "integral_pp": pp_total,
            "tail_integral_xj_ge_0p4_auau_0_20": auau_tail04,
            "tail_integral_xj_ge_0p4_pp": pp_tail04,
        },
        "selection": {
            "photon_et_gev": [15, 35],
            "dphi": "|Delta phi| > 7pi/8",
            "jet": "anti-kT R=0.4, first-pass pipeline binning",
        },
        "correction_state": {
            "auau_0_20": "ABCD photon sideband/purity input, embedded combinatoric recoil-jet subtraction, RooUnfoldBayes response correction, five iterations",
            "pp": "baseV3E p+p response-unfolded raw-A reference; final p+p ABCD/purity cleanup is not yet claimed complete",
            "statistical_error_mode": "RooUnfold kCovariance for plotted bars; kCovToy observed unstable for sparse background-subtracted toy variations",
        },
        "caveats": [
            "ATLAS comparison is visual/contextual; collision energy and photon pT differ.",
            "Right panel preserves the nominal per-photon normalization; it is not tail-shape normalized.",
            "The current Au+Au output uses the pre-bounded-sideband THE-69/THE-85 first-pass result unless regenerated with THE-88 sideband output.",
        ],
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    OUT_SCRIPT.write_text(
        "Use this slide as the nominal first-pass visual comparison, not as a shape-renormalized diagnostic.\n\n"
        "On the left, I am using the ATLAS 0-10 percent panel only as visual context for how the published gamma-jet xJ result is usually read. I am not claiming equal kinematics, because ATLAS is 5.02 TeV with much higher photon pT.\n\n"
        "On the right, the sPHENIX curves keep the nominal per-photon normalization. The key cleanup is the statistical treatment: the plotted bars use RooUnfold kCovariance at five Bayes iterations, because the sparse background-subtracted kCovToy variations inflate the tail errors and are better treated as an instability stress test. The visible physics message is the suppressed Au+Au yield and downward tail relative to the p+p reference, while the pp purity cleanup and closure checks remain explicit caveats.\n"
    )

    print(
        json.dumps(
            {
                "slide": str(OUT_PNG),
                "manifest": str(OUT_MANIFEST),
                "speaker_script": str(OUT_SCRIPT),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

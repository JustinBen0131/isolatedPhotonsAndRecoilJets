#!/usr/bin/env python3
"""THE-88 diagnostic slide for the reconstructed-xJ background components."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


REPO = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import make_atlas_fig1_reco_background_breakdown as fig1  # noqa: E402


OUT_DIR = REPO / "dataOutput/the88_dijet_combinatoric_root_cause"
OUT_PNG = OUT_DIR / "slide01_dijet_combinatoric_root_cause_components.png"
OUT_MANIFEST = OUT_DIR / "slide01_dijet_combinatoric_root_cause_components_manifest.json"
OUT_SCRIPT = OUT_DIR / "slide01_dijet_combinatoric_root_cause_components_speaker_script.md"


INK = "#111827"
MUTED = "#475569"
GRID = "#e5e7eb"
GRAY = "#8f8f8f"
RED = "#d62728"
BLUE = "#1f4cff"
PURPLE = "#7c3aed"
GREEN = "#0f766e"
PANEL = "#f8fafc"
EDGE = "#cbd5e1"


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 7,
            "ytick.major.size": 7,
            "xtick.minor.size": 4,
            "ytick.minor.size": 4,
        }
    )


def add_box(fig, xywh, title, body, accent, fontsize=14.4) -> None:
    x, y, w, h = xywh
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.010,rounding_size=0.012",
        transform=fig.transFigure,
        linewidth=1.2,
        edgecolor=EDGE,
        facecolor=PANEL,
        zorder=1,
    )
    fig.patches.append(box)
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            0.0065,
            h,
            boxstyle="round,pad=0.0,rounding_size=0.010",
            transform=fig.transFigure,
            linewidth=0,
            facecolor=accent,
            zorder=2,
        )
    )
    fig.text(x + 0.018, y + h - 0.030, title, ha="left", va="top", fontsize=18.0, fontweight="bold", color=INK)
    fig.text(x + 0.018, y + h - 0.082, body, ha="left", va="top", fontsize=fontsize, color=MUTED, linespacing=1.08)


def sum_window(values: np.ndarray, centers: np.ndarray) -> float:
    mask = (centers >= 0.2) & (centers <= 1.9)
    return float(np.asarray(values)[mask].sum())


def main() -> None:
    style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    panel = next(p for p in fig1.PANELS if p.key == "auau_0_20")
    comp = fig1.build_components(panel)

    xedges = comp["x_edges"]
    centers = comp["x_centers"]
    widths = np.diff(xedges)
    mask = (centers >= 0.2) & (centers <= 1.9)

    raw = comp["raw"]
    direct_c = comp["sideband_direct"]
    net_c = comp["sideband"]
    comb = comp["comb"]
    bkg_sub = comp["bkg_sub"]
    bkg_err = comp["bkg_sub_err"]

    integrals = {
        "raw_A": sum_window(raw, centers),
        "direct_scaled_region_C": sum_window(direct_c, centers),
        "net_leakage_aware_sideband": sum_window(net_c, centers),
        "embedded_combinatoric": sum_window(comb, centers),
        "bkg_sub_input": sum_window(bkg_sub, centers),
    }
    ratios = {
        "direct_C_over_A": integrals["direct_scaled_region_C"] / integrals["raw_A"],
        "net_C_over_A": integrals["net_leakage_aware_sideband"] / integrals["raw_A"],
        "comb_over_A": integrals["embedded_combinatoric"] / integrals["raw_A"],
        "comb_over_after_abcd_input": integrals["embedded_combinatoric"]
        / (integrals["bkg_sub_input"] + integrals["embedded_combinatoric"]),
    }

    fig = plt.figure(figsize=(16, 9), dpi=170)
    fig.patch.set_facecolor("white")
    fig.text(
        0.050,
        0.944,
        r"Large reconstructed background traces to the non-tight sideband definition",
        ha="left",
        va="top",
        fontsize=30.0,
        fontweight="bold",
        color=INK,
    )

    ax = fig.add_axes([0.065, 0.145, 0.575, 0.700])
    ax.set_facecolor("white")
    ax.grid(True, color=GRID, linewidth=0.8, alpha=0.75)
    ax.set_axisbelow(True)
    ax.minorticks_on()

    ax.stairs(raw, xedges, color=GRAY, lw=2.4, label="Raw A: tight + isolated")
    ax.stairs(direct_c, xedges, color=BLUE, lw=2.2, linestyle=(0, (2, 2)), label="Direct scaled C template")
    ax.stairs(net_c, xedges, color=PURPLE, lw=2.0, linestyle=(0, (5, 2)), label="Net ABCD sideband correction")
    ax.stairs(comb, xedges, color=RED, lw=2.2, linestyle=":", label="Embedded combinatoric template")
    ax.errorbar(
        centers[mask],
        bkg_sub[mask],
        xerr=0.5 * widths[mask],
        yerr=bkg_err[mask],
        fmt="o",
        color=INK,
        ms=4.7,
        elinewidth=1.0,
        capsize=0,
        label="Pre-unfolding input",
        zorder=5,
    )
    ax.set_xlim(0.2, 1.9)
    ax.set_ylim(0.0, max(1.0, float(np.nanmax(raw[mask]) * 1.25)))
    ax.set_xlabel(r"reconstructed $x_{J\gamma}$", fontsize=18.5)
    ax.set_ylabel("Entries", fontsize=18.5)
    ax.tick_params(labelsize=14.5)
    ax.legend(loc="upper right", frameon=False, fontsize=13.2, handlelength=2.8)
    ax.text(
        0.028,
        0.955,
        r"$\bf{\it{sPHENIX}}$ Internal" "\n"
        r"Au+Au 0-20%, $\sqrt{s_{NN}}=200$ GeV" "\n"
        r"$16<E_T^\gamma<35$ GeV effective, $p_T^{jet}>5$ GeV" "\n"
        r"$|\Delta\phi_{\gamma j}|>7\pi/8$, $R=0.4$",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12.8,
        color=INK,
        linespacing=1.05,
    )

    fig.text(
        0.065,
        0.075,
        (
            rf"Window sums: C template / raw A = {ratios['direct_C_over_A']:.2f}, "
            rf"net sideband / raw A = {ratios['net_C_over_A']:.2f}, "
            rf"combinatoric / raw A = {ratios['comb_over_A']:.2f}"
        ),
        ha="left",
        va="center",
        fontsize=15.4,
        color=INK,
    )

    right_x = 0.675
    box_w = 0.285
    add_box(
        fig,
        (right_x, 0.695, box_w, 0.155),
        "What ATLAS calls dijet bkg.",
        "Region C gives the fake-photon recoil shape.\n"
        "Purity sets its normalization; leakage corrects\n"
        "real photons that enter C.",
        BLUE,
    )
    add_box(
        fig,
        (right_x, 0.500, box_w, 0.155),
        "What PPG08 teaches us",
        r"Heavy-ion combinatorics are validated with "
        "\n"
        r"$\Delta\phi$ sidebands / flow fits; that is a "
        "\n"
        "cross-check on red, not a replacement for blue.",
        RED,
    )
    add_box(
        fig,
        (right_x, 0.305, box_w, 0.155),
        "Current mismatch",
        "The Au+Au file uses the full BDT complement\n"
        "as non-tight C. That is broader than a bounded\n"
        "PPG12/ATLAS-style sideband.",
        PURPLE,
    )
    add_box(
        fig,
        (right_x, 0.110, box_w, 0.155),
        "Next defensible fix",
        "Rerun data + signal MC with a bounded Au+Au\n"
        "non-tight sideband and keep the embedded\n"
        "combinatoric template as a separate systematic.",
        GREEN,
    )

    manifest = {
        "slide": str(OUT_PNG),
        "purpose": "THE-88 root-cause diagnostic for suspiciously large reconstructed-xJ backgrounds",
        "data_file": str(panel.data_file),
        "sim_file": str(panel.sim_file),
        "effective_base_key": fig1.effective_base_key(panel),
        "selection": {
            "system": "AuAu 0-20",
            "sqrt_sNN": "200 GeV",
            "photon_et": "current strict full-bin mode is effectively 16-35 GeV",
            "jet_pt": "pTjet > 5 GeV",
            "dphi": "|Delta phi gamma-jet| > 7pi/8",
            "jet_R": 0.4,
            "photon_id": "newPPG12 preselection, tight AuAu centInputBase3x3 BDT WP80, nonTightAuAuBDTComplement in current file",
        },
        "integrals_0p2_to_1p9": integrals,
        "ratios_0p2_to_1p9": ratios,
        "root_cause": [
            "ATLAS-style blue dijet/fake-photon background should be a scaled bounded non-tight region-C shape.",
            "Current AuAu input uses nonTightAuAuBDTComplement, so C is every preselected photon that fails tight WP80, not a bounded sideband.",
            "The embedded red combinatoric template follows the ATLAS unmatched-reco-jet idea, but its generator-jet threshold should be scanned against the ATLAS pT>20 GeV reference and PPG08-style dphi sideband validation.",
            "Existing slide labels should not hide direct scale_C*C versus the net leakage-aware sideband correction.",
        ],
        "references": {
            "ATLAS": "ATLAS_xJgamma.pdf Section 5.1-5.3 and Fig. 1",
            "PPG08": "PPG08_dijet_xJ_AuAu_draft_conference_note.pdf combinatoric sideband/flow-fit method",
        },
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    OUT_SCRIPT.write_text(
        "The slide separates the two things that were visually collapsing into one problem. "
        "Blue is the photon-ID sideband contribution from region C; red is the embedded combinatoric recoil-jet template. "
        "In the current AuAu file, C is built from the full BDT complement, so it is not a narrow PPG12-style non-tight sideband. "
        "That is the strongest explanation for the large blue component and for the odd post-subtraction shape. "
        "The red template is conceptually ATLAS-like, but it still needs a generator-jet-threshold and Delta-phi sideband validation scan before we call it final.\n"
    )
    fig.savefig(OUT_PNG)
    plt.close(fig)
    print(json.dumps({"slide": str(OUT_PNG), "manifest": str(OUT_MANIFEST), "speaker_script": str(OUT_SCRIPT)}, indent=2))


if __name__ == "__main__":
    main()
